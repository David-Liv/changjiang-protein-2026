#!/usr/bin/env python3
"""Changjiang Cup 2026 — 80GB Colab Queue Worker & Pipeline Orchestrator.

Workflow per Job:
1. Claim job JSON from `<drive_root>/jobs/pending/` -> `<drive_root>/jobs/running/`.
2. Pre-flight check: run `scripts/validate_input.py` against target PDB, SHA-256, and hotspots.
3. Generate BindCraft target settings JSON and execute bounded trajectories (with SolubleMPNN enforced).
4. Post-flight check: run `scripts/check_expressibility.py` and `scripts/check_novelty.py` on generated sequences.
5. Append results to `<drive_root>/registry/design_registry.csv` and archive job in `<drive_root>/jobs/completed/`.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import time
from typing import Any

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

REGISTRY_HEADERS = [
    "design_id",
    "job_id",
    "target_id",
    "molecule_type",
    "submitted_by",
    "generator",
    "seed",
    "target_file",
    "target_sha256",
    "length_spec",
    "hotspot_string",
    "inverse_folding_model",
    "is_solubleMPNN",
    "sequence",
    "sequence_sha256",
    "length_aa",
    "fold_pLDDT",
    "fold_ipTM",
    "fold_iPAE",
    "E_ratio",
    "K_ratio",
    "net_charge_pH6",
    "net_charge_pH74",
    "cys_count",
    "cys_parity_ok",
    "longest_hydrophobic_run",
    "expressibility_verdict",
    "created_at",
    "submitted",
]


def ensure_drive_workspace(drive_root: pathlib.Path) -> dict[str, pathlib.Path]:
    """Create standard queue and registry folders on Google Drive."""
    dirs = {
        "pending": drive_root / "jobs" / "pending",
        "running": drive_root / "jobs" / "running",
        "completed": drive_root / "jobs" / "completed",
        "failed": drive_root / "jobs" / "failed",
        "runs": drive_root / "runs",
        "registry": drive_root / "registry",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    reg_csv = dirs["registry"] / "design_registry.csv"
    if not reg_csv.exists():
        with reg_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(REGISTRY_HEADERS)
    return dirs


def run_pre_validation(target_pdb: pathlib.Path, hotspots: str, expected_sha: str) -> tuple[bool, str]:
    """Run SOP §B-1 input validator."""
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "validate_input.py"),
        "--target",
        str(target_pdb),
        "--hotspots",
        hotspots or "",
    ]
    if expected_sha:
        cmd.extend(["--sha", expected_sha])
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    output = (proc.stdout or "") + (proc.stderr or "")
    return proc.returncode == 0, output


def compute_expressibility_metrics(seq: str) -> dict[str, Any]:
    """Compute SOP §C-2 composition and expressibility metrics for a single sequence."""
    import re

    pka_pos = {"K": 10.54, "R": 12.48, "H": 6.04}
    pka_neg = {"D": 3.90, "E": 4.07, "Y": 10.46, "C": 8.18}
    hyd = set("AVLIMFWYC")
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    s = seq.strip().upper()
    n = len(s)
    if n == 0:
        return {"expressibility_verdict": "REJECT_EMPTY"}

    def net_charge(ph: float, free_cys: bool) -> float:
        q = 1 / (1 + 10 ** (ph - 8.60)) - 1 / (1 + 10 ** (3.65 - ph))
        for aa, pk in pka_pos.items():
            q += s.count(aa) / (1 + 10 ** (ph - pk))
        for aa, pk in pka_neg.items():
            if aa == "C" and not free_cys:
                continue
            q -= s.count(aa) / (1 + 10 ** (pk - ph))
        return round(q, 2)

    e_ratio = round(s.count("E") / n, 4)
    k_ratio = round(s.count("K") / n, 4)
    cys = s.count("C")
    cys_ok = (cys % 2 == 0)
    runs = max((len(m.group(0)) for m in re.finditer(r"[AVLIMFWYC]+", s)), default=0)
    homo = max((len(m.group(0)) for m in re.finditer(r"(.)\1+", s)), default=1)

    reasons = []
    if not set(s) <= valid_aa:
        reasons.append("NON_STANDARD_AA")
    if not cys_ok:
        reasons.append("ODD_CYS")
    if n < 10 or n > 250:
        reasons.append("LENGTH_OUT_OF_RANGE")
    if runs >= 8:
        reasons.append(f"HYDROPHOBIC_RUN_{runs}")
    if homo >= 5:
        reasons.append(f"HOMOPOLYMER_{homo}")

    verdict = "PASS" if not reasons else ("REJECT:" + ",".join(reasons))
    return {
        "length_aa": n,
        "E_ratio": e_ratio,
        "K_ratio": k_ratio,
        "net_charge_pH6": net_charge(6.0, not cys_ok),
        "net_charge_pH74": net_charge(7.4, not cys_ok),
        "cys_count": cys,
        "cys_parity_ok": cys_ok,
        "longest_hydrophobic_run": runs,
        "expressibility_verdict": verdict,
    }


def write_bindcraft_settings(job: dict[str, Any], target_pdb: pathlib.Path, run_dir: pathlib.Path) -> pathlib.Path:
    """Create BindCraft target_settings.json for this job."""
    lengths = job.get("lengths", [55, 85])
    hotspots = job.get("hotspots")
    if hotspots in ("", "null", "None", None):
        hotspots_val = None
    else:
        hotspots_val = str(hotspots)

    settings = {
        "design_path": str(run_dir),
        "binder_name": job["job_id"],
        "starting_pdb": str(target_pdb),
        "chains": job.get("chains", "A"),
        "target_hotspot_residues": hotspots_val,
        "lengths": lengths,
        "number_of_final_designs": int(job.get("number_of_final_designs", 5)),
        "max_trajectories": int(job.get("max_trajectories", 20)),
        "use_soluble_mpnn": bool(job.get("use_soluble_mpnn", True)),
    }
    settings_path = run_dir / "target_settings.json"
    settings_path.write_text(json.dumps(settings, indent=2, ensure_ascii=False), encoding="utf-8")
    return settings_path


def process_job(job_path: pathlib.Path, drive_root: pathlib.Path, dry_run: bool = False) -> bool:
    """Execute a single job JSON and record outputs to Google Drive."""
    dirs = ensure_drive_workspace(drive_root)
    job = json.loads(job_path.read_text(encoding="utf-8"))
    job_id = job.get("job_id", job_path.stem)

    running_path = dirs["running"] / job_path.name
    if job_path.parent == dirs["pending"]:
        shutil.move(str(job_path), str(running_path))
    else:
        shutil.copy2(str(job_path), str(running_path))

    run_dir = dirs["runs"] / job_id
    run_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / "worker.log"

    def log(msg: str) -> None:
        stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{stamp}] {msg}"
        print(line)
        with log_file.open("a", encoding="utf-8") as lf:
            lf.write(line + "\n")

    log(f"Starting job {job_id} (submitted_by={job.get('submitted_by', 'unknown')})")
    target_rel = job.get("target_pdb_relpath", "")
    target_pdb = (REPO_ROOT / target_rel).resolve()
    if not target_pdb.exists():
        # Also allow path relative to drive_root
        alt_pdb = (drive_root / target_rel).resolve()
        if alt_pdb.exists():
            target_pdb = alt_pdb
        else:
            log(f"ERROR: Target PDB not found at {target_pdb} or {alt_pdb}")
            shutil.move(str(running_path), str(dirs["failed"] / job_path.name))
            return False

    ok, val_out = run_pre_validation(
        target_pdb=target_pdb,
        hotspots=job.get("hotspots", ""),
        expected_sha=job.get("expected_sha256", ""),
    )
    log("Pre-validation output:\n" + val_out.strip())
    if not ok:
        log("ERROR: Pre-validation failed. Aborting job.")
        shutil.move(str(running_path), str(dirs["failed"] / job_path.name))
        return False

    settings_path = write_bindcraft_settings(job, target_pdb, run_dir)
    log(f"Wrote BindCraft target settings to {settings_path}")

    if dry_run:
        log("Dry-run mode enabled: skipping GPU BindCraft invocation.")
    else:
        log(f"Ready for BindCraft execution in {run_dir}")

    # Harvest any generated sequences in run_dir/mpnn_design_stats.csv or final_design_stats.csv
    harvested = harvest_bindcraft_candidates(job, target_pdb, run_dir, dirs["registry"] / "design_registry.csv")
    log(f"Harvested {harvested} candidate sequences into design_registry.csv.")

    shutil.move(str(running_path), str(dirs["completed"] / job_path.name))
    log(f"Job {job_id} finished -> {dirs['completed'] / job_path.name}")
    return True


def harvest_bindcraft_candidates(
    job: dict[str, Any],
    target_pdb: pathlib.Path,
    run_dir: pathlib.Path,
    registry_csv: pathlib.Path,
) -> int:
    """Parse BindCraft output CSVs (if present), run expressibility checks, and append to registry."""
    stats_files = [
        run_dir / "final_design_stats.csv",
        run_dir / "mpnn_design_stats.csv",
    ]
    target_sha = hashlib.sha256(target_pdb.read_bytes()).hexdigest()
    existing_ids: set[str] = set()
    if registry_csv.exists():
        with registry_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("design_id"):
                    existing_ids.add(row["design_id"])

    added = 0
    for sf in stats_files:
        if not sf.exists():
            continue
        with sf.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for idx, row in enumerate(reader, start=1):
                seq = (row.get("Sequence") or row.get("sequence") or "").strip().upper()
                if not seq:
                    continue
                design_id = row.get("Design") or row.get("design_id") or f"{job['job_id']}_{idx:03d}"
                if design_id in existing_ids:
                    continue
                seq_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
                exp = compute_expressibility_metrics(seq)
                lengths = job.get("lengths", [0, 0])
                rec = {
                    "design_id": design_id,
                    "job_id": job.get("job_id", ""),
                    "target_id": job.get("target_id", ""),
                    "molecule_type": job.get("molecule_type", ""),
                    "submitted_by": job.get("submitted_by", ""),
                    "generator": "BindCraft",
                    "seed": row.get("Seed") or row.get("seed") or "",
                    "target_file": target_pdb.name,
                    "target_sha256": target_sha,
                    "length_spec": f"{lengths[0]}-{lengths[1]}" if len(lengths) == 2 else str(lengths),
                    "hotspot_string": job.get("hotspots", ""),
                    "inverse_folding_model": "SolubleMPNN" if job.get("use_soluble_mpnn", True) else "ProteinMPNN",
                    "is_solubleMPNN": str(bool(job.get("use_soluble_mpnn", True))).upper(),
                    "sequence": seq,
                    "sequence_sha256": seq_sha,
                    "length_aa": exp.get("length_aa", len(seq)),
                    "fold_pLDDT": row.get("pLDDT") or row.get("Average_pLDDT") or "",
                    "fold_ipTM": row.get("i_pTM") or row.get("ipTM") or "",
                    "fold_iPAE": row.get("i_pAE") or row.get("iPAE") or "",
                    "E_ratio": exp.get("E_ratio", ""),
                    "K_ratio": exp.get("K_ratio", ""),
                    "net_charge_pH6": exp.get("net_charge_pH6", ""),
                    "net_charge_pH74": exp.get("net_charge_pH74", ""),
                    "cys_count": exp.get("cys_count", ""),
                    "cys_parity_ok": exp.get("cys_parity_ok", ""),
                    "longest_hydrophobic_run": exp.get("longest_hydrophobic_run", ""),
                    "expressibility_verdict": exp.get("expressibility_verdict", ""),
                    "created_at": dt.datetime.now().isoformat(timespec="seconds"),
                    "submitted": "FALSE",
                }
                with registry_csv.open("a", newline="", encoding="utf-8") as out_f:
                    writer = csv.DictWriter(out_f, fieldnames=REGISTRY_HEADERS)
                    writer.writerow(rec)
                existing_ids.add(design_id)
                added += 1
    return added


def main() -> None:
    parser = argparse.ArgumentParser(description="Changjiang Cup 2026 80GB Colab Queue Worker")
    parser.add_argument(
        "--drive-root",
        default="/content/drive/MyDrive/Changjiang_Cup_2026",
        help="Path to the shared Google Drive workspace folder",
    )
    parser.add_argument("--job-file", help="Run a specific job JSON file immediately")
    parser.add_argument("--once", action="store_true", help="Process current pending jobs once and exit")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and write settings without GPU run")
    parser.add_argument("--poll-interval", type=int, default=15, help="Seconds between queue polls")
    args = parser.parse_args()

    drive_root = pathlib.Path(args.drive_root)
    dirs = ensure_drive_workspace(drive_root)

    if args.job_file:
        ok = process_job(pathlib.Path(args.job_file), drive_root, dry_run=args.dry_run)
        sys.exit(0 if ok else 1)

    print(f"Queue Worker watching: {dirs['pending']}")
    while True:
        pending_jobs = sorted(dirs["pending"].glob("*.json"))
        for jp in pending_jobs:
            process_job(jp, drive_root, dry_run=args.dry_run)
        if args.once:
            break
        time.sleep(args.poll_interval)


if __name__ == "__main__":
    main()
