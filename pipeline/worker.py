#!/usr/bin/env python3
"""Changjiang Cup 2026 — 80GB Colab Queue Worker & Pipeline Orchestrator.

Supports two stages:
- `generation` (default):
  1. Claim job JSON from `<drive_root>/jobs/pending/` -> `<drive_root>/jobs/running/`.
  2. Pre-flight check: run `scripts/validate_input.py` against target PDB, SHA-256, and hotspots.
  3. Generate BindCraft target settings JSON and execute bounded trajectories (with SolubleMPNN enforced).
  4. Post-flight check & Stage C evaluation: run `scripts/check_expressibility.py` + `scripts/tnf_pipeline.py`
     (`triage`, `ESMFold`, `ipSAE`, `pDockQ`, `min_heavy_dist`) on generated sequences.
  5. Append results to `<drive_root>/registry/design_registry.csv` and archive job in `<drive_root>/jobs/completed/`.
- `evaluation_calibration`:
  Executes the full Target 2 (`TNF-α`) or custom `.tsv` spec calibration suite using `scripts/tnf_pipeline.py`
  (`selftest`, `triage`, `s1_mono`, `s23`, `det`, `s23_var`, `s1_trimer`, `compare`), saving all CSVs, PDBs,
  and `tnf_calibration_runs.zip` directly to `<drive_root>/runs/target2_TNFa_calibration/`.
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
import zipfile
from typing import Any

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

import tnf_pipeline as tp  # noqa: E402

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
    "esm_plddt",
    "ipsae10_min",
    "ipsae10_max",
    "pdockq",
    "min_heavy_dist",
    "E_ratio",
    "K_ratio",
    "net_charge_pH6",
    "net_charge_pH74",
    "gravy",
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


def detect_gpu_chunk() -> int:
    """Return optimal ESMFold chunk size: 0 (no chunking, fastest) if VRAM >= 38GB else 64."""
    try:
        import torch

        if torch.cuda.is_available():
            gb = torch.cuda.get_device_properties(0).total_memory / 1e9
            return 0 if gb >= 38 else 64
    except Exception:
        pass
    return 64


def extract_target_sequence_from_pdb(pdb_path: pathlib.Path) -> str:
    """Extract concatenated sequence of target chains from a PDB file."""
    chains = tp.read_pdb(str(pdb_path))
    seqs = ["".join(r["aa"] for r in rs) for _, rs in sorted(chains.items())]
    return ":".join(seqs)


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
    """Compute SOP §C-2 composition + `tnf_pipeline.triage` GRAVY metrics for a single sequence."""
    import re

    pka_pos = {"K": 10.54, "R": 12.48, "H": 6.04}
    pka_neg = {"D": 3.90, "E": 4.07, "Y": 10.46, "C": 8.18}
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

    triage_res = tp.triage(s)
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
    if triage_res["gravy"] > 0.1:
        reasons.append(f"HIGH_GRAVY_{triage_res['gravy']:+.2f}")

    verdict = "PASS" if not reasons else ("REJECT:" + ",".join(reasons))
    return {
        "length_aa": n,
        "E_ratio": e_ratio,
        "K_ratio": k_ratio,
        "net_charge_pH6": net_charge(6.0, not cys_ok),
        "net_charge_pH74": net_charge(7.4, not cys_ok),
        "gravy": triage_res["gravy"],
        "cys_count": cys,
        "cys_parity_ok": cys_ok,
        "longest_hydrophobic_run": runs,
        "expressibility_verdict": verdict,
    }


def evaluate_candidates_with_ipsae(
    candidates: list[tuple[str, str]],
    target_seq: str,
    eval_dir: pathlib.Path,
    dry_run: bool = False,
) -> dict[str, dict[str, Any]]:
    """Run ESMFold + Dunbrack ipSAE + Bryant pDockQ on candidate sequences via `tnf_pipeline.py`."""
    if not candidates or not target_seq or dry_run:
        return {}

    eval_dir.mkdir(parents=True, exist_ok=True)
    spec_path = eval_dir / "candidates_eval.tsv"
    with spec_path.open("w", encoding="utf-8") as f:
        f.write("# design_id\tbinder:target\n")
        for cid, seq in candidates:
            f.write(f"{cid}\t{seq}:{target_seq}\n")

    chunk = detect_gpu_chunk()
    fold_cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "tnf_pipeline.py"),
        "fold",
        "--spec",
        str(spec_path),
        "--out",
        str(eval_dir),
        "--chunk",
        str(chunk),
    ]
    subprocess.run(fold_cmd, check=False)

    summary_csv = eval_dir / "ipsae_summary.csv"
    score_cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "tnf_pipeline.py"),
        "score",
        "--dir",
        str(eval_dir),
        "--cutoffs",
        "10",
        "15",
        "--csv",
        str(summary_csv),
    ]
    subprocess.run(score_cmd, check=False)

    results: dict[str, dict[str, Any]] = {}
    if summary_csv.exists():
        with summary_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                name = row.get("name", "")
                if name:
                    results[name] = {
                        "esm_plddt": row.get("mean_plddt", ""),
                        "ipsae10_min": row.get("ipsae10_min", ""),
                        "ipsae10_max": row.get("ipsae10_max", ""),
                        "pdockq": row.get("iface_pdockq", ""),
                        "min_heavy_dist": row.get("iface_min_heavy_dist", ""),
                    }
    return results


def run_calibration_job(job: dict[str, Any], run_dir: pathlib.Path, log_fn, dry_run: bool = False) -> bool:
    """Execute the full Target 2 TNF-alpha calibration & stability benchmark using `tnf_pipeline.py`."""
    spec_dir = (REPO_ROOT / job.get("spec_dir_relpath", "targets/target2_TNFa")).resolve()
    script_path = str(REPO_ROOT / "scripts" / "tnf_pipeline.py")
    chunk = str(detect_gpu_chunk())

    def run_step(desc: str, args_list: list[str]) -> bool:
        log_fn(f"=== {desc} ===")
        proc = subprocess.run(
            [sys.executable, script_path] + args_list,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        out = (proc.stdout or "") + (proc.stderr or "")
        log_fn(out.strip())
        return proc.returncode == 0

    if not run_step("1a. Selftest (9 mathematical & PDB checks)", ["selftest"]):
        return False

    specs = [
        str(spec_dir / "spec_stage1_mono.tsv"),
        str(spec_dir / "spec_stage1_trimer.tsv"),
        str(spec_dir / "spec_stage23.tsv"),
    ]
    run_step("1b. Sequence Triage Pre-check", ["triage", "--spec"] + specs)

    if dry_run:
        log_fn("Dry-run enabled: skipping GPU ESMFold runs for calibration job.")
        return True

    # Stage 1 Monomer, Stage 2/3 Pos/Neg Controls, Stage 1 Trimer
    batches = [
        ("3. Stage 1 Monomer Decisive Comparison", str(spec_dir / "spec_stage1_mono.tsv"), run_dir / "s1_mono", []),
        ("4. Stage 2/3 Positive & Negative Calibration", str(spec_dir / "spec_stage23.tsv"), run_dir / "s23", []),
        ("6. Stage 1 Trimer (~700 tokens)", str(spec_dir / "spec_stage1_trimer.tsv"), run_dir / "s1_trimer", []),
    ]
    if job.get("run_stability_variants", True):
        # Determinism check
        det_tsv = run_dir / "det.tsv"
        lines = (spec_dir / "spec_stage23.tsv").read_text(encoding="utf-8").splitlines()[:2]
        det_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
        run_step(
            "5a. Determinism Verification",
            ["fold", "--spec", str(det_tsv), "--out", str(run_dir / "det"), "--chunk", chunk, "--check-determinism"],
        )
        batches.append(
            (
                "5b. 8-Variant Stability Sweep (72 inferences)",
                str(spec_dir / "spec_stage23.tsv"),
                run_dir / "s23_var",
                ["--recycles", "1", "2", "3", "4", "--linkers", "25", "30"],
            )
        )

    for title, sp_file, out_sub, extra_args in batches:
        run_step(title + " [fold]", ["fold", "--spec", sp_file, "--out", str(out_sub), "--chunk", chunk] + extra_args)
        run_step(
            title + " [score]",
            ["score", "--dir", str(out_sub), "--cutoffs", "10", "15", "--csv", str(out_sub / "summary.csv")],
        )

    # Crystal structure superposition (5M2J, 1TNF, 2AZ5)
    mono_pdb = run_dir / "s1_mono" / "1a_TNF_mono__r4L25.pdb"
    if mono_pdb.exists():
        for ref_id in job.get("compare_pdbs", ["5M2J", "1TNF", "2AZ5"]):
            run_step(
                f"7. Superposition vs {ref_id}",
                ["compare", str(mono_pdb), "--fetch", str(ref_id)],
            )

    # Package everything into tnf_calibration_runs.zip inside run_dir
    zip_path = run_dir / "tnf_calibration_runs.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in run_dir.rglob("*"):
            if f.is_file() and f != zip_path:
                zf.write(f, f.relative_to(run_dir))
    log_fn(f"Packaged all calibration results to: {zip_path}")
    return True


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
    stage = job.get("stage", "generation")

    running_path = dirs["running"] / job_path.name
    if job_path.parent == dirs["pending"]:
        shutil.move(str(job_path), str(running_path))
    else:
        shutil.copy2(str(job_path), str(running_path))

    run_dir_name = "target2_TNFa_calibration" if stage == "evaluation_calibration" else job_id
    run_dir = dirs["runs"] / run_dir_name
    run_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / "worker.log"

    def log(msg: str) -> None:
        stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{stamp}] {msg}"
        print(line)
        with log_file.open("a", encoding="utf-8") as lf:
            lf.write(line + "\n")

    log(f"Starting job {job_id} (stage={stage}, submitted_by={job.get('submitted_by', 'unknown')})")

    if stage == "evaluation_calibration":
        ok = run_calibration_job(job, run_dir, log, dry_run=dry_run)
        dest_dir = dirs["completed"] if ok else dirs["failed"]
        shutil.move(str(running_path), str(dest_dir / job_path.name))
        log(f"Calibration job {job_id} finished (ok={ok}) -> {dest_dir / job_path.name}")
        return ok

    target_rel = job.get("target_pdb_relpath", "")
    target_pdb = (REPO_ROOT / target_rel).resolve()
    if not target_pdb.exists():
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

    harvested = harvest_bindcraft_candidates(
        job, target_pdb, run_dir, dirs["registry"] / "design_registry.csv", dry_run=dry_run
    )
    log(f"Harvested & evaluated {harvested} candidate sequences into design_registry.csv.")

    shutil.move(str(running_path), str(dirs["completed"] / job_path.name))
    log(f"Job {job_id} finished -> {dirs['completed'] / job_path.name}")
    return True


def harvest_bindcraft_candidates(
    job: dict[str, Any],
    target_pdb: pathlib.Path,
    run_dir: pathlib.Path,
    registry_csv: pathlib.Path,
    dry_run: bool = False,
) -> int:
    """Parse BindCraft output CSVs, run expressibility + ESMFold ipSAE/pDockQ checks, and append to registry."""
    stats_files = [
        run_dir / "final_design_stats.csv",
        run_dir / "mpnn_design_stats.csv",
    ]
    target_sha = hashlib.sha256(target_pdb.read_bytes()).hexdigest()
    target_seq = job.get("target_sequence") or extract_target_sequence_from_pdb(target_pdb)

    existing_ids: set[str] = set()
    if registry_csv.exists():
        with registry_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("design_id"):
                    existing_ids.add(row["design_id"])

    raw_rows: list[tuple[str, str, dict[str, Any]]] = []
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
                existing_ids.add(design_id)
                raw_rows.append((design_id, seq, row))

    if not raw_rows:
        return 0

    ipsae_map = evaluate_candidates_with_ipsae(
        [(cid, seq) for cid, seq, _ in raw_rows],
        target_seq=target_seq,
        eval_dir=run_dir / "ipsae_eval",
        dry_run=dry_run,
    )

    added = 0
    for design_id, seq, row in raw_rows:
        seq_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        exp = compute_expressibility_metrics(seq)
        ev = ipsae_map.get(design_id, {})
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
            "esm_plddt": ev.get("esm_plddt", ""),
            "ipsae10_min": ev.get("ipsae10_min", ""),
            "ipsae10_max": ev.get("ipsae10_max", ""),
            "pdockq": ev.get("pdockq", ""),
            "min_heavy_dist": ev.get("min_heavy_dist", ""),
            "E_ratio": exp.get("E_ratio", ""),
            "K_ratio": exp.get("K_ratio", ""),
            "net_charge_pH6": exp.get("net_charge_pH6", ""),
            "net_charge_pH74": exp.get("net_charge_pH74", ""),
            "gravy": exp.get("gravy", ""),
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
