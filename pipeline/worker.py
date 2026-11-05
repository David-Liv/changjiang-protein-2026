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
import os
import pathlib
import shutil
import subprocess
import sys
import time
import traceback
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


def inspect_gpu_environment() -> dict[str, Any]:
    """Inspect active GPU hardware, compute capability, and estimated Colab CU burn rate (1,000 CU budget)."""
    info: dict[str, Any] = {
        "cuda_available": False,
        "gpu_name": "CPU",
        "vram_gb": 0.0,
        "sm_arch": "cpu",
        "major": 0,
        "minor": 0,
        "est_cu_per_hour": 0.1,
        "bindcraft_compatible": True,
        "warning": "",
    }
    try:
        import torch

        if torch.cuda.is_available():
            p = torch.cuda.get_device_properties(0)
            gb = round(p.total_memory / 1e9, 1)
            major, minor = int(p.major), int(p.minor)
            sm = f"sm_{major}{minor}"
            name_up = p.name.upper()
            if "H100" in name_up or major == 9:
                cu_rate = 18.0
                warn = "H100 burns ~18 CU/hr (~55 hrs on 1,000 CU). A100 is recommended for better CU efficiency."
                bc_compat = True
            elif major >= 12 or "BLACKWELL" in name_up or "RTX PRO 6000" in name_up or "G4" in name_up:
                cu_rate = 9.5
                warn = (
                    f"Detected Blackwell {p.name} ({sm}). BindCraft's CUDA 12.4 / JAX stack does NOT support {sm}. "
                    "Please switch Colab runtime to A100 for BindCraft generation jobs!"
                )
                bc_compat = False
            elif "A100" in name_up or (major == 8 and minor == 0):
                cu_rate = 7.5 if gb >= 70 else 5.4
                warn = ""
                bc_compat = True
            else:
                cu_rate = 4.0
                warn = ""
                bc_compat = True

            info.update(
                {
                    "cuda_available": True,
                    "gpu_name": p.name,
                    "vram_gb": gb,
                    "sm_arch": sm,
                    "major": major,
                    "minor": minor,
                    "est_cu_per_hour": cu_rate,
                    "bindcraft_compatible": bc_compat,
                    "warning": warn,
                }
            )
    except Exception:
        pass
    return info


def detect_gpu_chunk() -> int:
    """Return optimal ESMFold chunk size: 0 (no chunking, fastest) if VRAM >= 38GB else 64."""
    info = inspect_gpu_environment()
    return 0 if info["vram_gb"] >= 38 else 64


def sync_filesystem() -> None:
    """Flush Python streams and OS filesystem buffers to ensure Google Drive FUSE writes complete."""
    sys.stdout.flush()
    sys.stderr.flush()
    if hasattr(os, "sync"):
        try:
            os.sync()
        except Exception:
            pass


def flush_drive_and_unassign_gpu(drive_root: pathlib.Path, unassign: bool = True, countdown_sec: int = 15) -> None:
    """Verify results are saved on Google Drive, flush Drive mount, and unassign Colab GPU runtime."""
    sync_filesystem()
    if not drive_root.exists():
        print(
            f"[SAFEGUARD] Drive root {drive_root} does not exist! "
            "Skipping GPU unassign so local files are not lost."
        )
        return

    saved_files = [p for p in (drive_root / "runs").rglob("*") if p.is_file()] if (drive_root / "runs").exists() else []
    print(f"\n=== [Auto-Save Verification] {len(saved_files)} files persisted under {drive_root / 'runs'} ===")
    if len(saved_files) == 0:
        print("[SAFEGUARD] No output files found in runs/ — skipping GPU unassign!")
        return

    if unassign and countdown_sec > 0:
        print(
            f"\n⏳ 任務已全部完成並寫入 Google Drive！將在 {countdown_sec} 秒後自動釋放 GPU 以節省 1,000 CU 額度..."
            "\n   （若想留在 Notebook 查看結果，請在倒數結束前按上方『停止 ■』鍵即可保留連線）",
            flush=True,
        )
        try:
            for remaining in range(countdown_sec, 0, -5):
                print(f"   倒數 {remaining} 秒...", flush=True)
                time.sleep(min(5, remaining))
        except KeyboardInterrupt:
            print("\n✋ 已取消自動斷開 GPU，當前 Colab 連線保留。")
            return

    try:
        from google.colab import drive  # type: ignore

        print("Flushing Google Drive FUSE buffers (drive.flush_and_unmount)...", flush=True)
        drive.flush_and_unmount()
        print("Google Drive safely flushed and unmounted ✓", flush=True)
    except Exception as e:
        print(f"[Info] Google Drive unmount skipped (not in Colab or already unmounted: {e})", flush=True)

    if unassign:
        try:
            from google.colab import runtime  # type: ignore

            print("Releasing Colab GPU runtime now (runtime.unassign())...", flush=True)
            time.sleep(1)
            runtime.unassign()
        except Exception as e:
            print(f"[Info] runtime.unassign() skipped (not running inside Google Colab: {e})", flush=True)


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


def check_calibration_completeness(
    run_dir: pathlib.Path, require_variants: bool = True
) -> tuple[bool, list[str], list[str]]:
    """Verify whether all stages of Target 2 TNF-alpha calibration are truly completed on Google Drive."""
    stages = [
        ("s1_mono (單體 3 條)", run_dir / "s1_mono", 3),
        ("s23 (陽/陰性校準 8 條)", run_dir / "s23", 8),
        ("s1_trimer (三聚體 2 條)", run_dir / "s1_trimer", 2),
    ]
    if require_variants:
        stages.append(("s23_var (穩定度 8×8=64 變體)", run_dir / "s23_var", 64))

    done_list: list[str] = []
    missing_list: list[str] = []
    for label, sdir, expected_n in stages:
        pdbs = list(sdir.glob("*.pdb")) if sdir.exists() else []
        csv_ok = (sdir / "summary.csv").exists() and (sdir / "summary.csv").stat().st_size > 0
        if len(pdbs) >= expected_n and csv_ok:
            done_list.append(f"{label}: ✅ {len(pdbs)}/{expected_n} PDB + summary.csv")
        else:
            missing_list.append(f"{label}: ⏳ {len(pdbs)}/{expected_n} PDB (summary.csv={'有' if csv_ok else '無'})")

    zip_path = run_dir / "tnf_calibration_runs.zip"
    if zip_path.exists() and zip_path.stat().st_size > 0 and not missing_list:
        done_list.append(f"tnf_calibration_runs.zip: ✅ ({zip_path.stat().st_size // 1024} KB)")
    else:
        missing_list.append("tnf_calibration_runs.zip: ⏳ 待打包")

    return len(missing_list) == 0, done_list, missing_list


def print_drive_diagnostic(drive_root: pathlib.Path) -> None:
    """Print detailed diagnostic of Google Drive queue and Target 2 calibration progress."""
    dirs = ensure_drive_workspace(drive_root)
    print("\n================ [Google Drive 狀態與進度診斷] ================")
    for k in ("running", "pending", "completed", "failed"):
        files = [p.name for p in sorted(dirs[k].glob("*.json"))]
        print(f"  📂 jobs/{k:<10}: {len(files)} 個 -> {files}")

    cal_dir = dirs["runs"] / "target2_TNFa_calibration"
    if cal_dir.exists():
        is_done, done_items, missing_items = check_calibration_completeness(cal_dir)
        print(f"\n  📊 Target 2 (TNF-α) 校準進度 ({'🎉 100% 全部完成' if is_done else '🔄 部分完成，可斷點續跑'}):")
        for item in done_items:
            print(f"     {item}")
        for item in missing_items:
            print(f"     {item}")
        manifest_p = cal_dir / "run_manifest.json"
        if manifest_p.exists():
            try:
                m = json.loads(manifest_p.read_text(encoding="utf-8"))
                print(
                    f"  🕒 上次運行記錄: status={m.get('status')}, 耗時={m.get('elapsed_seconds')}s, "
                    f"預估消耗={m.get('estimated_cu_used')} CU, GPU={m.get('gpu_name')}"
                )
            except Exception:
                pass
        log_p = cal_dir / "worker.log"
        if log_p.exists():
            lines = log_p.read_text(encoding="utf-8", errors="replace").splitlines()
            print("\n  📜 worker.log 最後 12 行記錄:")
            for ln in lines[-12:]:
                print(f"     {ln}")
    print("===============================================================\n", flush=True)


def run_calibration_job(job: dict[str, Any], run_dir: pathlib.Path, log_fn, dry_run: bool = False) -> bool:
    """Execute the full Target 2 TNF-alpha calibration & stability benchmark using `tnf_pipeline.py`."""
    spec_dir = (REPO_ROOT / job.get("spec_dir_relpath", "targets/target2_TNFa")).resolve()
    script_path = str(REPO_ROOT / "scripts" / "tnf_pipeline.py")
    chunk = str(detect_gpu_chunk())

    def run_step(desc: str, args_list: list[str]) -> bool:
        log_fn(f"=== {desc} ===")
        proc = subprocess.Popen(
            [sys.executable, "-u", script_path] + args_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            clean = line.rstrip("\r\n")
            if clean:
                log_fn(clean)
        proc.wait()
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

    all_steps_ok = True
    for title, sp_file, out_sub, extra_args in batches:
        fold_ok = run_step(
            title + " [fold]", ["fold", "--spec", sp_file, "--out", str(out_sub), "--chunk", chunk] + extra_args
        )
        score_ok = run_step(
            title + " [score]",
            ["score", "--dir", str(out_sub), "--cutoffs", "10", "15", "--csv", str(out_sub / "summary.csv")],
        )
        all_steps_ok = all_steps_ok and fold_ok and score_ok

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

    is_done, _, missing = check_calibration_completeness(
        run_dir, require_variants=bool(job.get("run_stability_variants", True))
    )
    if not is_done:
        log_fn(f"WARNING: Calibration job incomplete, missing items: {missing}")
    return all_steps_ok and is_done


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
    """Execute a single job JSON and record outputs to Google Drive with guaranteed persistence."""
    t_start = time.time()
    dirs = ensure_drive_workspace(drive_root)
    job = json.loads(job_path.read_text(encoding="utf-8"))
    job_id = job.get("job_id", job_path.stem)
    stage = job.get("stage", "generation")

    running_path = dirs["running"] / job_path.name
    is_resumed = job_path.resolve() == running_path.resolve()
    if job_path.parent == dirs["pending"]:
        shutil.move(str(job_path), str(running_path))
    elif not is_resumed:
        shutil.copy2(str(job_path), str(running_path))

    run_dir_name = "target2_TNFa_calibration" if stage == "evaluation_calibration" else job_id
    run_dir = dirs["runs"] / run_dir_name
    run_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / "worker.log"

    def log(msg: str) -> None:
        stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{stamp}] {msg}"
        print(line, flush=True)
        with log_file.open("a", encoding="utf-8") as lf:
            lf.write(line + "\n")

    gpu_info = inspect_gpu_environment()
    action_str = "Resuming interrupted job" if is_resumed else "Starting job"
    log(
        f"{action_str} {job_id} (stage={stage}, submitted_by={job.get('submitted_by', 'unknown')}) | "
        f"GPU={gpu_info['gpu_name']} ({gpu_info['vram_gb']} GB, {gpu_info['sm_arch']}, "
        f"~{gpu_info['est_cu_per_hour']} CU/hr on 1,000 CU Ultra 5x quota)"
    )
    if gpu_info["warning"]:
        log(f"[GPU WARNING] {gpu_info['warning']}")

    ok = False
    interrupted = False
    try:
        if stage == "generation" and not dry_run and not gpu_info["bindcraft_compatible"]:
            log(
                f"ERROR: Aborting BindCraft generation on incompatible GPU {gpu_info['gpu_name']} "
                f"({gpu_info['sm_arch']}) to prevent wasting Compute Units. Switch runtime to A100!"
            )
            return False

        if stage == "evaluation_calibration":
            ok = run_calibration_job(job, run_dir, log, dry_run=dry_run)
            return ok

        target_rel = job.get("target_pdb_relpath", "")
        target_pdb = (REPO_ROOT / target_rel).resolve()
        if not target_pdb.exists():
            alt_pdb = (drive_root / target_rel).resolve()
            if alt_pdb.exists():
                target_pdb = alt_pdb
            else:
                log(f"ERROR: Target PDB not found at {target_pdb} or {alt_pdb}")
                return False

        pre_ok, val_out = run_pre_validation(
            target_pdb=target_pdb,
            hotspots=job.get("hotspots", ""),
            expected_sha=job.get("expected_sha256", ""),
        )
        log("Pre-validation output:\n" + val_out.strip())
        if not pre_ok:
            log("ERROR: Pre-validation failed. Aborting job.")
            return False

        settings_path = write_bindcraft_settings(job, target_pdb, run_dir)
        log(f"Wrote BindCraft target settings to {settings_path}")

        if dry_run:
            log("Dry-run mode enabled: skipping GPU BindCraft invocation.")
        else:
            bindcraft_script = pathlib.Path("/content/BindCraft/bindcraft.py")
            if bindcraft_script.exists():
                log(f"Launching BindCraft trajectories in {run_dir}...")
                bc_cmd = [
                    "conda",
                    "run",
                    "-n",
                    "BindCraft",
                    "python",
                    "-u",
                    str(bindcraft_script),
                    "--settings",
                    str(settings_path),
                    "--filters",
                    "/content/BindCraft/settings_filters/default_filters.json",
                    "--advanced",
                    "/content/BindCraft/settings_advanced/default_4stage_multimer.json",
                ]
                proc = subprocess.Popen(
                    bc_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", bufsize=1
                )
                assert proc.stdout is not None
                for line in proc.stdout:
                    clean = line.rstrip("\r\n")
                    if clean:
                        log(clean)
                proc.wait()
                if proc.returncode != 0:
                    log(f"WARNING: BindCraft exited with return code {proc.returncode}; harvesting any partial outputs.")
            else:
                log(f"Ready for BindCraft execution in {run_dir} (/content/BindCraft not installed in this environment).")

        harvested = harvest_bindcraft_candidates(
            job, target_pdb, run_dir, dirs["registry"] / "design_registry.csv", dry_run=dry_run
        )
        log(f"Harvested & evaluated {harvested} candidate sequences into design_registry.csv.")
        ok = True
        return True
    except KeyboardInterrupt:
        interrupted = True
        log(f"Job {job_id} interrupted by user (KeyboardInterrupt). Keeping in jobs/running/ for auto-resume.")
        raise
    except Exception as exc:
        log(f"FATAL EXCEPTION in job {job_id}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
        ok = False
        return False
    finally:
        elapsed_sec = round(time.time() - t_start, 1)
        est_cu_used = round((elapsed_sec / 3600.0) * float(gpu_info["est_cu_per_hour"]), 3)
        manifest = {
            "job_id": job_id,
            "stage": stage,
            "status": "INTERRUPTED_RESUMABLE" if interrupted else ("COMPLETED" if ok else "FAILED"),
            "gpu_name": gpu_info["gpu_name"],
            "vram_gb": gpu_info["vram_gb"],
            "sm_arch": gpu_info["sm_arch"],
            "elapsed_seconds": elapsed_sec,
            "estimated_cu_used": est_cu_used,
            "finished_at": dt.datetime.now().isoformat(timespec="seconds"),
        }
        try:
            (run_dir / "run_manifest.json").write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        except Exception:
            pass
        if not interrupted:
            dest_dir = dirs["completed"] if ok else dirs["failed"]
            if running_path.exists():
                shutil.move(str(running_path), str(dest_dir / job_path.name))
            log(
                f"Job {job_id} finalized (ok={ok}, elapsed={elapsed_sec}s, est_CU={est_cu_used}) -> "
                f"{dest_dir / job_path.name}"
            )
        sync_filesystem()


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
    parser.add_argument(
        "--unassign-on-exit",
        action="store_true",
        help="Flush Google Drive and immediately release Colab GPU runtime (runtime.unassign()) after finishing",
    )
    parser.add_argument(
        "--diagnose",
        action="store_true",
        help="Print Google Drive queue and Target 2 calibration completeness diagnostic and exit",
    )
    parser.add_argument("--poll-interval", type=int, default=15, help="Seconds between queue polls")
    args = parser.parse_args()

    drive_root = pathlib.Path(args.drive_root)
    dirs = ensure_drive_workspace(drive_root)

    if args.diagnose:
        print_drive_diagnostic(drive_root)
        return

    exit_code = 0
    jobs_completed = 0
    any_failed = False
    try:
        if args.job_file:
            ok = process_job(pathlib.Path(args.job_file), drive_root, dry_run=args.dry_run)
            if ok:
                jobs_completed += 1
            else:
                any_failed = True
                exit_code = 1
            return

        print(f"Queue Worker watching: {dirs['pending']} (and auto-resuming from {dirs['running']})", flush=True)
        auto_recovered_once = False
        while True:
            resumed_jobs = sorted(dirs["running"].glob("*.json"))
            if resumed_jobs:
                print(
                    f"[AUTO-RESUME] 檢測到 {len(resumed_jobs)} 個上次中斷的任務於 {dirs['running']}，正在從斷點繼續執行...",
                    flush=True,
                )
            pending_jobs = sorted(dirs["pending"].glob("*.json"))
            all_queue = resumed_jobs + pending_jobs

            if not all_queue and not auto_recovered_once:
                auto_recovered_once = True
                # Self-healing check: did an earlier run move job_target2_TNFa_calibration.json to completed/failed prematurely?
                cal_dir = dirs["runs"] / "target2_TNFa_calibration"
                comp_job = dirs["completed"] / "job_target2_TNFa_calibration.json"
                fail_job = dirs["failed"] / "job_target2_TNFa_calibration.json"
                tmpl_job = REPO_ROOT / "jobs" / "templates" / "job_target2_TNFa_calibration.json"
                if cal_dir.exists() or comp_job.exists() or fail_job.exists():
                    is_done, _, missing_items = check_calibration_completeness(cal_dir)
                    if not is_done and tmpl_job.exists():
                        recovered = dirs["running"] / "job_target2_TNFa_calibration.json"
                        shutil.copy2(str(tmpl_job), str(recovered))
                        for old_p in (comp_job, fail_job):
                            if old_p.exists():
                                old_p.unlink()
                        print(
                            f"[AUTO-RECOVER] 檢測到上次 TNF-α 校準任務未跑完全部階段（待完成：{missing_items}）！\n"
                            f"               已自動將任務從 completed/failed 恢復至 {recovered} 並從斷點繼續執行...",
                            flush=True,
                        )
                        all_queue = [recovered]
                    elif is_done and args.once:
                        print_drive_diagnostic(drive_root)
                        print(
                            "🎉 [100% 已完成] Target 2 (TNF-α) 全部校準實驗已在上次運行中 100% 完成並保存在 Google Drive！無需重複消耗 GPU 額度。",
                            flush=True,
                        )

            if not all_queue and args.once:
                print(
                    f"[提示] 當前隊列為空 ({dirs['pending']} 與 {dirs['running']} 均無待執行任務)。"
                    "若需投遞新任務，請運行『第 3 格』！",
                    flush=True,
                )
            for jp in all_queue:
                ok = process_job(jp, drive_root, dry_run=args.dry_run)
                if ok:
                    jobs_completed += 1
                else:
                    any_failed = True
            if args.once:
                break
            time.sleep(args.poll_interval)
    except KeyboardInterrupt:
        any_failed = True
        print("\n[中斷] 收到手動停止信號 (KeyboardInterrupt)，任務已保留在 jobs/running/ 以便下次斷點續跑。", flush=True)
    finally:
        if args.unassign_on_exit:
            if not args.dry_run and jobs_completed > 0 and not any_failed:
                flush_drive_and_unassign_gpu(drive_root, unassign=True, countdown_sec=15)
            else:
                sync_filesystem()
                print(
                    f"\n[SAFEGUARD] 未觸發自動斷開 GPU (completed={jobs_completed}, failed_or_stopped={any_failed}, dry_run={args.dry_run})。"
                    "\n   當前 Colab GPU 連線已保留，方便你檢查日誌或重新投遞任務！",
                    flush=True,
                )

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
