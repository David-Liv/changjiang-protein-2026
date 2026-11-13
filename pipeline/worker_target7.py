#!/usr/bin/env python3
"""Changjiang Cup 2026 — Isolated Target 7 (IL-7Rα / CD127) Queue Worker & Harvester.

Designed to run concurrently in a second Google Colab session WITHOUT touching `main`'s
active Target 5 queue (`jobs/`) or shared registry (`registry/design_registry.csv`):

- Isolated Queue Root   : `<drive_root>/jobs_target7_IL7Ra/{pending,running,completed,failed}/`
- Isolated Registry CSV : `<drive_root>/registry/design_registry_T7_IL7Ra.csv`
- Local HF Weight Cache : `/content/hf_cache` (seeded read-only from `<drive_root>/cache/huggingface`
                          to eliminate Google Drive FUSE `.locks/` contention across concurrent Colabs)
- Shared Read-Only Tar  : `<drive_root>/cache/bindcraft_colab_cache.tar.zst` (45-second restore)

Enforces all Target 7 competition & biophysical gates:
1. SOP §B-1 Pre-flight: `scripts/validate_input.py --strict-sidechains` + SHA-256 verification.
2. SOP §C-4 Inverse folding: `SolubleMPNN` (`mpnn_weights="soluble"`, `omit_AAs="CM"` / `"CXM"`).
3. Official format & SOP §C-2 Expressibility: `40–100 aa` hard limit, 20 standard AAs, even Cys,
   `GRAVY <= 0.10`, hydrophobic run `< 8`, homopolymer `< 5`, no `WW`, `Ala <= 20%`.
4. SOP §C-1 Novelty hard gate: `<= 70.0%` sequence identity (when alignment covers `>= 80%`)
   against `targets/target7_IL7Ra/targets/known_binders_IL7Ra_SEQRES.fasta` (`7OPB_D`, `IL-7`, `TSLP`).
5. 3D `3DI3` Crystal Audit: Kabsch superposition onto `3DI3` Chain B (`2.15 Å`) to reject clashes
   with crystallographic N-glycans (`Chains C/D/E` at `Asn29/Asn45/Asn131`, `< 2.5 Å`) and verify
   direct orthosteric competition with `IL-7` (`Chain A`, `< 3.5 Å`).
6. Stage-C Complex Screening & Export: `ESMFold + ipSAE + pDockQ` against the full official 219-aa
   `P16871 (Glu21–Asp239)` ECD construct, plus `1:1` AlphaFold 3 Server JSON export.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
import traceback
from typing import Any

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT, REPO_ROOT / "pipeline", REPO_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import check_novelty  # noqa: E402
import tnf_pipeline as tp  # noqa: E402

SHARED_DRIVE_FOLDER_URL = "https://drive.google.com/drive/folders/1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg?usp=drive_link"
SHARED_DRIVE_FOLDER_ID = "1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg"

DEFAULT_T7_QUEUE_SUBDIR = "jobs_target7_IL7Ra"
DEFAULT_T7_REGISTRY_FILENAME = "design_registry_T7_IL7Ra.csv"
OFFICIAL_IL7RA_ECD_21_239 = (
    "ESGYAQNGDLEDAELDDYSFSCYSQLEVNGSQHSLTCAFEDPDVNITNLEFEICGALVEVKCLNFRKLQEIYFIETKKFLLIGKSNICVKVGEKSLTCKK"
    "IDLTTIVKPEAPFDLSVVYREGANDFVVTFNTSHLQKKYVKVLMHDVAYRQEKDENKWTHVNLSSTKLTLLQRKLQPAAMYEIKVRSIPDHYFKGFWSEW"
    "SPSYYFRTPEINNSSGEMD"
)

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

AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",
}
GLYCAN_RESNAMES = frozenset({"NAG", "NDG", "BMA", "MAN", "FUC", "GAL", "SIA", "GLC"})
BACKBONE_ATOMS = ("N", "CA", "C", "O")
CORE_ANCHOR_AAS = frozenset("FYWLIV")


def _is_hydrogen(atom_name: str, elem: str) -> bool:
    return (
        elem == "H"
        or atom_name.startswith("H")
        or (len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == "H")
    )


def _extract_drive_folder_id(spec: str) -> str | None:
    s = str(spec).strip()
    m = re.search(r"folders/([a-zA-Z0-9_-]{15,})", s)
    if m:
        return m.group(1)
    if re.fullmatch(r"[a-zA-Z0-9_-]{25,45}", s):
        return s
    if s.rstrip("/") == "/content/drive/MyDrive/Changjiang_Cup_2026":
        return SHARED_DRIVE_FOLDER_ID
    return None


def resolve_shared_drive_workspace(drive_spec: str | pathlib.Path) -> pathlib.Path:
    """Resolve the shared Google Drive folder (ID: 1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg) without creating private folders."""
    spec_str = str(drive_spec).strip()
    folder_id = _extract_drive_folder_id(spec_str)
    colab_drive = pathlib.Path("/content/drive")

    if not colab_drive.exists() and not spec_str.startswith("/content/drive") and not spec_str.startswith("http"):
        return pathlib.Path(spec_str)

    if folder_id:
        if (pathlib.Path("/content").exists() or "google.colab" in sys.modules) and not (colab_drive / "MyDrive").exists():
            try:
                from google.colab import drive as _colab_drive  # type: ignore

                print("🔄 偵測到 Colab 重連後尚未掛載 Google Drive，正在自動掛載 /content/drive ...", flush=True)
                _colab_drive.mount("/content/drive")
            except Exception as _mount_exc:
                raise RuntimeError(
                    f"❌ Google Drive 尚未掛載，且自動掛載失敗 ({_mount_exc})！請先執行 drive.mount('/content/drive')。"
                ) from _mount_exc

        mydrive = colab_drive / "MyDrive"
        if mydrive.exists():
            try:
                list(mydrive.iterdir())
            except Exception:
                pass

        id_target = colab_drive / ".shortcut-targets-by-id" / folder_id
        resolved: pathlib.Path | None = None
        if id_target.exists():
            subdirs = [p for p in id_target.iterdir() if p.is_dir()]
            resolved = subdirs[0] if len(subdirs) == 1 else id_target
        elif mydrive.exists():
            for item in mydrive.iterdir():
                try:
                    if item.is_symlink() and folder_id in str(item.resolve()):
                        resolved = item.resolve()
                        break
                except Exception:
                    continue

        if resolved is None or not resolved.exists():
            raise FileNotFoundError(
                f"\n❌ 未能在當前掛載的 Google Drive 中定位到團隊共享資料夾 (ID: {folder_id})！\n"
                f"   共享連結：{SHARED_DRIVE_FOLDER_URL}\n"
                "   請先在 Google Drive 將該資料夾『新增捷徑至我的雲端硬碟』後重試。"
            )
        return resolved

    p = pathlib.Path(spec_str)
    if str(p).startswith("/content/drive") and not p.exists():
        raise FileNotFoundError(f"❌ 指定的 Google Drive 路徑不存在：{p}")
    return p


def ensure_t7_workspace(
    drive_root: pathlib.Path,
    queue_subdir: str = DEFAULT_T7_QUEUE_SUBDIR,
    registry_filename: str = DEFAULT_T7_REGISTRY_FILENAME,
    allow_shared_queue: bool = False,
) -> dict[str, pathlib.Path]:
    """Create isolated Target 7 queue & registry folders without touching T5's `jobs/` or `design_registry.csv`."""
    if not allow_shared_queue:
        if queue_subdir.strip("/\\") == "jobs":
            raise ValueError(
                "Refusing to use shared `jobs/` queue in worker_target7.py so active T5 runs are never disturbed! "
                f"Use `{DEFAULT_T7_QUEUE_SUBDIR}`."
            )
        if registry_filename.strip() == "design_registry.csv":
            raise ValueError(
                "Refusing to write to shared `design_registry.csv` while T5 is active! "
                f"Use `{DEFAULT_T7_REGISTRY_FILENAME}`."
            )

    drive_root = resolve_shared_drive_workspace(drive_root)
    if not str(drive_root).startswith("/content/drive"):
        drive_root.mkdir(parents=True, exist_ok=True)
    elif not drive_root.exists():
        raise FileNotFoundError(f"Shared Google Drive root does not exist: {drive_root}")

    dirs = {
        "root": drive_root,
        "queue_root": drive_root / queue_subdir,
        "pending": drive_root / queue_subdir / "pending",
        "running": drive_root / queue_subdir / "running",
        "completed": drive_root / queue_subdir / "completed",
        "failed": drive_root / queue_subdir / "failed",
        "runs": drive_root / "runs",
        "registry": drive_root / "registry",
        "registry_csv": drive_root / "registry" / registry_filename,
        "cache": drive_root / "cache",
    }
    for k, d in dirs.items():
        if k not in ("root", "registry_csv"):
            d.mkdir(parents=True, exist_ok=True)

    # Use local NVMe /content/hf_cache on Colab (seeded read-only from Drive cache if available)
    # so concurrent T5 and T7 Colabs NEVER contend on Google Drive FUSE `.locks/`!
    if str(drive_root).startswith("/content/drive") and pathlib.Path("/content").exists():
        local_hf = pathlib.Path("/content/hf_cache")
        drive_hf = dirs["cache"] / "huggingface"
        if not local_hf.exists() and drive_hf.exists():
            try:
                print("⚡ [HF Cache] 正在將網盤已快取的 ESMFold 權重唯讀複製至本地 NVMe (/content/hf_cache)，避免與 T5 爭搶網盤鎖...", flush=True)
                shutil.copytree(
                    str(drive_hf),
                    str(local_hf),
                    ignore=shutil.ignore_patterns("*.lock", ".locks"),
                    dirs_exist_ok=True,
                )
            except Exception as exc:
                print(f"[HF Cache Info] 本地預載略過 ({exc})，將在本地 /content/hf_cache 自動下載。", flush=True)
        local_hf.mkdir(parents=True, exist_ok=True)
        os.environ["HF_HOME"] = str(local_hf)

    reg_csv = dirs["registry_csv"]
    if not reg_csv.exists():
        with reg_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(REGISTRY_HEADERS)
    return dirs


def inspect_gpu_environment() -> dict[str, Any]:
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
                warn = "H100 burns ~18 CU/hr. A100 (~5.4-7.5 CU/hr) is recommended for better CU efficiency."
                bc_compat = True
            elif major >= 12 or "BLACKWELL" in name_up or "RTX PRO 6000" in name_up or "G4" in name_up:
                cu_rate = 9.5
                warn = f"Detected Blackwell {p.name} ({sm}). BindCraft CUDA 12.4/JAX does not support {sm}; switch to A100!"
                bc_compat = False
            elif "A100" in name_up or (major == 8 and minor == 0):
                cu_rate = 7.5 if gb >= 70 else 5.4
                warn = ""
                bc_compat = True
            else:
                cu_rate = 2.0
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
    info = inspect_gpu_environment()
    return 0 if info["vram_gb"] >= 38 else 64


def sync_filesystem() -> None:
    sys.stdout.flush()
    sys.stderr.flush()
    if hasattr(os, "sync"):
        try:
            os.sync()
        except Exception:
            pass


def _has_active_ipython_kernel() -> bool:
    try:
        import IPython  # type: ignore

        ip = IPython.get_ipython()
        return ip is not None and getattr(ip, "kernel", None) is not None
    except Exception:
        return False


def flush_drive_and_unassign_gpu(drive_root: pathlib.Path, unassign: bool = True, countdown_sec: int = 15) -> None:
    sync_filesystem()
    if not drive_root.exists():
        print(f"[SAFEGUARD] Drive root {drive_root} does not exist! Skipping GPU unassign.")
        return

    t7_runs = [p for p in (drive_root / "runs").glob("job_t7_*") if p.is_dir()] if (drive_root / "runs").exists() else []
    saved_files = [f for d in t7_runs for f in d.rglob("*") if f.is_file()]
    print(f"\n=== [Auto-Save Verification] {len(saved_files)} Target 7 files persisted across {len(t7_runs)} run(s) ===")
    if len(saved_files) == 0:
        print("[SAFEGUARD] No Target 7 output files found in runs/job_t7_* — skipping GPU unassign!")
        return

    if unassign and not _has_active_ipython_kernel():
        print("[Info] Deferring drive.flush_and_unmount() and runtime.unassign() to parent Colab Notebook cell.", flush=True)
        return

    if unassign and countdown_sec > 0:
        print(
            f"\n⏳ Target 7 任務已全部完成並寫入 Google Drive！將在 {countdown_sec} 秒後自動釋放當前 Colab 的 GPU..."
            "\n   （不影響另一個正在跑 T5 的 Colab 連線；若想保留當前連線查看結果，請按『停止 ■』鍵）",
            flush=True,
        )
        try:
            for remaining in range(countdown_sec, 0, -5):
                print(f"   倒數 {remaining} 秒...", flush=True)
                time.sleep(min(5, remaining))
        except KeyboardInterrupt:
            print("\n✋ 已取消自動斷開 GPU，當前 Colab 連線保留。")
            return

    sync_filesystem()
    if unassign:
        try:
            from google.colab import drive as _colab_drive  # type: ignore

            _colab_drive.flush_and_unmount()
        except Exception:
            pass
        try:
            from google.colab import runtime  # type: ignore

            time.sleep(2)
            runtime.unassign()
        except Exception:
            pass


def print_t7_diagnostic(
    drive_root: pathlib.Path,
    queue_subdir: str = DEFAULT_T7_QUEUE_SUBDIR,
    registry_filename: str = DEFAULT_T7_REGISTRY_FILENAME,
) -> None:
    dirs = ensure_t7_workspace(drive_root, queue_subdir=queue_subdir, registry_filename=registry_filename)
    print("\n================ [Target 7 (IL-7Rα) 獨立隊列與註冊表診斷] ================")
    print(f"  🔒 隔離隊列目錄 : {dirs['queue_root']} （與 T5 的 jobs/ 100% 物理隔離）")
    print(f"  📒 隔離註冊表   : {dirs['registry_csv']}")
    for k in ("running", "pending", "completed", "failed"):
        files = [p.name for p in sorted(dirs[k].glob("*.json"))]
        print(f"  📂 {queue_subdir}/{k:<10}: {len(files)} 個 -> {files}")

    reg_csv = dirs["registry_csv"]
    if reg_csv.exists():
        with reg_csv.open("r", newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        n_pass = sum(1 for r in rows if r.get("expressibility_verdict") == "PASS")
        print(f"  🧬 已登錄 Target 7 序列總數: {len(rows)} 條（其中全部門控 PASS = {n_pass} 條）")

    t7_runs = sorted((dirs["runs"]).glob("job_t7_*"))
    if t7_runs:
        print("  🏃 Target 7 運行目錄摘要:")
        for rdir in t7_runs[-6:]:
            n_traj = len(list((rdir / "Trajectory" / "Relaxed").glob("*.pdb"))) if (rdir / "Trajectory" / "Relaxed").exists() else 0
            n_acc = len(list((rdir / "Accepted").rglob("*.pdb"))) if (rdir / "Accepted").exists() else 0
            n_af3 = len(list((rdir / "af3_inputs").glob("*.json"))) if (rdir / "af3_inputs").exists() else 0
            status = "IN_PROGRESS"
            mf = rdir / "run_manifest.json"
            if mf.exists():
                try:
                    status = json.loads(mf.read_text(encoding="utf-8")).get("status", status)
                except Exception:
                    pass
            print(f"     - {rdir.name}: status={status}, Relaxed={n_traj} PDB, Accepted={n_acc} PDB, AF3_JSON={n_af3}")
    print("==========================================================================\n", flush=True)


def run_pre_validation(target_pdb: pathlib.Path, hotspots: str, expected_sha: str) -> tuple[bool, str]:
    """Run SOP §B-1 input validator with `--strict-sidechains`."""
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "validate_input.py"),
        "--target",
        str(target_pdb),
        "--hotspots",
        hotspots or "",
        "--strict-sidechains",
    ]
    if expected_sha:
        cmd.extend(["--sha", expected_sha])
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    output = (proc.stdout or "") + (proc.stderr or "")
    return proc.returncode == 0, output


def parse_pdb_residues_and_glycans(
    pdb_path: pathlib.Path,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, np.ndarray]]:
    """Parse heavy-atom protein residues grouped by chain AND non-water glycan HETATM heavy atoms."""
    by_chain: dict[str, list[dict[str, Any]]] = {}
    cur_map: dict[tuple[str, str], dict[str, Any]] = {}
    glycan_pts: dict[str, list[np.ndarray]] = {}

    for ln in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        if len(ln) < 54 or ln[16] not in (" ", "A"):
            continue
        atom_name = ln[12:16].strip()
        elem = ln[76:78].strip() if len(ln) >= 78 else ""
        if _is_hydrogen(atom_name, elem):
            continue
        res3 = ln[17:20].strip()
        ch = ln[21]
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])], dtype=float)

        if ln.startswith("ATOM  ") and res3 in AA3_TO_1:
            res_uid = ln[22:27]
            resseq = int(ln[22:26])
            key = (ch, res_uid)
            rec = cur_map.get(key)
            if rec is None:
                rec = {
                    "chain": ch,
                    "resseq": resseq,
                    "aa": AA3_TO_1[res3],
                    "bb": {},
                    "heavy_atoms": [],
                }
                cur_map[key] = rec
                by_chain.setdefault(ch, []).append(rec)
            rec["heavy_atoms"].append(xyz)
            if atom_name in BACKBONE_ATOMS and atom_name not in rec["bb"]:
                rec["bb"][atom_name] = xyz
        elif ln.startswith("HETATM") and res3 in GLYCAN_RESNAMES:
            glycan_pts.setdefault(ch, []).append(xyz)

    glycan_arrays = {c: np.vstack(pts) for c, pts in glycan_pts.items() if pts}
    return by_chain, glycan_arrays


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    P0 = P - cP
    Q0 = Q - cQ
    U, _, Vt = np.linalg.svd(P0.T @ Q0)
    d = float(np.linalg.det(Vt.T @ U.T))
    S = np.diag([1.0, 1.0, -1.0 if d < 0.0 else 1.0])
    R = Vt.T @ S @ U.T
    aligned = P0 @ R.T
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - Q0) ** 2, axis=1))))
    return R, cP, cQ, rmsd


def audit_3di3_glycan_and_il7(
    complex_pdb: pathlib.Path,
    binder_seq: str,
    ref_3di3_pdb: pathlib.Path,
    ref_receptor_chain: str = "B",
    ref_ligand_chain: str = "A",
) -> dict[str, Any]:
    """Superimpose candidate complex onto `3DI3` receptor Chain B and measure N-glycan clashes & IL-7 overlap."""
    out = {
        "aligned_ok": False,
        "receptor_rmsd": 999.0,
        "glycan_clashes": 0,
        "min_glycan_dist": 999.0,
        "il7_overlap_atoms": -1,
    }
    if not complex_pdb.exists() or not ref_3di3_pdb.exists():
        return out

    cmp_chains, _ = parse_pdb_residues_and_glycans(complex_pdb)
    ref_chains, ref_glycans = parse_pdb_residues_and_glycans(ref_3di3_pdb)
    if len(cmp_chains) < 2 or ref_receptor_chain not in ref_chains:
        return out

    chain_seqs = {c: "".join(r["aa"] for r in rl if "CA" in r["bb"]) for c, rl in cmp_chains.items()}
    binder_c = next((c for c, s in chain_seqs.items() if s == binder_seq), None)
    if binder_c is None:
        binder_c = min(chain_seqs.keys(), key=lambda c: abs(len(chain_seqs[c]) - len(binder_seq)))

    target_candidates = [c for c in cmp_chains.keys() if c != binder_c]
    if not target_candidates:
        return out
    target_c = max(target_candidates, key=lambda c: len(cmp_chains[c]))

    tgt_ca_res = [r for r in cmp_chains[target_c] if "CA" in r["bb"]]
    ref_ca_res = [r for r in ref_chains[ref_receptor_chain] if "CA" in r["bb"]]
    seq_tgt = "".join(r["aa"] for r in tgt_ca_res)
    seq_ref = "".join(r["aa"] for r in ref_ca_res)
    if len(seq_tgt) < 80 or len(seq_ref) < 80:
        return out

    aln_t, aln_r = check_novelty.nw(seq_tgt, seq_ref)
    P_pts: list[np.ndarray] = []
    Q_pts: list[np.ndarray] = []
    it, ir = 0, 0
    for xt, xr in zip(aln_t, aln_r):
        if xt != "-" and xr != "-":
            P_pts.append(tgt_ca_res[it]["bb"]["CA"])
            Q_pts.append(ref_ca_res[ir]["bb"]["CA"])
        if xt != "-":
            it += 1
        if xr != "-":
            ir += 1

    if len(P_pts) < 120:
        return out

    R, cP, cQ, rmsd = kabsch_fit(np.vstack(P_pts), np.vstack(Q_pts))
    out["receptor_rmsd"] = round(rmsd, 3)
    if rmsd > 3.5:
        return out
    out["aligned_ok"] = True

    binder_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in cmp_chains[binder_c]])
    binder_in_3di3 = (binder_heavy - cP) @ R.T + cQ

    if ref_glycans:
        all_glycan_pts = np.vstack(list(ref_glycans.values()))
        d_gly = np.linalg.norm(binder_in_3di3[:, None, :] - all_glycan_pts[None, :, :], axis=2)
        out["glycan_clashes"] = int((d_gly.min(axis=1) < 2.5).sum())
        out["min_glycan_dist"] = round(float(d_gly.min()), 2)

    if ref_ligand_chain in ref_chains:
        il7_pts = np.vstack([np.asarray(r["heavy_atoms"]) for r in ref_chains[ref_ligand_chain]])
        d_il7 = np.linalg.norm(binder_in_3di3[:, None, :] - il7_pts[None, :, :], axis=2)
        out["il7_overlap_atoms"] = int((d_il7.min(axis=1) < 3.5).sum())

    return out


def check_sequence_novelty_t7(
    seq: str,
    known_fasta: pathlib.Path,
    max_identity_pct: float = 70.0,
    min_cov_ratio: float = 0.80,
) -> tuple[bool, float, str]:
    """Verify Changjiang Cup novelty gate (`<= 70%` identity when alignment covers `>= 80%`).

    Checks BOTH global Needleman-Wunsch alignment AND sliding-window fragments on longer cytokines
    (`IL-7`, `TSLP`) so grafted cytokine helices are also caught.
    """
    s = seq.strip().upper()
    if not s or not known_fasta.exists():
        return True, 0.0, ""

    refs = check_novelty.read_fasta(str(known_fasta))
    worst_pct = 0.0
    worst_name = ""

    for ref_name, ref_seq_raw in refs.items():
        ref_seq = ref_seq_raw.strip().upper()
        if not ref_seq:
            continue

        # 1) Global Needleman-Wunsch alignment
        same, aln, pct, _ = check_novelty.ident(s, ref_seq)
        cov = aln / max(1, min(len(s), len(ref_seq)))
        if cov >= min_cov_ratio and pct > worst_pct:
            worst_pct = pct
            worst_name = ref_name

        # 2) Ungapped + local window check if reference is longer than candidate (e.g. IL-7 152 aa, TSLP 126 aa)
        if len(ref_seq) > len(s):
            win = len(s)
            for start in range(0, len(ref_seq) - win + 1):
                sub = ref_seq[start : start + win]
                ungapped_pct = 100.0 * sum(1 for a, b in zip(s, sub) if a == b) / win
                if ungapped_pct > worst_pct:
                    worst_pct = ungapped_pct
                    worst_name = ref_name

    is_novel = worst_pct <= max_identity_pct
    return is_novel, round(worst_pct, 2), worst_name


def compute_expressibility_metrics_t7(
    seq: str,
    min_len: int = 40,
    max_len: int = 100,
    max_ala_ratio: float = 0.20,
) -> dict[str, Any]:
    """Compute SOP §C-2 composition + official `40-100 aa` minibinder gate."""
    pka_pos = {"K": 10.54, "R": 12.48, "H": 6.04}
    pka_neg = {"D": 3.90, "E": 4.07, "Y": 10.46, "C": 8.18}
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    s = seq.strip().upper()
    n = len(s)
    if n == 0:
        return {"expressibility_verdict": "REJECT:EMPTY"}

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
    ala_ratio = s.count("A") / n
    cys = s.count("C")
    cys_ok = (cys % 2 == 0)
    runs = max((len(m.group(0)) for m in re.finditer(r"[AVLIMFWYC]+", s)), default=0)
    homo = max((len(m.group(0)) for m in re.finditer(r"(.)\1+", s)), default=1)

    reasons: list[str] = []
    if not set(s) <= valid_aa:
        reasons.append("NON_STANDARD_AA")
    if n < min_len or n > max_len:
        reasons.append(f"LENGTH_OUT_OF_{min_len}_{max_len}")
    if not cys_ok:
        reasons.append("ODD_CYS")
    if runs >= 8:
        reasons.append(f"HYDROPHOBIC_RUN_{runs}")
    if homo >= 5:
        reasons.append(f"HOMOPOLYMER_{homo}")
    if "WW" in s:
        reasons.append("ADJACENT_WW")
    if ala_ratio > max_ala_ratio:
        reasons.append(f"HIGH_ALA_{ala_ratio:.2f}")
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
    """Fold `binder : 219-aa official IL-7Rα` with ESMFold and compute Dunbrack ipSAE + pDockQ."""
    if not candidates or not target_seq or dry_run:
        return {}

    eval_dir.mkdir(parents=True, exist_ok=True)
    spec_path = eval_dir / "candidates_eval.tsv"
    with spec_path.open("w", encoding="utf-8") as f:
        f.write("# design_id\tbinder:target\n")
        for cid, seq in candidates:
            f.write(f"{cid}\t{seq}:{target_seq}\n")

    chunk = detect_gpu_chunk()
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "tnf_pipeline.py"),
            "fold",
            "--spec",
            str(spec_path),
            "--out",
            str(eval_dir),
            "--chunk",
            str(chunk),
        ],
        check=False,
    )

    summary_csv = eval_dir / "ipsae_summary.csv"
    subprocess.run(
        [
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
        ],
        check=False,
    )

    results: dict[str, dict[str, Any]] = {}
    if summary_csv.exists():
        with summary_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                name = row.get("name", "")
                if name:
                    results[name] = {
                        "esm_plddt": row.get("mean_plddt", ""),
                        "ipsae10_min": row.get("ipsae10_min", ""),
                        "ipsae10_max": row.get("ipsae10_merged_max") or row.get("ipsae10_max", ""),
                        "pdockq": row.get("iface_pdockq", ""),
                        "min_heavy_dist": row.get("iface_min_heavy_dist", ""),
                    }
    return results


def export_af3_and_submission_preview(
    records: list[dict[str, Any]],
    job: dict[str, Any],
    run_dir: pathlib.Path,
    batch_size: int = 20,
    suffix_tag: str = "",
) -> list[pathlib.Path]:
    """Export PASS candidates as `1x Binder + 1x Official 219-aa IL-7Rα` AlphaFold 3 Server JSONs + Drive preview CSV."""
    passing = [r for r in records if r.get("expressibility_verdict") == "PASS"]
    if not passing:
        return []

    af3_tgt_seq = (job.get("af3_target_sequence") or OFFICIAL_IL7RA_ECD_21_239).strip().upper()
    af3_copies = int(job.get("af3_target_copies", 1))

    def _sort_key(r: dict[str, Any]) -> tuple[float, float]:
        try:
            ips = float(r.get("ipsae10_max") or 0.0)
        except ValueError:
            ips = 0.0
        try:
            iptm = float(r.get("fold_ipTM") or 0.0)
        except ValueError:
            iptm = 0.0
        return (-ips, -iptm)

    passing.sort(key=_sort_key)
    out_dir = run_dir / "af3_inputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    job_id = job.get("job_id", "job_t7_il7ra")
    written_files: list[pathlib.Path] = []

    for part_idx in range(0, len(passing), batch_size):
        chunk = passing[part_idx : part_idx + batch_size]
        payload = []
        for rec in chunk:
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(rec["design_id"]))
            payload.append(
                {
                    "name": safe_name,
                    "modelSeeds": [],
                    "sequences": [
                        {"proteinChain": {"sequence": rec["sequence"], "count": 1}},
                        {"proteinChain": {"sequence": af3_tgt_seq, "count": af3_copies}},
                    ],
                    "dialect": "alphafoldserver",
                    "version": 1,
                }
            )
        part_s = f"_part{part_idx // batch_size + 1}" if len(passing) > batch_size else ""
        out_p = out_dir / f"af3_server_{job_id}{suffix_tag}{part_s}.json"
        out_p.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        written_files.append(out_p)

    # Also write an official-format Drive-only top-20 preview CSV (`name,sequence`) inside Google Drive `runs/<job_id>/`
    prev_dir = run_dir / "submission_preview"
    prev_dir.mkdir(parents=True, exist_ok=True)
    prev_csv = prev_dir / f"IL7Ra_candidates_{job_id}{suffix_tag}_top20.csv"
    with prev_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "sequence"])
        for rec in passing[:20]:
            w.writerow([rec["design_id"], rec["sequence"]])
    return written_files


def write_bindcraft_settings(
    job: dict[str, Any], target_pdb: pathlib.Path, run_dir: pathlib.Path
) -> tuple[pathlib.Path, pathlib.Path]:
    lengths = job.get("lengths", [55, 85])
    hotspots = job.get("hotspots")
    hotspots_val = "" if hotspots in ("", "null", "None", None) else str(hotspots)

    settings = {
        "design_path": str(run_dir).rstrip("/") + "/",
        "binder_name": job["job_id"],
        "starting_pdb": str(target_pdb),
        "chains": job.get("chains", "B"),
        "target_hotspot_residues": hotspots_val,
        "lengths": lengths,
        "number_of_final_designs": int(job.get("number_of_final_designs", 15)),
    }
    settings_path = run_dir / "target_settings.json"
    settings_path.write_text(json.dumps(settings, indent=2, ensure_ascii=False), encoding="utf-8")

    default_adv = pathlib.Path("/content/BindCraft/settings_advanced/default_4stage_multimer.json")
    adv_data: dict[str, Any] = {}
    if default_adv.exists():
        try:
            adv_data = json.loads(default_adv.read_text(encoding="utf-8"))
        except Exception:
            adv_data = {}
    adv_data["max_trajectories"] = int(job.get("max_trajectories", 35))
    adv_data["enable_mpnn"] = True
    adv_data["mpnn_weights"] = "soluble" if bool(job.get("use_soluble_mpnn", True)) else "original"
    if pathlib.Path("/content/BindCraft/params").exists():
        adv_data["af_params_dir"] = "/content/BindCraft/params/"
    adv_overrides = job.get("advanced_overrides")
    if isinstance(adv_overrides, dict):
        adv_data.update(adv_overrides)

    advanced_path = run_dir / "advanced_settings.json"
    advanced_path.write_text(json.dumps(adv_data, indent=2, ensure_ascii=False), encoding="utf-8")
    return settings_path, advanced_path


def find_candidate_complex_pdb(run_dir: pathlib.Path, design_id: str) -> pathlib.Path | None:
    search_dirs = [
        run_dir / "Accepted" / "Ranked",
        run_dir / "Accepted",
        run_dir / "MPNN" / "Relaxed",
        run_dir / "Trajectory" / "Relaxed",
    ]
    for sdir in search_dirs:
        if not sdir.exists():
            continue
        exact = sdir / f"{design_id}_model1.pdb"
        if exact.exists():
            return exact
        matches = sorted(sdir.glob(f"*{design_id}*.pdb"))
        if matches:
            return matches[0]
    return None


def ensure_protein_mpnn_repo() -> pathlib.Path:
    candidates = [
        pathlib.Path("/content/ProteinMPNN"),
        REPO_ROOT.parent / "ProteinMPNN",
    ]
    for c in candidates:
        if (c / "protein_mpnn_run.py").exists() and (c / "soluble_model_weights" / "v_48_020.pt").exists():
            return c
    dest = pathlib.Path("/content/ProteinMPNN") if pathlib.Path("/content").exists() else (REPO_ROOT / ".cache_protein_mpnn")
    print(f"[SolubleMPNN] Cloning official dauparas/ProteinMPNN into {dest} ...", flush=True)
    subprocess.run(["git", "clone", "--depth", "1", "https://github.com/dauparas/ProteinMPNN.git", str(dest)], check=True)
    return dest


def _select_t7_core_anchors(
    binder_res: list[dict[str, Any]],
    target_res: list[dict[str, Any]],
    target_labels: list[str],
    hotspots: set[str],
) -> list[int]:
    tgt_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in target_res])
    hs_list = [np.asarray(r["heavy_atoms"]) for i, r in enumerate(target_res) if target_labels[i] in hotspots]
    hs_heavy = np.vstack(hs_list) if hs_list else None

    cands: list[tuple[float, int, str]] = []
    for idx, bres in enumerate(binder_res):
        aa = bres["aa"]
        if aa not in CORE_ANCHOR_AAS:
            continue
        bh = np.asarray(bres["heavy_atoms"])
        d_tgt = float(np.sqrt(np.min(np.sum((bh[:, None, :] - tgt_heavy[None, :, :]) ** 2, axis=2))))
        if d_tgt > 4.0:
            continue
        d_hs = (
            float(np.sqrt(np.min(np.sum((bh[:, None, :] - hs_heavy[None, :, :]) ** 2, axis=2))))
            if hs_heavy is not None
            else d_tgt
        )
        cands.append((d_tgt + 0.15 * min(d_hs, 8.0), idx + 1, aa))

    cands.sort(key=lambda x: x[0])
    max_n = max(1, min(9, len(binder_res) // 6))
    chosen: list[int] = []
    w_pos: set[int] = set()
    for _, pos1, aa in cands:
        if len(chosen) >= max_n:
            break
        if aa == "W":
            if len(w_pos) >= 2 or (pos1 - 1) in w_pos or (pos1 + 1) in w_pos:
                continue
            w_pos.add(pos1)
        chosen.append(pos1)
    return sorted(chosen)


def harvest_relaxed_t7_trajectories(
    run_dir: pathlib.Path,
    job: dict[str, Any],
    registry_csv: pathlib.Path,
    min_iptm: float = 0.72,
    min_plddt: float = 0.78,
    num_seqs_core_fix: int = 6,
    num_seqs_full: int = 6,
    omit_aas: str = "CXM",
    max_per_backbone: int = 2,
    skip_esmfold: bool = False,
) -> int:
    """Harvest PyRosetta-relaxed single-chain IL-7Rα trajectories via SolubleMPNN + 3DI3 glycan/IL-7 + ESMFold gates."""
    traj_csv = run_dir / "trajectory_stats.csv"
    relaxed_dir = run_dir / "Trajectory" / "Relaxed"
    if not traj_csv.exists() or not relaxed_dir.exists():
        raise FileNotFoundError(f"Missing trajectory_stats.csv or Trajectory/Relaxed under {run_dir}")

    target_pdb = (REPO_ROOT / job["target_pdb_relpath"]).resolve()
    target_sha = hashlib.sha256(target_pdb.read_bytes()).hexdigest()
    tgt_chains, _ = parse_pdb_residues_and_glycans(target_pdb)
    orig_ch = sorted(tgt_chains.keys())[0]
    tgt_template_res = tgt_chains[orig_ch]
    target_labels = [f"{orig_ch}{r['resseq']}" for r in tgt_template_res]

    ref_3di3 = (REPO_ROOT / job.get("glycan_reference_pdb_relpath", "targets/target7_IL7Ra/evidence/3DI3.pdb")).resolve()
    ref_chains_3di3, _ = parse_pdb_residues_and_glycans(ref_3di3)
    wt_by_resseq = {r["resseq"]: ("I" if r["resseq"] == 46 else r["aa"]) for r in ref_chains_3di3.get("B", [])}
    wt_target_seq = "".join(wt_by_resseq.get(r["resseq"], r["aa"]) for r in tgt_template_res)
    masked_indices = {i for i, r in enumerate(tgt_template_res) if wt_by_resseq.get(r["resseq"], r["aa"]) != r["aa"]}

    known_fasta = (REPO_ROOT / job.get("known_binders_fasta_relpath", "targets/target7_IL7Ra/targets/known_binders_IL7Ra_SEQRES.fasta")).resolve()
    hotspots = {h.strip() for h in (job.get("hotspots") or "").split(",") if h.strip()}
    hard_lim = job.get("lengths_hard_limit", [40, 100])
    min_len_h, max_len_h = int(hard_lim[0]), int(hard_lim[1])

    traj_rows = list(csv.DictReader(traj_csv.open("r", newline="", encoding="utf-8")))
    passing_bbs: list[dict[str, Any]] = []
    backbone_ca_map: dict[str, np.ndarray] = {}

    for r in traj_rows:
        did = r.get("Design", "")
        iptm = float(r.get("i_pTM") or 0.0)
        plddt = float(r.get("pLDDT") or 0.0)
        ipae = float(r.get("i_pAE") or 1.0)
        pdb_p = relaxed_dir / f"{did}.pdb"
        if not did or iptm < min_iptm or plddt < min_plddt or not pdb_p.exists():
            continue

        parsed, _ = parse_pdb_residues_and_glycans(pdb_p)
        if "A" not in parsed or "B" not in parsed:
            continue
        t_res, b_res = parsed["A"], parsed["B"]
        if len(t_res) != len(tgt_template_res) or not (min_len_h <= len(b_res) <= max_len_h):
            continue
        if any(any(a not in x["bb"] for a in BACKBONE_ATOMS) for x in t_res + b_res):
            continue

        b_seq_raw = "".join(x["aa"] for x in b_res)
        audit3d = audit_3di3_glycan_and_il7(pdb_p, b_seq_raw, ref_3di3)
        if not audit3d["aligned_ok"] or audit3d["glycan_clashes"] > 2 or audit3d["il7_overlap_atoms"] == 0:
            continue

        b_heavy = np.vstack([np.asarray(x["heavy_atoms"]) for x in b_res])
        hs_hits: list[str] = []
        near_mask = False
        for idx, tres in enumerate(t_res):
            th = np.asarray(tres["heavy_atoms"])
            min_d = float(np.sqrt(np.min(np.sum((th[:, None, :] - b_heavy[None, :, :]) ** 2, axis=2))))
            if min_d <= 5.0 and target_labels[idx] in hotspots:
                hs_hits.append(target_labels[idx])
            if idx in masked_indices and min_d < 8.0:
                near_mask = True
        if near_mask or (hotspots and not hs_hits):
            continue

        mpnn_entry = {
            "name": did,
            "seq_chain_A": wt_target_seq,
            "coords_chain_A": {f"{an}_chain_A": [x["bb"][an].tolist() for x in t_res] for an in BACKBONE_ATOMS},
            "seq_chain_B": b_seq_raw,
            "coords_chain_B": {f"{an}_chain_B": [x["bb"][an].tolist() for x in b_res] for an in BACKBONE_ATOMS},
        }
        backbone_ca_map[did] = np.array([x["bb"]["CA"] for x in b_res], dtype=float)
        passing_bbs.append(
            {
                "design_id": did,
                "i_pTM": iptm,
                "pLDDT": plddt,
                "i_pAE": ipae,
                "hs_hits": hs_hits,
                "core_anchors": _select_t7_core_anchors(b_res, t_res, target_labels, hotspots),
                "mpnn_entry": mpnn_entry,
            }
        )

    print(
        f"[T7 Relaxed Harvest] Audited {len(traj_rows)} trajectories -> {len(passing_bbs)} passed "
        f"3DI3 N-glycan clearance, IL-7 orthosteric overlap, and hotspot gates.",
        flush=True,
    )
    if not passing_bbs:
        return 0

    mpnn_repo = ensure_protein_mpnn_repo()
    work_dir = run_dir / "harvest_soluble_mpnn_t7"
    work_dir.mkdir(parents=True, exist_ok=True)

    raw_cands: list[dict[str, Any]] = []
    bb_by_id = {b["design_id"]: b for b in passing_bbs}

    for mode_tag, n_seqs, temp, fix_core, seed in [
        ("core_fix", num_seqs_core_fix, 0.12, True, 101),
        ("full", num_seqs_full, 0.15, False, 73),
    ]:
        if n_seqs <= 0:
            continue
        mdir = work_dir / f"mpnn_{mode_tag}"
        mdir.mkdir(parents=True, exist_ok=True)
        parsed_j = mdir / "parsed.jsonl"
        chains_j = mdir / "chains.jsonl"
        fixed_j = mdir / "fixed.jsonl"
        c_dict: dict[str, Any] = {}
        f_dict: dict[str, Any] = {}
        active_ids: list[str] = []

        with parsed_j.open("w", encoding="utf-8") as pf:
            for b in passing_bbs:
                if fix_core and not b["core_anchors"]:
                    continue
                did = b["design_id"]
                active_ids.append(did)
                pf.write(json.dumps(b["mpnn_entry"]) + "\n")
                c_dict[did] = [["B"], ["A"]]
                if fix_core:
                    f_dict[did] = {"A": [], "B": b["core_anchors"]}

        if not active_ids:
            continue
        chains_j.write_text(json.dumps(c_dict) + "\n", encoding="utf-8")
        cmd = [
            sys.executable,
            str((mpnn_repo / "protein_mpnn_run.py").resolve()),
            "--jsonl_path", str(parsed_j),
            "--chain_id_jsonl", str(chains_j),
            "--out_folder", str(mdir.resolve()) + "/",
            "--path_to_model_weights", str((mpnn_repo / "soluble_model_weights").resolve()) + "/",
            "--use_soluble_model",
            "--model_name", "v_48_020",
            "--omit_AAs", omit_aas,
            "--num_seq_per_target", str(n_seqs),
            "--batch_size", "1",
            "--sampling_temp", str(temp),
            "--seed", str(seed),
            "--suppress_print", "1",
        ]
        if fix_core:
            fixed_j.write_text(json.dumps(f_dict) + "\n", encoding="utf-8")
            cmd.extend(["--fixed_positions_jsonl", str(fixed_j)])

        env = os.environ.copy()
        env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
        subprocess.run(cmd, env=env, check=True)

        for did in active_ids:
            fa_p = mdir / "seqs" / f"{did}.fa"
            if not fa_p.exists():
                continue
            lines = [ln.strip() for ln in fa_p.read_text(encoding="utf-8").splitlines() if ln.strip()]
            sample_idx = 0
            for i in range(2, len(lines), 2):
                if i + 1 > len(lines):
                    break
                seq = lines[i + 1].split("/")[-1].strip().upper()
                sample_idx += 1
                raw_cands.append(
                    {
                        "design_id": f"{did}_smpnn_{mode_tag}_{sample_idx:02d}",
                        "backbone_id": did,
                        "sequence": seq,
                    }
                )

    # Filter by expressibility + <=70% novelty + deduplicate by sequence
    seen_seqs: set[str] = set()
    if registry_csv.exists():
        with registry_csv.open("r", newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("sequence"):
                    seen_seqs.add(r["sequence"].strip().upper())

    filtered: list[dict[str, Any]] = []
    for c in raw_cands:
        seq = c["sequence"]
        if seq in seen_seqs:
            continue
        exp = compute_expressibility_metrics_t7(seq, min_len=min_len_h, max_len=max_len_h)
        if exp["expressibility_verdict"] != "PASS":
            continue
        novel_ok, max_id, _ = check_sequence_novelty_t7(seq, known_fasta)
        if not novel_ok:
            continue
        seen_seqs.add(seq)
        c["exp"] = exp
        c["max_known_id"] = max_id
        filtered.append(c)

    if not filtered:
        print("[T7 Relaxed Harvest] 0 new SolubleMPNN sequences passed expressibility + novelty gates.", flush=True)
        return 0

    # Optional unbound monomer scRMSD gate via ESMFold on CUDA
    if not skip_esmfold:
        try:
            import torch

            if torch.cuda.is_available():
                chunk = detect_gpu_chunk()
                tok, model = tp.load_hf("cuda", chunk, fp16=True)
                by_len: dict[int, list[dict[str, Any]]] = {}
                for c in filtered:
                    by_len.setdefault(len(c["sequence"]), []).append(c)
                for slen, grp in sorted(by_len.items()):
                    for st in range(0, len(grp), 8):
                        sub = grp[st : st + 8]
                        outs = tp.infer_hf_batch(
                            model, tok, [x["sequence"] for x in sub], [np.arange(slen, dtype=np.int64) for _ in sub], "cuda", recycles=3
                        )
                        for item, (_, pl_arr, pdb_str) in zip(sub, outs):
                            cas = [
                                [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                                for ln in pdb_str.splitlines()
                                if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA"
                            ]
                            _, _, _, sc_rmsd = kabsch_fit(np.array(cas, dtype=float), backbone_ca_map[item["backbone_id"]])
                            item["sc_rmsd"] = round(sc_rmsd, 3)
                            item["mono_plddt"] = round(float(np.mean(pl_arr)), 2)
                del model
                torch.cuda.empty_cache()
                filtered = [x for x in filtered if x.get("sc_rmsd", 0.0) <= 2.2 and x.get("mono_plddt", 99.0) >= 75.0]
        except Exception as exc:
            print(f"[T7 Relaxed Harvest] Monomer scRMSD check skipped: {exc}", flush=True)

    # Cap per backbone to preserve structural diversity (SOP §D-2)
    kept_per_bb: dict[str, int] = {}
    selected: list[dict[str, Any]] = []
    filtered.sort(key=lambda x: (x.get("sc_rmsd", 1.0), -bb_by_id[x["backbone_id"]]["i_pTM"]))
    for c in filtered:
        bid = c["backbone_id"]
        if kept_per_bb.get(bid, 0) >= max_per_backbone:
            continue
        kept_per_bb[bid] = kept_per_bb.get(bid, 0) + 1
        selected.append(c)

    if not selected:
        return 0

    official_seq = (job.get("target_sequence") or OFFICIAL_IL7RA_ECD_21_239).strip().upper()
    ipsae_map = evaluate_candidates_with_ipsae(
        [(x["design_id"], x["sequence"]) for x in selected],
        target_seq=official_seq,
        eval_dir=run_dir / "ipsae_eval_harvest",
        dry_run=skip_esmfold,
    )

    lengths = job.get("lengths", [55, 85])
    new_records: list[dict[str, Any]] = []
    for c in selected:
        bb = bb_by_id[c["backbone_id"]]
        exp = c["exp"]
        ev = ipsae_map.get(c["design_id"], {})
        rec = {
            "design_id": c["design_id"],
            "job_id": job.get("job_id", ""),
            "target_id": job.get("target_id", "target7_IL7Ra"),
            "molecule_type": job.get("molecule_type", "minibinder"),
            "submitted_by": job.get("submitted_by", ""),
            "generator": "BindCraft+SolubleMPNN_T7Harvest",
            "seed": "101",
            "target_file": target_pdb.name,
            "target_sha256": target_sha,
            "length_spec": f"{lengths[0]}-{lengths[1]}" if len(lengths) == 2 else str(lengths),
            "hotspot_string": job.get("hotspots", ""),
            "inverse_folding_model": "SolubleMPNN",
            "is_solubleMPNN": "TRUE",
            "sequence": c["sequence"],
            "sequence_sha256": hashlib.sha256(c["sequence"].encode("utf-8")).hexdigest(),
            "length_aa": exp["length_aa"],
            "fold_pLDDT": round(bb["pLDDT"], 4),
            "fold_ipTM": round(bb["i_pTM"], 4),
            "fold_iPAE": round(bb["i_pAE"], 4),
            "esm_plddt": ev.get("esm_plddt", ""),
            "ipsae10_min": ev.get("ipsae10_min", ""),
            "ipsae10_max": ev.get("ipsae10_max", ""),
            "pdockq": ev.get("pdockq", ""),
            "min_heavy_dist": ev.get("min_heavy_dist", ""),
            "E_ratio": exp["E_ratio"],
            "K_ratio": exp["K_ratio"],
            "net_charge_pH6": exp["net_charge_pH6"],
            "net_charge_pH74": exp["net_charge_pH74"],
            "gravy": exp["gravy"],
            "cys_count": exp["cys_count"],
            "cys_parity_ok": exp["cys_parity_ok"],
            "longest_hydrophobic_run": exp["longest_hydrophobic_run"],
            "expressibility_verdict": "PASS",
            "created_at": dt.datetime.now().isoformat(timespec="seconds"),
            "submitted": "FALSE",
        }
        with registry_csv.open("a", newline="", encoding="utf-8") as out_f:
            csv.DictWriter(out_f, fieldnames=REGISTRY_HEADERS).writerow(rec)
        new_records.append(rec)

    af3_files = export_af3_and_submission_preview(new_records, job, run_dir, suffix_tag="_harvest")
    if af3_files:
        print(
            f"✅ [T7 AF3 Export] 已匯出 {len(af3_files)} 個 1:1 (Minibinder + 219aa IL-7Rα) AlphaFold 3 JSON -> "
            + ", ".join(str(p) for p in af3_files),
            flush=True,
        )
    return len(new_records)


def harvest_bindcraft_candidates_t7(
    job: dict[str, Any],
    target_pdb: pathlib.Path,
    run_dir: pathlib.Path,
    registry_csv: pathlib.Path,
    dry_run: bool = False,
) -> int:
    """Harvest BindCraft designs, enforce 40-100aa + SolubleMPNN + <=70% novelty + 3DI3 glycan/IL-7 + ipSAE gates."""
    stats_files = [
        run_dir / "final_design_stats.csv",
        run_dir / "mpnn_design_stats.csv",
    ]
    target_sha = hashlib.sha256(target_pdb.read_bytes()).hexdigest()
    target_seq = (job.get("target_sequence") or OFFICIAL_IL7RA_ECD_21_239).strip().upper()
    hard_lim = job.get("lengths_hard_limit", [40, 100])
    min_len_h, max_len_h = int(hard_lim[0]), int(hard_lim[1])

    ref_3di3 = (REPO_ROOT / job.get("glycan_reference_pdb_relpath", "targets/target7_IL7Ra/evidence/3DI3.pdb")).resolve()
    known_fasta = (
        REPO_ROOT / job.get("known_binders_fasta_relpath", "targets/target7_IL7Ra/targets/known_binders_IL7Ra_SEQRES.fasta")
    ).resolve()

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
            for idx, row in enumerate(csv.DictReader(f), start=1):
                seq = (row.get("Sequence") or row.get("sequence") or "").strip().upper()
                if not seq:
                    continue
                design_id = row.get("Design") or row.get("design_id") or f"{job['job_id']}_{idx:03d}"
                if design_id in existing_ids:
                    continue
                existing_ids.add(design_id)
                raw_rows.append((design_id, seq, row))

    min_desired = int(job.get("number_of_final_designs", 10))
    if not raw_rows:
        traj_csv = run_dir / "trajectory_stats.csv"
        relaxed_dir = run_dir / "Trajectory" / "Relaxed"
        if not dry_run and traj_csv.exists() and relaxed_dir.exists() and any(relaxed_dir.glob("*.pdb")):
            print(
                "[T7 Fallback Harvest] Default MPNN CSV yielded 0 rows; triggering SolubleMPNN + 3DI3 glycan/IL-7 harvest...",
                flush=True,
            )
            return harvest_relaxed_t7_trajectories(run_dir=run_dir, job=job, registry_csv=registry_csv)
        return 0

    ipsae_map = evaluate_candidates_with_ipsae(
        [(cid, seq) for cid, seq, _ in raw_rows],
        target_seq=target_seq,
        eval_dir=run_dir / "ipsae_eval",
        dry_run=dry_run,
    )

    added = 0
    new_records: list[dict[str, Any]] = []
    for design_id, seq, row in raw_rows:
        seq_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        exp = compute_expressibility_metrics_t7(seq, min_len=min_len_h, max_len=max_len_h)
        verdict = str(exp.get("expressibility_verdict", "PASS"))

        # Gate 1: SOP §C-1 Novelty gate (<= 70% identity vs 7OPB_D, IL-7, TSLP)
        novel_ok, max_id, ref_hit = check_sequence_novelty_t7(seq, known_fasta)
        if not novel_ok:
            tag = f"NOVELTY_FAIL_{ref_hit}_{max_id:.1f}PCT"
            verdict = f"REJECT:{tag}" if verdict == "PASS" else f"{verdict},{tag}"

        # Gate 2: 3DI3 Crystal N-glycan clash & IL-7 orthosteric overlap check
        cmp_pdb = find_candidate_complex_pdb(run_dir, design_id)
        if cmp_pdb is not None and ref_3di3.exists():
            aud3d = audit_3di3_glycan_and_il7(cmp_pdb, seq, ref_3di3)
            if aud3d["aligned_ok"]:
                if aud3d["glycan_clashes"] > 2:
                    tag = f"GLYCAN_CLASH_{aud3d['glycan_clashes']}"
                    verdict = f"REJECT:{tag}" if verdict == "PASS" else f"{verdict},{tag}"
                if aud3d["il7_overlap_atoms"] == 0:
                    tag = "NON_COMPETITIVE_IL7_0"
                    verdict = f"REJECT:{tag}" if verdict == "PASS" else f"{verdict},{tag}"

        ev = ipsae_map.get(design_id, {})
        lengths = job.get("lengths", [55, 85])
        rec = {
            "design_id": design_id,
            "job_id": job.get("job_id", ""),
            "target_id": job.get("target_id", "target7_IL7Ra"),
            "molecule_type": job.get("molecule_type", "minibinder"),
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
            "fold_ipTM": row.get("i_pTM") or row.get("Average_i_pTM") or row.get("ipTM") or "",
            "fold_iPAE": row.get("i_pAE") or row.get("Average_i_pAE") or row.get("iPAE") or "",
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
            "expressibility_verdict": verdict,
            "created_at": dt.datetime.now().isoformat(timespec="seconds"),
            "submitted": "FALSE",
        }
        with registry_csv.open("a", newline="", encoding="utf-8") as out_f:
            csv.DictWriter(out_f, fieldnames=REGISTRY_HEADERS).writerow(rec)
        new_records.append(rec)
        added += 1

    # Top up with relaxed-trajectory SolubleMPNN harvest if default MPNN accepted fewer than half of desired designs
    n_pass = sum(1 for r in new_records if r["expressibility_verdict"] == "PASS")
    if not dry_run and n_pass < max(4, min_desired // 2):
        traj_csv = run_dir / "trajectory_stats.csv"
        relaxed_dir = run_dir / "Trajectory" / "Relaxed"
        if traj_csv.exists() and relaxed_dir.exists() and any(relaxed_dir.glob("*.pdb")):
            try:
                extra = harvest_relaxed_t7_trajectories(run_dir=run_dir, job=job, registry_csv=registry_csv)
                added += extra
            except Exception as exc:
                print(f"[T7 Top-Up Harvest Warning] {exc}", flush=True)

    af3_files = export_af3_and_submission_preview(new_records, job, run_dir)
    if af3_files:
        print(
            f"✅ [AF3 Export] 已自動匯出 {len(af3_files)} 個 Target 7 (1:1 IL-7Rα) AlphaFold 3 JSON -> "
            + ", ".join(str(p) for p in af3_files),
            flush=True,
        )
    return added


def process_t7_job(
    job_path: pathlib.Path,
    drive_root: pathlib.Path,
    queue_subdir: str = DEFAULT_T7_QUEUE_SUBDIR,
    registry_filename: str = DEFAULT_T7_REGISTRY_FILENAME,
    dry_run: bool = False,
) -> bool:
    t_start = time.time()
    dirs = ensure_t7_workspace(drive_root, queue_subdir=queue_subdir, registry_filename=registry_filename)
    job = json.loads(job_path.read_text(encoding="utf-8"))
    job_id = job.get("job_id", job_path.stem)

    if str(job.get("target_id", "")) != "target7_IL7Ra":
        print(f"[SKIP] Refusing non-Target-7 job {job_path.name} (target_id={job.get('target_id')}) in isolated T7 worker.", flush=True)
        return False

    running_path = dirs["running"] / job_path.name
    is_resumed = job_path.resolve() == running_path.resolve()
    if job_path.parent == dirs["pending"]:
        shutil.move(str(job_path), str(running_path))
    elif not is_resumed:
        shutil.copy2(str(job_path), str(running_path))

    run_dir = dirs["runs"] / job_id
    run_dir.mkdir(parents=True, exist_ok=True)
    log_file = run_dir / "worker_t7.log"

    def log(msg: str) -> None:
        stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{stamp}] {msg}"
        print(line, flush=True)
        with log_file.open("a", encoding="utf-8") as lf:
            lf.write(line + "\n")

    gpu_info = inspect_gpu_environment()
    action_str = "Resuming interrupted T7 job" if is_resumed else "Starting isolated T7 job"
    log(
        f"{action_str} {job_id} | Queue={queue_subdir} | Registry={registry_filename} | "
        f"GPU={gpu_info['gpu_name']} ({gpu_info['vram_gb']} GB, {gpu_info['sm_arch']})"
    )
    if gpu_info["warning"]:
        log(f"[GPU WARNING] {gpu_info['warning']}")

    ok = False
    interrupted = False
    try:
        if not dry_run and not gpu_info["bindcraft_compatible"]:
            log(f"ERROR: Incompatible GPU {gpu_info['gpu_name']} ({gpu_info['sm_arch']}). Switch Colab runtime to A100!")
            return False

        target_rel = job.get("target_pdb_relpath", "")
        target_pdb = (REPO_ROOT / target_rel).resolve()
        if not target_pdb.exists():
            log(f"ERROR: Target PDB not found at {target_pdb}")
            return False

        pre_ok, val_out = run_pre_validation(
            target_pdb=target_pdb,
            hotspots=job.get("hotspots", ""),
            expected_sha=job.get("expected_sha256", ""),
        )
        log("Pre-validation (SOP §B-1 --strict-sidechains) output:\n" + val_out.strip())
        if not pre_ok:
            log("ERROR: Pre-validation failed. Aborting job.")
            return False

        settings_path, advanced_path = write_bindcraft_settings(job, target_pdb, run_dir)
        log(f"Wrote BindCraft settings to {settings_path} and {advanced_path}")

        bc_returncode = 0
        if dry_run:
            log("Dry-run mode enabled: skipping GPU BindCraft execution.")
        else:
            bindcraft_script = pathlib.Path("/content/BindCraft/bindcraft.py")
            bindcraft_params = pathlib.Path("/content/BindCraft/params")
            bc_env_python = pathlib.Path("/content/miniforge3/envs/BindCraft/bin/python")
            misplaced_params = pathlib.Path("/content/params")
            if (
                bindcraft_script.exists()
                and not (bindcraft_params / "params_model_1_ptm.npz").exists()
                and (misplaced_params / "params_model_1_ptm.npz").exists()
            ):
                if bindcraft_params.exists():
                    shutil.rmtree(bindcraft_params)
                shutil.move(str(misplaced_params), str(bindcraft_params))

            for _bin_rel in ("functions/dssp", "functions/DAlphaBall.gcc"):
                _bin_p = pathlib.Path("/content/BindCraft") / _bin_rel
                if _bin_p.exists():
                    try:
                        _bin_p.chmod(0o755)
                    except Exception:
                        pass

            if not bindcraft_script.exists() or not (bindcraft_params / "params_model_1_ptm.npz").exists():
                interrupted = True
                log(
                    "❌ [ERROR] 尚未載入 BindCraft + SolubleMPNN + AF2 權重環境 (/content/BindCraft)！\n"
                    "   👉 請先執行 T7 Notebook『第 2 格（45 秒快取還原）』，再執行『第 4 格』！"
                )
                return False

            log(f"Launching BindCraft for Target 7 (IL-7Rα) in {run_dir} ...")
            bc_env = os.environ.copy()
            bc_env["PATH"] = "/content/miniforge3/envs/BindCraft/bin:/content/miniforge3/bin:" + bc_env.get("PATH", "")
            bc_env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.88"
            bc_env["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true"
            bc_cmd = (
                [
                    str(bc_env_python),
                    "-u",
                    str(bindcraft_script),
                    "--settings",
                    str(settings_path),
                    "--filters",
                    "/content/BindCraft/settings_filters/default_filters.json",
                    "--advanced",
                    str(advanced_path),
                ]
                if bc_env_python.exists()
                else [
                    "conda", "run", "--no-capture-output", "-n", "BindCraft", "python", "-u",
                    str(bindcraft_script),
                    "--settings", str(settings_path),
                    "--filters", "/content/BindCraft/settings_filters/default_filters.json",
                    "--advanced", str(advanced_path),
                ]
            )
            proc = subprocess.Popen(
                bc_cmd,
                cwd="/content/BindCraft",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                env=bc_env,
            )
            assert proc.stdout is not None
            for line in proc.stdout:
                clean = line.rstrip("\r\n")
                if clean:
                    log(clean)
            proc.wait()
            bc_returncode = int(proc.returncode or 0)

        harvested = harvest_bindcraft_candidates_t7(
            job, target_pdb, run_dir, dirs["registry_csv"], dry_run=dry_run
        )
        log(f"Harvested & evaluated {harvested} Target 7 candidate sequences into {dirs['registry_csv'].name}.")
        if not dry_run and bc_returncode != 0 and harvested == 0:
            log(f"❌ [ERROR] BindCraft exited with code {bc_returncode} and produced 0 designs.")
            ok = False
            return False

        ok = True
        return True
    except KeyboardInterrupt:
        interrupted = True
        log(f"Job {job_id} interrupted by user. Keeping in {queue_subdir}/running/ for auto-resume.")
        raise
    except Exception as exc:
        log(f"FATAL EXCEPTION in T7 job {job_id}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
        ok = False
        return False
    finally:
        elapsed_sec = round(time.time() - t_start, 1)
        est_cu_used = round((elapsed_sec / 3600.0) * float(gpu_info["est_cu_per_hour"]), 3)
        manifest = {
            "job_id": job_id,
            "target_id": "target7_IL7Ra",
            "queue_subdir": queue_subdir,
            "registry_filename": registry_filename,
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
            log(f"Job {job_id} finalized (ok={ok}, elapsed={elapsed_sec}s, est_CU={est_cu_used}) -> {dest_dir / job_path.name}")
        sync_filesystem()


def _find_t7_job_json(dirs: dict[str, pathlib.Path], job_id: str) -> dict[str, Any]:
    for folder in (
        dirs["running"],
        dirs["pending"],
        dirs["completed"],
        dirs["failed"],
        REPO_ROOT / "jobs" / "templates",
    ):
        if not folder.exists():
            continue
        for jp in sorted(folder.glob("*target7*.json")):
            try:
                d = json.loads(jp.read_text(encoding="utf-8"))
                if d.get("job_id") == job_id:
                    return d
            except Exception:
                continue
    # Fallback to primary T7 template with overridden job_id (e.g. auto-incremented _002, _003)
    default_t7 = REPO_ROOT / "jobs" / "templates" / "job_target7_IL7Ra_minibinder.json"
    d = json.loads(default_t7.read_text(encoding="utf-8"))
    d["job_id"] = job_id
    return d


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Changjiang Cup 2026 — Isolated Target 7 (IL-7Rα) Colab Worker")
    parser.add_argument("--drive-root", default=SHARED_DRIVE_FOLDER_URL)
    parser.add_argument("--queue-subdir", default=DEFAULT_T7_QUEUE_SUBDIR)
    parser.add_argument("--registry-filename", default=DEFAULT_T7_REGISTRY_FILENAME)
    parser.add_argument("--job-file", help="Run a specific Target 7 job JSON immediately")
    parser.add_argument("--harvest-job-id", help="Cell-5 standalone SolubleMPNN + ESMFold harvest on an existing T7 run_dir")
    parser.add_argument("--min-iptm", type=float, default=0.72)
    parser.add_argument("--num-seqs-core-fix", type=int, default=6)
    parser.add_argument("--num-seqs-full", type=int, default=6)
    parser.add_argument("--omit-aas", default="CXM")
    parser.add_argument("--max-per-backbone", type=int, default=2)
    parser.add_argument("--skip-esmfold", action="store_true")
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--unassign-on-exit", action="store_true")
    parser.add_argument("--diagnose", action="store_true")
    parser.add_argument("--poll-interval", type=int, default=15)
    args = parser.parse_args(argv)

    dirs = ensure_t7_workspace(
        pathlib.Path(args.drive_root),
        queue_subdir=args.queue_subdir,
        registry_filename=args.registry_filename,
    )
    drive_root = dirs["root"]

    if args.diagnose:
        print_t7_diagnostic(drive_root, queue_subdir=args.queue_subdir, registry_filename=args.registry_filename)
        return 2

    if args.harvest_job_id:
        jid = args.harvest_job_id.strip()
        run_dir = dirs["runs"] / jid
        job = _find_t7_job_json(dirs, jid)
        added = harvest_relaxed_t7_trajectories(
            run_dir=run_dir,
            job=job,
            registry_csv=dirs["registry_csv"],
            min_iptm=args.min_iptm,
            num_seqs_core_fix=args.num_seqs_core_fix,
            num_seqs_full=args.num_seqs_full,
            omit_aas=args.omit_aas,
            max_per_backbone=args.max_per_backbone,
            skip_esmfold=args.skip_esmfold,
        )
        print(f"🎉 [T7 Standalone Harvest] 新增 {added} 條 PASS 候選序列至 {dirs['registry_csv'].name}。", flush=True)
        return 0

    jobs_completed = 0
    any_failed = False
    try:
        if args.job_file:
            ok = process_t7_job(
                pathlib.Path(args.job_file),
                drive_root,
                queue_subdir=args.queue_subdir,
                registry_filename=args.registry_filename,
                dry_run=args.dry_run,
            )
            if ok:
                jobs_completed += 1
            else:
                any_failed = True
        else:
            print(
                f"Target 7 Isolated Worker watching: {dirs['pending']} "
                f"(auto-resuming from {dirs['running']}, T5 queue untouched)",
                flush=True,
            )
            auto_recovered_once = False
            while True:
                resumed_jobs = sorted(dirs["running"].glob("*.json"))
                if resumed_jobs:
                    print(f"[AUTO-RESUME] 檢測到 {len(resumed_jobs)} 個中斷的 T7 任務，正在從斷點續跑...", flush=True)
                pending_jobs = sorted(dirs["pending"].glob("*.json"))
                all_queue = resumed_jobs + pending_jobs

                if not all_queue and not auto_recovered_once:
                    auto_recovered_once = True
                    for t7_name in (
                        "job_target7_IL7Ra_minibinder.json",
                        "job_target7_IL7Ra_holo7OPB.json",
                        "job_target7_IL7Ra_apo3UP1.json",
                        "job_target7_IL7Ra_backMaskK.json",
                    ):
                        for src_folder in (dirs["completed"], dirs["failed"]):
                            cand_p = src_folder / t7_name
                            if cand_p.exists():
                                try:
                                    jdata = json.loads(cand_p.read_text(encoding="utf-8"))
                                    jid = jdata.get("job_id", cand_p.stem)
                                    rdir = dirs["runs"] / jid
                                    has_csv = (rdir / "final_design_stats.csv").exists() or (rdir / "trajectory_stats.csv").exists()
                                    has_pdb = any(rdir.rglob("*.pdb")) if rdir.exists() else False
                                    if not has_csv and not has_pdb:
                                        rec_p = dirs["running"] / t7_name
                                        shutil.move(str(cand_p), str(rec_p))
                                        print(f"[AUTO-RECOVER] 已自動恢復 0 產出的 T7 任務 {t7_name} 至 {rec_p}！", flush=True)
                                        all_queue.append(rec_p)
                                except Exception:
                                    pass

                if not all_queue and args.once:
                    print(
                        f"[提示] T7 獨立隊列為空 ({dirs['pending']} 與 {dirs['running']} 均無待執行任務)。"
                        "請先執行『第 3 格』投遞 Target 7 任務！",
                        flush=True,
                    )
                for jp in all_queue:
                    ok = process_t7_job(
                        jp,
                        drive_root,
                        queue_subdir=args.queue_subdir,
                        registry_filename=args.registry_filename,
                        dry_run=args.dry_run,
                    )
                    if ok:
                        jobs_completed += 1
                    else:
                        any_failed = True
                if args.once:
                    break
                time.sleep(args.poll_interval)
    except KeyboardInterrupt:
        any_failed = True
        print(f"\n[中斷] 收到停止信號，T7 任務已保留在 {args.queue_subdir}/running/ 以便下次續跑。", flush=True)
    finally:
        ready_to_unassign = not args.dry_run and not any_failed and jobs_completed > 0
        if args.unassign_on_exit:
            if ready_to_unassign:
                flush_drive_and_unassign_gpu(drive_root, unassign=True, countdown_sec=15)
            else:
                sync_filesystem()
                print(
                    f"\n[SAFEGUARD] 未觸發自動斷開 GPU (completed={jobs_completed}, failed_or_stopped={any_failed}, dry_run={args.dry_run})。",
                    flush=True,
                )

    if any_failed:
        return 1
    if not args.dry_run and jobs_completed > 0:
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
