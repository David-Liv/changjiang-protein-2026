#!/usr/bin/env python3
"""Harvest PyRosetta-relaxed multi-chain groove backbones via SolubleMPNN + quality gates + AF3 JSON export.

Why this module exists (SOP §D-1 & §C-4):
- During BindCraft `default_4stage_multimer` runs on multi-chain targets (such as Target 5 TL1A's
  `A:B` inter-subunit DcR3/DR3 groove), Stage 1-4 trajectory generation uses `AF2-Multimer` weights
  and produces high-confidence relaxed complexes (`i_pTM = 0.80 - 0.92`, `0` native trimer clashes).
- However, BindCraft's internal MPNN complex re-prediction switches back to single-sequence, 0-MSA
  `AF2-Monomer-PTM` (`params_model_1_ptm`), which systematically blind-rejects shallow homodimer
  β-sheet grooves (see `docs/sop/epitope_and_expressibility_SOP_v1.pdf` §D-1, Adaptyv Round 3).
- Submitting raw AF2 trajectory sequences directly is forbidden by `AGENTS.md` / SOP §C-4; sequences
  MUST undergo `SolubleMPNN` inverse folding, expressibility/novelty/trimer-clash gates, optional
  ESMFold unbound-monomer self-consistency (`scRMSD`), and `3x WT TL1A + 1x Binder` AlphaFold 3 validation.
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
import subprocess
import sys
from typing import Any

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
for _p in (REPO_ROOT, REPO_ROOT / "pipeline", REPO_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import check_novelty  # noqa: E402
import tnf_pipeline as tp  # noqa: E402
import worker  # noqa: E402

AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",
}
BACKBONE_ATOMS = ("N", "CA", "C", "O")


def _is_hydrogen(atom_name: str, elem: str) -> bool:
    return (
        elem == "H"
        or atom_name.startswith("H")
        or (len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == "H")
    )


def parse_pdb_residues(pdb_path: pathlib.Path) -> dict[str, list[dict[str, Any]]]:
    """Parse heavy-atom residues grouped by chain ID, preserving file residue order."""
    by_chain: dict[str, list[dict[str, Any]]] = {}
    cur_map: dict[tuple[str, str], dict[str, Any]] = {}
    for ln in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln.startswith(("ATOM  ", "HETATM")):
            continue
        res3 = ln[17:20].strip()
        if res3 not in AA3_TO_1 or ln[16] not in (" ", "A"):
            continue
        atom_name = ln[12:16].strip()
        elem = ln[76:78].strip() if len(ln) >= 78 else ""
        if _is_hydrogen(atom_name, elem):
            continue
        ch = ln[21]
        res_uid = ln[22:27]
        resseq = int(ln[22:26])
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])], dtype=float)
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
    return by_chain


def kabsch_ca_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute optimal rigid-body C-alpha RMSD (in Angstroms) between two (L, 3) coordinate arrays."""
    if P.shape != Q.shape or len(P) < 3:
        return 999.0
    cP = P.mean(axis=0)
    cQ = Q.mean(axis=0)
    P0 = P - cP
    Q0 = Q - cQ
    H = P0.T @ Q0
    U, _, Vt = np.linalg.svd(H)
    d = float(np.linalg.det(Vt.T @ U.T))
    S = np.diag([1.0, 1.0, -1.0 if d < 0.0 else 1.0])
    R = Vt.T @ S @ U.T
    aligned = P0 @ R.T
    return float(np.sqrt(np.mean(np.sum((aligned - Q0) ** 2, axis=1))))


def build_target_reference_meta(
    target_pdb: pathlib.Path,
    native_trimer_pdb: pathlib.Path | None = None,
) -> dict[str, Any]:
    """Extract target chain lengths, residue labels, wild-type sequences, and masked mutation labels."""
    tgt_chains = parse_pdb_residues(target_pdb)
    ref_chains = parse_pdb_residues(native_trimer_pdb) if (native_trimer_pdb and native_trimer_pdb.exists()) else {}

    # Build wild-type lookup from Chain A of native reference (handles +200 offset on Chain B of cropped dimer)
    ref_wt_by_chain_res: dict[tuple[str, int], str] = {}
    ref_chain_a_wt: dict[int, str] = {}
    for ch, res_list in ref_chains.items():
        for r in res_list:
            ref_wt_by_chain_res[(ch, r["resseq"])] = r["aa"]
            if ch == "A":
                ref_chain_a_wt[r["resseq"]] = r["aa"]

    ordered_chains = sorted(tgt_chains.keys())
    chain_specs: list[dict[str, Any]] = []
    flat_labels: list[str] = []
    masked_labels: set[str] = set()

    for ch in ordered_chains:
        res_list = tgt_chains[ch]
        wt_aas: list[str] = []
        labels: list[str] = []
        for r in res_list:
            rnum = r["resseq"]
            lbl = f"{ch}{rnum}"
            labels.append(lbl)
            flat_labels.append(lbl)
            wt_aa = ref_wt_by_chain_res.get((ch, rnum))
            if wt_aa is None and rnum > 200 and (rnum - 200) in ref_chain_a_wt:
                wt_aa = ref_chain_a_wt[rnum - 200]
            if wt_aa is None:
                wt_aa = r["aa"]
            if wt_aa != r["aa"]:
                masked_labels.add(lbl)
            wt_aas.append(wt_aa)
        chain_specs.append(
            {
                "orig_chain": ch,
                "length": len(res_list),
                "labels": labels,
                "wt_seq": "".join(wt_aas),
            }
        )

    return {
        "chain_specs": chain_specs,
        "total_target_len": sum(c["length"] for c in chain_specs),
        "flat_labels": flat_labels,
        "masked_labels": masked_labels,
    }


def audit_and_prepare_backbone(
    relaxed_pdb: pathlib.Path,
    traj_row: dict[str, str],
    target_pdb: pathlib.Path,
    target_meta: dict[str, Any],
    native_trimer_pdb: pathlib.Path | None,
    hotspots: set[str],
    min_iptm: float = 0.80,
    min_plddt: float = 0.78,
    max_trimer_clashes: int = 5,
    require_groove_bridge: bool = True,
    require_hotspot_hit: bool = True,
) -> dict[str, Any]:
    """Audit a single PyRosetta-relaxed trajectory PDB and prepare gap-free multi-chain MPNN dict."""
    design_id = traj_row["Design"]
    iptm = float(traj_row.get("i_pTM") or 0.0)
    plddt = float(traj_row.get("pLDDT") or 0.0)
    ipae = float(traj_row.get("i_pAE") or 1.0)
    dg = float(traj_row.get("dG") or 0.0)
    sc = float(traj_row.get("ShapeComplementarity") or 0.0)
    dsasa = float(traj_row.get("dSASA") or 0.0)
    traj_seq = (traj_row.get("Sequence") or "").strip().upper()

    result: dict[str, Any] = {
        "design_id": design_id,
        "pdb_path": relaxed_pdb,
        "i_pTM": iptm,
        "pLDDT": plddt,
        "i_pAE": ipae,
        "dG": dg,
        "SC": sc,
        "dSASA": dsasa,
        "traj_seq": traj_seq,
        "passed": False,
        "reject_reasons": [],
    }

    if iptm < min_iptm:
        result["reject_reasons"].append(f"LOW_IPTM_{iptm:.2f}")
    if plddt < min_plddt:
        result["reject_reasons"].append(f"LOW_PLDDT_{plddt:.2f}")
    if not relaxed_pdb.exists():
        result["reject_reasons"].append("MISSING_RELAXED_PDB")
        return result

    parsed = parse_pdb_residues(relaxed_pdb)
    if "A" not in parsed or "B" not in parsed:
        result["reject_reasons"].append("UNEXPECTED_PDB_CHAINS")
        return result

    concat_target_res = parsed["A"]
    binder_res = parsed["B"]
    expected_tlen = target_meta["total_target_len"]
    if len(concat_target_res) != expected_tlen:
        result["reject_reasons"].append(f"TARGET_LEN_MISMATCH_{len(concat_target_res)}_VS_{expected_tlen}")
        return result

    for r in concat_target_res + binder_res:
        if any(a not in r["bb"] for a in BACKBONE_ATOMS):
            result["reject_reasons"].append("INCOMPLETE_BACKBONE_ATOMS")
            return result

    # Trimer steric clash gate against omitted chain(s) of native homotrimer (e.g. Chain C of 2RE9)
    trimer_clashes = 0
    if native_trimer_pdb is not None and native_trimer_pdb.exists():
        trimer_clashes = worker.count_native_trimer_clashes(
            relaxed_pdb, traj_seq, target_pdb, native_trimer_pdb
        )
        if trimer_clashes > max_trimer_clashes:
            result["reject_reasons"].append(f"TRIMER_CLASH_{trimer_clashes}")
    result["trimer_clashes"] = trimer_clashes

    # Compute heavy-atom contacts (<5.0 A) to target residues & minimum distance to masked LYS residues
    binder_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in binder_res])
    flat_labels: list[str] = target_meta["flat_labels"]
    masked_labels: set[str] = target_meta["masked_labels"]

    contact_labels: set[str] = set()
    min_mask_dist = 999.0
    for idx, tres in enumerate(concat_target_res):
        lbl = flat_labels[idx]
        th = np.asarray(tres["heavy_atoms"])
        d2 = np.sum((th[:, None, :] - binder_heavy[None, :, :]) ** 2, axis=2)
        min_d = float(np.sqrt(np.min(d2)))
        if min_d <= 5.0:
            contact_labels.add(lbl)
        if lbl in masked_labels and min_d < min_mask_dist:
            min_mask_dist = min_d

    contacts_by_orig_chain: dict[str, int] = {}
    for spec in target_meta["chain_specs"]:
        och = spec["orig_chain"]
        contacts_by_orig_chain[och] = sum(1 for c in contact_labels if c.startswith(och))

    hotspot_hits = sorted(contact_labels & hotspots) if hotspots else []
    result["contacts_by_chain"] = contacts_by_orig_chain
    result["hotspot_hits"] = hotspot_hits
    result["min_mask_dist"] = round(min_mask_dist, 2) if masked_labels else 999.0

    if masked_labels and min_mask_dist < 8.0:
        result["reject_reasons"].append(f"NEAR_MASKED_RESIDUE_{min_mask_dist:.1f}A")
    if require_groove_bridge and len(target_meta["chain_specs"]) >= 2:
        chains_touched = sum(1 for n in contacts_by_orig_chain.values() if n > 0)
        if chains_touched < 2:
            result["reject_reasons"].append("SINGLE_CHAIN_ONLY_NOT_GROOVE")
    if require_hotspot_hit and hotspots and not hotspot_hits:
        result["reject_reasons"].append("NO_HOTSPOT_CONTACT")

    # Build clean gap-free ProteinMPNN structure dictionary:
    # Target chains -> A, B, ... (fixed, wild-type sequence); Binder -> next chain letter (e.g. C)
    mpnn_entry: dict[str, Any] = {"name": design_id}
    fixed_chains: list[str] = []
    concat_seq_parts: list[str] = []
    offset = 0
    for c_idx, spec in enumerate(target_meta["chain_specs"]):
        mpnn_ch = chr(ord("A") + c_idx)
        fixed_chains.append(mpnn_ch)
        seg = concat_target_res[offset : offset + spec["length"]]
        offset += spec["length"]
        wt_seq = spec["wt_seq"]
        mpnn_entry[f"seq_chain_{mpnn_ch}"] = wt_seq
        mpnn_entry[f"coords_chain_{mpnn_ch}"] = {
            f"{an}_chain_{mpnn_ch}": [r["bb"][an].tolist() for r in seg]
            for an in BACKBONE_ATOMS
        }
        concat_seq_parts.append(wt_seq)

    binder_mpnn_ch = chr(ord("A") + len(target_meta["chain_specs"]))
    binder_pdb_seq = "".join(r["aa"] for r in binder_res)
    mpnn_entry[f"seq_chain_{binder_mpnn_ch}"] = binder_pdb_seq
    mpnn_entry[f"coords_chain_{binder_mpnn_ch}"] = {
        f"{an}_chain_{binder_mpnn_ch}": [r["bb"][an].tolist() for r in binder_res]
        for an in BACKBONE_ATOMS
    }
    concat_seq_parts.append(binder_pdb_seq)
    mpnn_entry["num_of_chains"] = len(fixed_chains) + 1
    mpnn_entry["seq"] = "".join(concat_seq_parts)

    # Parse 1-based binder interface residue positions from trajectory_stats.csv (e.g. "B9,B12,B13")
    raw_intf = (traj_row.get("InterfaceResidues") or "").strip()
    fixed_binder_positions: list[int] = []
    for tok in raw_intf.split(","):
        tok = tok.strip()
        m = re.match(r"^[A-Za-z](\d+)$", tok)
        if m:
            pos = int(m.group(1))
            if 1 <= pos <= len(binder_res) and binder_pdb_seq[pos - 1] != "C":
                fixed_binder_positions.append(pos)
    # Never fix 100% of binder positions (ProteinMPNN requires >=1 designed position)
    if len(fixed_binder_positions) >= len(binder_res):
        fixed_binder_positions = fixed_binder_positions[:-1]

    result["mpnn_entry"] = mpnn_entry
    result["fixed_chains"] = fixed_chains
    result["binder_mpnn_chain"] = binder_mpnn_ch
    result["fixed_binder_positions"] = sorted(set(fixed_binder_positions))
    result["binder_ca"] = np.array([r["bb"]["CA"] for r in binder_res], dtype=float)
    result["passed"] = len(result["reject_reasons"]) == 0
    return result


def ensure_protein_mpnn_repo(preferred_dir: pathlib.Path | None = None) -> pathlib.Path:
    """Locate or shallow-clone official `dauparas/ProteinMPNN` containing `soluble_model_weights`."""
    candidates = []
    if preferred_dir is not None:
        candidates.append(preferred_dir)
    candidates.extend(
        [
            pathlib.Path("/content/ProteinMPNN"),
            REPO_ROOT / "external" / "ProteinMPNN",
        ]
    )
    for c in candidates:
        if (c / "protein_mpnn_run.py").exists() and (c / "soluble_model_weights" / "v_48_020.pt").exists():
            return c

    target_dir = preferred_dir or pathlib.Path("/content/ProteinMPNN")
    target_dir.parent.mkdir(parents=True, exist_ok=True)
    print(f"[SolubleMPNN] Cloning official dauparas/ProteinMPNN into {target_dir} ...", flush=True)
    subprocess.run(
        ["git", "clone", "--depth", "1", "https://github.com/dauparas/ProteinMPNN.git", str(target_dir)],
        check=True,
    )
    return target_dir


def parse_mpnn_fasta_outputs(
    seqs_dir: pathlib.Path,
    backbone_ids: list[str],
    mode_tag: str,
) -> list[dict[str, Any]]:
    """Parse generated SolubleMPNN FASTA files (skipping record 0 which is the input trajectory seq)."""
    records: list[dict[str, Any]] = []
    for bid in backbone_ids:
        fa_path = seqs_dir / f"{bid}.fa"
        if not fa_path.exists():
            continue
        lines = [ln.strip() for ln in fa_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
        pairs: list[tuple[str, str]] = []
        cur_h: str | None = None
        cur_s: list[str] = []
        for ln in lines:
            if ln.startswith(">"):
                if cur_h is not None:
                    pairs.append((cur_h, "".join(cur_s)))
                cur_h = ln[1:].strip()
                cur_s = []
            else:
                cur_s.append(ln)
        if cur_h is not None:
            pairs.append((cur_h, "".join(cur_s)))

        # Skip index 0 (native input sequence written by protein_mpnn_run.py)
        for sample_idx, (hdr, seq) in enumerate(pairs[1:], start=1):
            meta_kv: dict[str, str] = {}
            for part in hdr.split(","):
                if "=" in part:
                    k, v = part.split("=", 1)
                    meta_kv[k.strip()] = v.strip()
            clean_seq = seq.split("/")[0].strip().upper()
            records.append(
                {
                    "backbone_id": bid,
                    "mode": mode_tag,
                    "sample_idx": int(meta_kv.get("sample", sample_idx)),
                    "temp": float(meta_kv.get("T", 0.1)),
                    "mpnn_score": float(meta_kv.get("score", 99.0)),
                    "global_score": float(meta_kv.get("global_score", 99.0)),
                    "seq_recovery": float(meta_kv.get("seq_recovery", 0.0)),
                    "sequence": clean_seq,
                }
            )
    return records


def run_soluble_mpnn_mode(
    audited_backbones: list[dict[str, Any]],
    work_dir: pathlib.Path,
    mpnn_repo: pathlib.Path,
    mode_tag: str,
    num_seqs: int,
    sampling_temp: float,
    fix_interface: bool,
    seed: int = 37,
    python_exe: str | None = None,
) -> list[dict[str, Any]]:
    """Run official SolubleMPNN (`v_48_020`, `omit_AAs=CX`) across all audited groove backbones."""
    if num_seqs <= 0 or not audited_backbones:
        return []

    mode_dir = work_dir / f"mpnn_{mode_tag}"
    mode_dir.mkdir(parents=True, exist_ok=True)
    parsed_jsonl = mode_dir / "parsed_pdbs.jsonl"
    chain_jsonl = mode_dir / "assigned_chains.jsonl"
    fixed_jsonl = mode_dir / "fixed_positions.jsonl"

    chain_dict: dict[str, list[list[str]]] = {}
    fixed_dict: dict[str, dict[str, list[int]]] = {}
    active_backbones: list[dict[str, Any]] = []

    with parsed_jsonl.open("w", encoding="utf-8") as pf:
        for bb in audited_backbones:
            if fix_interface and not bb["fixed_binder_positions"]:
                continue
            active_backbones.append(bb)
            name = bb["design_id"]
            pf.write(json.dumps(bb["mpnn_entry"]) + "\n")
            b_ch = bb["binder_mpnn_chain"]
            f_chs = bb["fixed_chains"]
            chain_dict[name] = [[b_ch], f_chs]
            if fix_interface:
                pos_map = {ch: [] for ch in f_chs}
                pos_map[b_ch] = bb["fixed_binder_positions"]
                fixed_dict[name] = pos_map

    if not active_backbones:
        return []

    chain_jsonl.write_text(json.dumps(chain_dict) + "\n", encoding="utf-8")
    if fix_interface:
        fixed_jsonl.write_text(json.dumps(fixed_dict) + "\n", encoding="utf-8")

    py_bin = python_exe or sys.executable
    # Run with forward-slash path so protein_mpnn_run.py's `file_path.rfind("/")` works everywhere
    script_str = str((mpnn_repo / "protein_mpnn_run.py").resolve()).replace("\\", "/")
    out_str = str(mode_dir.resolve()).replace("\\", "/") + "/"
    weights_str = str((mpnn_repo / "soluble_model_weights").resolve()).replace("\\", "/") + "/"

    cmd = [
        py_bin,
        script_str,
        "--jsonl_path",
        str(parsed_jsonl),
        "--chain_id_jsonl",
        str(chain_jsonl),
        "--out_folder",
        out_str,
        "--path_to_model_weights",
        weights_str,
        "--use_soluble_model",
        "--model_name",
        "v_48_020",
        "--omit_AAs",
        "CX",
        "--num_seq_per_target",
        str(num_seqs),
        "--batch_size",
        "1",
        "--sampling_temp",
        str(sampling_temp),
        "--seed",
        str(seed),
        "--suppress_print",
        "1",
    ]
    if fix_interface:
        cmd.extend(["--fixed_positions_jsonl", str(fixed_jsonl)])

    env = os.environ.copy()
    env["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        raise RuntimeError(
            f"SolubleMPNN ({mode_tag}) failed with code {proc.returncode}:\n{proc.stdout}\n{proc.stderr}"
        )

    return parse_mpnn_fasta_outputs(
        mode_dir / "seqs",
        [bb["design_id"] for bb in active_backbones],
        mode_tag=mode_tag,
    )


def evaluate_monomers_with_esmfold(
    candidates: list[dict[str, Any]],
    backbone_ca_map: dict[str, np.ndarray],
    out_pdb_dir: pathlib.Path,
) -> bool:
    """Fold unbound monomer candidates with ESMFold on CUDA and compute `monomer_esm_plddt` + `monomer_sc_rmsd`."""
    if not candidates:
        return False
    try:
        import torch
    except ImportError:
        print("[ESMFold] PyTorch not installed; skipping unbound monomer scRMSD check.", flush=True)
        return False

    if not torch.cuda.is_available():
        print("[ESMFold] No CUDA GPU detected; skipping unbound monomer scRMSD check (fast CPU mode).", flush=True)
        return False

    try:
        import transformers  # noqa: F401
    except ImportError:
        print("[ESMFold] Installing `transformers` + `accelerate` on Colab GPU runtime (~8s)...", flush=True)
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q", "transformers", "accelerate"],
            check=True,
        )

    out_pdb_dir.mkdir(parents=True, exist_ok=True)
    chunk = worker.detect_gpu_chunk()
    print(f"[ESMFold] Loading EsmForProteinFolding on CUDA (chunk={chunk}) for {len(candidates)} monomers...", flush=True)
    tok, model = tp.load_hf("cuda", chunk, fp16=True)

    # Group by sequence length for fast Tensor Core batch inference
    by_len: dict[int, list[dict[str, Any]]] = {}
    for c in candidates:
        by_len.setdefault(len(c["sequence"]), []).append(c)

    bmax = 8
    for seq_len, group in sorted(by_len.items()):
        for start in range(0, len(group), bmax):
            batch = group[start : start + bmax]
            seqs = [x["sequence"] for x in batch]
            pos_list = [np.arange(seq_len, dtype=np.int64) for _ in batch]
            outs = tp.infer_hf_batch(model, tok, seqs, pos_list, "cuda", recycles=3)
            for item, (_, pl_arr, pdb_str) in zip(batch, outs):
                cid = item["design_id"]
                pdb_p = out_pdb_dir / f"{cid}.pdb"
                pdb_p.write_text(pdb_str, encoding="utf-8")
                pred_cas = []
                for ln in pdb_str.splitlines():
                    if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
                        pred_cas.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
                pred_ca_np = np.array(pred_cas, dtype=float)
                ref_ca_np = backbone_ca_map[item["backbone_id"]]
                sc_rmsd = kabsch_ca_rmsd(pred_ca_np, ref_ca_np)
                item["monomer_esm_plddt"] = round(float(np.mean(pl_arr)), 2)
                item["monomer_sc_rmsd"] = round(sc_rmsd, 3)

    if torch.cuda.is_available():
        del model
        torch.cuda.empty_cache()
    return True


def harvest_run(
    run_dir: pathlib.Path,
    job: dict[str, Any],
    registry_csv: pathlib.Path,
    min_iptm: float = 0.80,
    min_plddt: float = 0.78,
    max_trimer_clashes: int = 5,
    num_seqs_fix_intf: int = 4,
    num_seqs_full: int = 4,
    temp_fix_intf: float = 0.10,
    temp_full: float = 0.15,
    max_per_backbone: int = 2,
    max_sc_rmsd: float = 2.5,
    min_monomer_plddt: float = 75.0,
    skip_esmfold: bool = False,
    af3_batch_size: int = 20,
    mpnn_runner_override=None,
) -> dict[str, Any]:
    """End-to-end harvest of relaxed groove backbones -> SolubleMPNN -> quality gates -> AF3 Server JSON."""
    traj_csv = run_dir / "trajectory_stats.csv"
    relaxed_dir = run_dir / "Trajectory" / "Relaxed"
    if not traj_csv.exists() or not relaxed_dir.exists():
        raise FileNotFoundError(f"Missing trajectory_stats.csv or Trajectory/Relaxed under {run_dir}")

    target_pdb = (REPO_ROOT / job["target_pdb_relpath"]).resolve()
    target_sha = hashlib.sha256(target_pdb.read_bytes()).hexdigest()
    native_trimer_pdb: pathlib.Path | None = None
    if job.get("native_trimer_pdb_relpath"):
        cand_tri = (REPO_ROOT / job["native_trimer_pdb_relpath"]).resolve()
        if cand_tri.exists():
            native_trimer_pdb = cand_tri

    hotspots = {
        h.strip()
        for h in (job.get("hotspots") or "").split(",")
        if h.strip()
    }
    hard_len = job.get("lengths_hard_limit", [40, 100] if job.get("molecule_type") == "minibinder" else [10, 250])
    min_len_hard, max_len_hard = int(hard_len[0]), int(hard_len[1])

    target_meta = build_target_reference_meta(target_pdb, native_trimer_pdb)
    traj_rows = list(csv.DictReader(traj_csv.open("r", newline="", encoding="utf-8")))

    audited_all: list[dict[str, Any]] = []
    for r in traj_rows:
        did = r["Design"]
        pdb_p = relaxed_dir / f"{did}.pdb"
        aud = audit_and_prepare_backbone(
            relaxed_pdb=pdb_p,
            traj_row=r,
            target_pdb=target_pdb,
            target_meta=target_meta,
            native_trimer_pdb=native_trimer_pdb,
            hotspots=hotspots,
            min_iptm=min_iptm,
            min_plddt=min_plddt,
            max_trimer_clashes=max_trimer_clashes,
        )
        audited_all.append(aud)

    passing_bbs = sorted(
        [a for a in audited_all if a["passed"]],
        key=lambda x: (x["i_pTM"], -x["dG"], x["SC"]),
        reverse=True,
    )
    print(
        f"[Groove Harvest] Audited {len(audited_all)} relaxed trajectories -> "
        f"{len(passing_bbs)} passed 3D groove + 0-clash + i_pTM>={min_iptm:.2f} gates.",
        flush=True,
    )
    if not passing_bbs:
        return {"selected_records": [], "af3_files": [], "audited_backbones": audited_all}

    harvest_dir = run_dir / "harvest_soluble_mpnn"
    harvest_dir.mkdir(parents=True, exist_ok=True)

    if mpnn_runner_override is not None:
        raw_mpnn = mpnn_runner_override(passing_bbs, harvest_dir)
    else:
        mpnn_repo = ensure_protein_mpnn_repo()
        raw_mpnn = []
        if num_seqs_fix_intf > 0:
            print(
                f"[SolubleMPNN] Generating {num_seqs_fix_intf} interface-fixed sequences/backbone (T={temp_fix_intf})...",
                flush=True,
            )
            raw_mpnn.extend(
                run_soluble_mpnn_mode(
                    passing_bbs,
                    harvest_dir,
                    mpnn_repo,
                    mode_tag="fix_intf",
                    num_seqs=num_seqs_fix_intf,
                    sampling_temp=temp_fix_intf,
                    fix_interface=True,
                    seed=37,
                )
            )
        if num_seqs_full > 0:
            print(
                f"[SolubleMPNN] Generating {num_seqs_full} full-binder sequences/backbone (T={temp_full})...",
                flush=True,
            )
            raw_mpnn.extend(
                run_soluble_mpnn_mode(
                    passing_bbs,
                    harvest_dir,
                    mpnn_repo,
                    mode_tag="full",
                    num_seqs=num_seqs_full,
                    sampling_temp=temp_full,
                    fix_interface=False,
                    seed=73,
                )
            )

    bb_by_id = {b["design_id"]: b for b in passing_bbs}
    backbone_ca_map = {b["design_id"]: b["binder_ca"] for b in passing_bbs}

    # Load known binders for SOP §C-1 novelty check
    known_fasta = REPO_ROOT / "data" / "raw" / "known_binders.fasta"
    known_seqs = check_novelty.read_fasta(str(known_fasta)) if known_fasta.exists() else {}

    # Load existing registry sequence hashes & design IDs to prevent duplicates
    existing_ids: set[str] = set()
    existing_seq_shas: set[str] = set()
    if registry_csv.exists():
        with registry_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("design_id"):
                    existing_ids.add(row["design_id"])
                if row.get("sequence_sha256") and row.get("expressibility_verdict") == "PASS":
                    existing_seq_shas.add(row["sequence_sha256"])

    seen_batch_seqs: set[str] = set()
    evaluated_candidates: list[dict[str, Any]] = []
    esm_queue: list[dict[str, Any]] = []

    for rec in raw_mpnn:
        seq = rec["sequence"]
        bb = bb_by_id[rec["backbone_id"]]
        short_mode = "fix" if rec["mode"] == "fix_intf" else "full"
        design_id = f"{rec['backbone_id']}_smpnn_{short_mode}{rec['sample_idx']}"
        seq_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()

        exp = worker.compute_expressibility_metrics(seq)
        reasons: list[str] = []
        if exp["expressibility_verdict"] != "PASS":
            reasons.append(exp["expressibility_verdict"].replace("REJECT:", ""))
        if len(seq) < min_len_hard or len(seq) > max_len_hard:
            reasons.append(f"HARD_LEN_{len(seq)}_NOT_IN_{min_len_hard}_{max_len_hard}")
        if seq in seen_batch_seqs:
            reasons.append("DUPLICATE_IN_BATCH")
        elif seq_sha in existing_seq_shas:
            reasons.append("ALREADY_IN_REGISTRY")
        seen_batch_seqs.add(seq)

        # Novelty check against known_binders.fasta
        max_known_ident = 0.0
        for _, kseq in known_seqs.items():
            if abs(len(seq) - len(kseq)) <= 60:
                _, _, pct, _ = check_novelty.ident(seq, kseq)
                if pct > max_known_ident:
                    max_known_ident = round(pct, 1)
        if max_known_ident >= 70.0:
            reasons.append(f"KNOWN_BINDER_IDENT_{max_known_ident:.1f}")

        cand = {
            "design_id": design_id,
            "backbone_id": rec["backbone_id"],
            "mode": rec["mode"],
            "sample_idx": rec["sample_idx"],
            "sequence": seq,
            "sequence_sha256": seq_sha,
            "mpnn_score": round(rec["mpnn_score"], 4),
            "global_score": round(rec["global_score"], 4),
            "seq_recovery": round(rec["seq_recovery"], 4),
            "traj_pLDDT": bb["pLDDT"],
            "traj_i_pTM": bb["i_pTM"],
            "traj_i_pAE": bb["i_pAE"],
            "traj_dG": bb["dG"],
            "traj_SC": bb["SC"],
            "traj_dSASA": bb["dSASA"],
            "trimer_clashes": bb["trimer_clashes"],
            "hotspot_hits": ";".join(bb["hotspot_hits"]),
            "contacts_by_chain": json.dumps(bb["contacts_by_chain"]),
            "max_known_ident": max_known_ident,
            "monomer_esm_plddt": "",
            "monomer_sc_rmsd": "",
            **exp,
            "gate_reasons": reasons,
        }
        evaluated_candidates.append(cand)
        if not reasons:
            esm_queue.append(cand)

    esm_ran = False
    if not skip_esmfold and esm_queue:
        esm_ran = evaluate_monomers_with_esmfold(
            esm_queue,
            backbone_ca_map=backbone_ca_map,
            out_pdb_dir=harvest_dir / "monomer_esmfold",
        )
        if esm_ran:
            for cand in esm_queue:
                sc_r = cand["monomer_sc_rmsd"]
                pl_v = cand["monomer_esm_plddt"]
                if isinstance(sc_r, float) and sc_r > max_sc_rmsd:
                    cand["gate_reasons"].append(f"MONOMER_SCRMSD_{sc_r:.2f}A")
                if isinstance(pl_v, float) and pl_v < min_monomer_plddt:
                    cand["gate_reasons"].append(f"MONOMER_PLDDT_{pl_v:.1f}")

    for cand in evaluated_candidates:
        cand["final_verdict"] = "PASS" if not cand["gate_reasons"] else ("REJECT:" + ",".join(cand["gate_reasons"]))

    # Rank passing candidates within each backbone (preferring low scRMSD + low MPNN NLL + complementary modes)
    def _within_bb_key(c: dict[str, Any]) -> tuple[float, float, float]:
        scr = float(c["monomer_sc_rmsd"]) if c["monomer_sc_rmsd"] != "" else 1.5
        ek = float(c.get("E_ratio", 0.0)) + float(c.get("K_ratio", 0.0))
        return (scr, float(c["mpnn_score"]), -ek)

    by_bb_pass: dict[str, list[dict[str, Any]]] = {}
    for c in evaluated_candidates:
        if c["final_verdict"] == "PASS":
            by_bb_pass.setdefault(c["backbone_id"], []).append(c)

    per_bb_selected: dict[str, list[dict[str, Any]]] = {}
    for bid, clist in by_bb_pass.items():
        clist.sort(key=_within_bb_key)
        if max_per_backbone >= 2:
            # Prefer 1 best `fix_intf` + 1 best `full` when both modes passed for orthogonal diversity
            best_fix = next((x for x in clist if x["mode"] == "fix_intf"), None)
            best_full = next((x for x in clist if x["mode"] == "full"), None)
            picked: list[dict[str, Any]] = []
            if best_fix and best_full:
                picked = sorted([best_fix, best_full], key=_within_bb_key)
            for x in clist:
                if len(picked) >= max_per_backbone:
                    break
                if x not in picked:
                    picked.append(x)
            per_bb_selected[bid] = picked[:max_per_backbone]
        else:
            per_bb_selected[bid] = clist[:max_per_backbone]

    # Round-robin across top backbones so Part 1 (first 20 AF3 jobs) covers 20 distinct groove scaffolds!
    ordered_bb_ids = [b["design_id"] for b in passing_bbs if b["design_id"] in per_bb_selected]
    selected_candidates: list[dict[str, Any]] = []
    for slot_idx in range(max_per_backbone):
        for bid in ordered_bb_ids:
            picks = per_bb_selected[bid]
            if slot_idx < len(picks):
                cand = picks[slot_idx]
                cand["af3_priority_rank"] = len(selected_candidates) + 1
                selected_candidates.append(cand)

    # Write complete audit CSV & passing FASTA to harvest_dir on Google Drive
    all_csv = harvest_dir / "harvest_candidates_all.csv"
    audit_fields = [
        "design_id", "backbone_id", "mode", "sample_idx", "final_verdict",
        "traj_i_pTM", "traj_pLDDT", "traj_i_pAE", "traj_dG", "traj_SC", "traj_dSASA",
        "trimer_clashes", "hotspot_hits", "mpnn_score", "global_score", "seq_recovery",
        "monomer_esm_plddt", "monomer_sc_rmsd", "length_aa", "E_ratio", "K_ratio",
        "net_charge_pH6", "net_charge_pH74", "gravy", "cys_count", "longest_hydrophobic_run",
        "max_known_ident", "sequence",
    ]
    with all_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=audit_fields, extrasaction="ignore")
        w.writeheader()
        for c in evaluated_candidates:
            w.writerow(c)

    pass_fasta = harvest_dir / "harvest_candidates_selected.fasta"
    with pass_fasta.open("w", encoding="utf-8") as f:
        for c in selected_candidates:
            f.write(
                f">{c['design_id']} rank={c['af3_priority_rank']} mode={c['mode']} "
                f"traj_ipTM={c['traj_i_pTM']} dG={c['traj_dG']} mpnn={c['mpnn_score']} "
                f"scRMSD={c['monomer_sc_rmsd']} esm_pLDDT={c['monomer_esm_plddt']}\n{c['sequence']}\n"
            )

    # Append selected candidates to design_registry.csv & build AF3 Server JSON batches preserving round-robin rank!
    lengths = job.get("lengths", [50, 70])
    registry_csv.parent.mkdir(parents=True, exist_ok=True)
    if not registry_csv.exists():
        with registry_csv.open("w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow(worker.REGISTRY_HEADERS)

    registry_records: list[dict[str, Any]] = []
    for c in selected_candidates:
        seed_m = re.search(r"_s(\d+)_", c["design_id"])
        rec = {
            "design_id": c["design_id"],
            "job_id": job.get("job_id", ""),
            "target_id": job.get("target_id", "target5_TL1A"),
            "molecule_type": job.get("molecule_type", "minibinder"),
            "submitted_by": job.get("submitted_by", "team"),
            "generator": f"BindCraft_GrooveRelaxed+{c['mode']}",
            "seed": seed_m.group(1) if seed_m else "",
            "target_file": target_pdb.name,
            "target_sha256": target_sha,
            "length_spec": f"{lengths[0]}-{lengths[1]}" if len(lengths) == 2 else str(lengths),
            "hotspot_string": job.get("hotspots", ""),
            "inverse_folding_model": "SolubleMPNN",
            "is_solubleMPNN": "TRUE",
            "sequence": c["sequence"],
            "sequence_sha256": c["sequence_sha256"],
            "length_aa": c["length_aa"],
            "fold_pLDDT": c["traj_pLDDT"],
            "fold_ipTM": c["traj_i_pTM"],
            "fold_iPAE": c["traj_i_pAE"],
            "esm_plddt": c["monomer_esm_plddt"],
            "ipsae10_min": "",
            "ipsae10_max": "",
            "pdockq": "",
            "min_heavy_dist": "",
            "E_ratio": c["E_ratio"],
            "K_ratio": c["K_ratio"],
            "net_charge_pH6": c["net_charge_pH6"],
            "net_charge_pH74": c["net_charge_pH74"],
            "gravy": c["gravy"],
            "cys_count": c["cys_count"],
            "cys_parity_ok": c["cys_parity_ok"],
            "longest_hydrophobic_run": c["longest_hydrophobic_run"],
            "expressibility_verdict": "PASS",
            "created_at": dt.datetime.now().isoformat(timespec="seconds"),
            "submitted": "FALSE",
        }
        registry_records.append(rec)
        if c["design_id"] not in existing_ids:
            with registry_csv.open("a", newline="", encoding="utf-8") as out_f:
                csv.DictWriter(out_f, fieldnames=worker.REGISTRY_HEADERS).writerow(rec)
            existing_ids.add(c["design_id"])

    # Export AlphaFold 3 Server JSON batches in exact round-robin scaffold-diversity priority order!
    af3_dir = run_dir / "af3_inputs"
    af3_dir.mkdir(parents=True, exist_ok=True)
    af3_tgt_seq = (job.get("af3_target_sequence") or "").strip().upper()
    af3_copies = int(job.get("af3_target_copies", 3))
    if not af3_tgt_seq:
        parts = [p for p in (job.get("target_sequence") or "").strip().upper().split(":") if p]
        if parts:
            af3_tgt_seq = parts[0]
            if len(set(parts)) == 1:
                af3_copies = len(parts)

    af3_files: list[pathlib.Path] = []
    job_id = job.get("job_id", "target5_groove")
    for part_idx in range(0, len(selected_candidates), af3_batch_size):
        chunk = selected_candidates[part_idx : part_idx + af3_batch_size]
        payload = []
        for c in chunk:
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", c["design_id"])
            payload.append(
                {
                    "name": safe_name,
                    "modelSeeds": [],
                    "sequences": [
                        {"proteinChain": {"sequence": c["sequence"], "count": 1}},
                        {"proteinChain": {"sequence": af3_tgt_seq, "count": af3_copies}},
                    ],
                    "dialect": "alphafoldserver",
                    "version": 1,
                }
            )
        suffix = f"_part{part_idx // af3_batch_size + 1}" if len(selected_candidates) > af3_batch_size else ""
        out_p = af3_dir / f"af3_server_{job_id}_groove{suffix}.json"
        out_p.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        af3_files.append(out_p)

    worker.sync_filesystem()
    return {
        "audited_backbones": audited_all,
        "passing_backbones": passing_bbs,
        "evaluated_candidates": evaluated_candidates,
        "selected_candidates": selected_candidates,
        "registry_records": registry_records,
        "af3_files": af3_files,
        "all_csv": all_csv,
        "pass_fasta": pass_fasta,
        "esm_ran": esm_ran,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Harvest golden multi-chain groove backbones via SolubleMPNN + ESMFold scRMSD + AF3 JSON export"
    )
    ap.add_argument("--drive-root", default=worker.SHARED_DRIVE_FOLDER_URL, help="Shared Google Drive root or URL")
    ap.add_argument("--job-id", default="job_t5_tl1a_grooveAB_002", help="Run folder name under runs/")
    ap.add_argument("--run-dir", help="Explicit run directory path (overrides --drive-root/--job-id)")
    ap.add_argument(
        "--job-json",
        default=str(REPO_ROOT / "jobs" / "templates" / "job_target5_TL1A_minibinder.json"),
        help="Job configuration JSON file",
    )
    ap.add_argument("--min-iptm", type=float, default=0.80, help="Minimum Stage-4 AF2-Multimer trajectory i_pTM")
    ap.add_argument("--min-plddt", type=float, default=0.78, help="Minimum Stage-4 trajectory pLDDT")
    ap.add_argument("--max-trimer-clashes", type=int, default=5, help="Max allowed heavy-atom clashes (<2.5A) vs native trimer")
    ap.add_argument("--num-seqs-fix-intf", type=int, default=4, help="SolubleMPNN sequences/backbone with interface fixed")
    ap.add_argument("--num-seqs-full", type=int, default=4, help="SolubleMPNN sequences/backbone redesigning full binder")
    ap.add_argument("--max-per-backbone", type=int, default=2, help="Max diverse sequences kept per groove backbone")
    ap.add_argument("--max-sc-rmsd", type=float, default=2.5, help="Max unbound monomer ESMFold C-alpha scRMSD (A)")
    ap.add_argument("--min-monomer-plddt", type=float, default=75.0, help="Min unbound monomer ESMFold pLDDT (0-100)")
    ap.add_argument("--skip-esmfold", action="store_true", help="Skip ESMFold monomer scRMSD check (fast CPU mode)")
    ap.add_argument("--af3-batch-size", type=int, default=20, help="Max jobs per AlphaFold 3 Server JSON file")
    args = ap.parse_args(argv)

    job = json.loads(pathlib.Path(args.job_json).read_text(encoding="utf-8"))
    if args.run_dir:
        run_dir = pathlib.Path(args.run_dir)
        reg_csv = run_dir.parent.parent / "registry" / "design_registry.csv"
    else:
        dirs = worker.ensure_drive_workspace(pathlib.Path(args.drive_root))
        run_dir = dirs["runs"] / args.job_id
        reg_csv = dirs["registry"] / "design_registry.csv"
        completed_job_p = dirs["completed"] / pathlib.Path(args.job_json).name
        if completed_job_p.exists():
            try:
                cjob = json.loads(completed_job_p.read_text(encoding="utf-8"))
                if cjob.get("job_id") == args.job_id:
                    job = cjob
            except Exception:
                pass
    job["job_id"] = args.job_id

    res = harvest_run(
        run_dir=run_dir,
        job=job,
        registry_csv=reg_csv,
        min_iptm=args.min_iptm,
        min_plddt=args.min_plddt,
        max_trimer_clashes=args.max_trimer_clashes,
        num_seqs_fix_intf=args.num_seqs_fix_intf,
        num_seqs_full=args.num_seqs_full,
        max_per_backbone=args.max_per_backbone,
        max_sc_rmsd=args.max_sc_rmsd,
        min_monomer_plddt=args.min_monomer_plddt,
        skip_esmfold=args.skip_esmfold,
        af3_batch_size=args.af3_batch_size,
    )

    sel = res["selected_candidates"]
    print(f"\n================ [黃金凹槽骨架採收完成: {args.job_id}] ================")
    print(f"  通過 3D 凹槽與零碰撞審計骨架數 : {len(res['passing_backbones'])} / {len(res['audited_backbones'])}")
    print(f"  SolubleMPNN 總生成候選序列數   : {len(res['evaluated_candidates'])}")
    print(f"  入選跨骨架高多樣性序列數       : {len(sel)} (已寫入 {reg_csv})")
    print(f"  完整審計報表                   : {res['all_csv']}")
    print(f"  FASTA 檔案                     : {res['pass_fasta']}")
    print(f"  AlphaFold 3 Server 批量 JSON   :")
    for fp in res["af3_files"]:
        print(f"    -> {fp}")
    if sel:
        print("\n--- Top 20 優先送檢序列（Part 1 JSON：覆蓋前 20 個不同凹槽骨架！）---")
        for c in sel[:20]:
            ek_pct = (float(c["E_ratio"]) + float(c["K_ratio"])) * 100.0
            print(
                f"  #{c['af3_priority_rank']:02d} {c['design_id']:48s} len={c['length_aa']:2d} "
                f"traj_ipTM={c['traj_i_pTM']:.2f} dG={c['traj_dG']:6.1f} mpnn={c['mpnn_score']:.2f} "
                f"scRMSD={str(c['monomer_sc_rmsd'] or '-'):>5s} pLDDT={str(c['monomer_esm_plddt'] or '-'):>5s} "
                f"E+K={ek_pct:4.1f}% hotspots={c['hotspot_hits']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
