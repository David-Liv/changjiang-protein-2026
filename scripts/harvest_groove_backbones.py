#!/usr/bin/env python3
"""Harvest PyRosetta-relaxed multi-chain groove backbones via SolubleMPNN + 3D/physics gates + AF3 JSON export.

Why this module exists (SOP §D-1, §D-2 & §C-4):
- During BindCraft `default_4stage_multimer` runs on multi-chain targets (such as Target 5 TL1A's
  `A:B` inter-subunit DcR3/DR3 groove), Stage 1-4 trajectory generation uses `AF2-Multimer` weights
  and produces high-confidence relaxed complexes (`i_pTM = 0.75 - 0.92`, `0` native trimer clashes).
- However, BindCraft's internal MPNN complex re-prediction switches back to single-sequence, 0-MSA
  `AF2-Monomer-PTM` (`params_model_1_ptm`), which systematically blind-rejects shallow homodimer
  β-sheet grooves (see `docs/sop/epitope_and_expressibility_SOP_v1.pdf` §D-1, Adaptyv Round 3).
- Submitting raw AF2 trajectory sequences directly is forbidden by `AGENTS.md` / SOP §C-4; sequences
  MUST undergo `SolubleMPNN` inverse folding (`--omit_AAs CXM`), 3-mode design (`core_fix` pocket
  anchors + `full` redesign + guarded `fix_intf` on low-Unsat backbones), composition/expressibility
  gates (`Ala% <= 16%`, no `WW`, `E+K >= 26%`, balanced charge), ESMFold unbound-monomer `scRMSD`,
  3D pose clustering, N/C-terminus BLI clearance, and `3x WT TL1A + 1x Binder` AlphaFold 3 export.
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
CORE_ANCHOR_AAS = frozenset("FYWLIV")


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


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Compute rigid-body rotation R and centroids (cP, cQ) aligning P -> Q, plus RMSD."""
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
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - Q0) ** 2, axis=1))))
    return R, cP, cQ, rmsd


def kabsch_ca_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute optimal rigid-body C-alpha RMSD (in Angstroms) between two (L, 3) coordinate arrays."""
    if P.shape != Q.shape or len(P) < 3:
        return 999.0
    _, _, _, rmsd = kabsch_fit(P, Q)
    return rmsd


def build_receptor_groove_reference(
    target_pdb: pathlib.Path,
    native_trimer_pdb: pathlib.Path | None,
    receptor_complex_pdb: pathlib.Path | None,
) -> dict[str, Any] | None:
    """Superimpose a ligand-receptor crystal complex (e.g. 3K51 TL1A-DcR3) onto the active target groove.

    Returns reference Chain-A C-alpha coordinates (`ref_A_ca`, matching `target_pdb` Chain A residues)
    and aligned receptor heavy-atom coordinates (`receptor_heavy`) in the `native_trimer_pdb` frame.
    """
    if native_trimer_pdb is None or not native_trimer_pdb.exists():
        return None
    tgt_chains = parse_pdb_residues(target_pdb)
    ref_chains = parse_pdb_residues(native_trimer_pdb)
    if "A" not in tgt_chains or "A" not in ref_chains:
        return None

    tgt_A_resnums = [r["resseq"] for r in tgt_chains["A"] if "CA" in r["bb"]]
    ref_A_map = {r["resseq"]: r["bb"]["CA"] for r in ref_chains["A"] if "CA" in r["bb"]}
    if len(tgt_A_resnums) < 30 or any(rn not in ref_A_map for rn in tgt_A_resnums):
        return None
    ref_A_ca = np.array([ref_A_map[rn] for rn in tgt_A_resnums], dtype=float)

    receptor_heavy: np.ndarray | None = None
    if receptor_complex_pdb is not None and receptor_complex_pdb.exists():
        rc_chains = parse_pdb_residues(receptor_complex_pdb)
        if len(rc_chains) >= 2:
            # Match receptor_complex_pdb ligand chain to native_trimer_pdb Chain A by residue-number offset
            # (robust to missing internal loop density in lower-resolution receptor complexes such as 3K51)
            ref_by_rnum = {r["resseq"]: r["aa"] for r in ref_chains["A"] if "CA" in r["bb"]}
            lig_ch = None
            best_pairs: list[tuple[int, int]] = []
            for ch, rlist in rc_chains.items():
                rc_by_rnum = {r["resseq"]: (idx, r["aa"]) for idx, r in enumerate(rlist) if "CA" in r["bb"]}
                delta_votes: dict[int, int] = {}
                for r_ref, aa_ref in ref_by_rnum.items():
                    for r_rc, (_, aa_rc) in rc_by_rnum.items():
                        if aa_ref == aa_rc:
                            d_off = r_rc - r_ref
                            delta_votes[d_off] = delta_votes.get(d_off, 0) + 1
                for d_off, votes in sorted(delta_votes.items(), key=lambda kv: kv[1], reverse=True)[:5]:
                    if votes < 60:
                        break
                    overlap = [rn for rn in ref_by_rnum if (rn + d_off) in rc_by_rnum]
                    if len(overlap) < 60:
                        continue
                    matched = [
                        (rn, rc_by_rnum[rn + d_off][0])
                        for rn in overlap
                        if ref_by_rnum[rn] == rc_by_rnum[rn + d_off][1]
                    ]
                    if len(matched) / len(overlap) >= 0.85 and len(matched) > len(best_pairs):
                        best_pairs = matched
                        lig_ch = ch

            if lig_ch is not None and len(best_pairs) >= 60:
                rec_atoms_raw = np.vstack(
                    [np.asarray(r["heavy_atoms"]) for ch, rl in rc_chains.items() if ch != lig_ch for r in rl]
                )
                lig_rlist = rc_chains[lig_ch]
                active_tgt_chains = sorted(set(tgt_chains.keys()) & set(ref_chains.keys()))
                active_pts_by_ch = {
                    ch: np.vstack([np.asarray(r["heavy_atoms"]) for r in ref_chains[ch]])
                    for ch in active_tgt_chains
                }
                best_contacts = -1
                for sym_ch, sym_rlist in ref_chains.items():
                    sym_ca_map = {r["resseq"]: r["bb"]["CA"] for r in sym_rlist if "CA" in r["bb"]}
                    valid_pairs = [(rnum, lidx) for rnum, lidx in best_pairs if rnum in sym_ca_map]
                    if len(valid_pairs) < 60:
                        continue
                    P_lig = np.array([lig_rlist[lidx]["bb"]["CA"] for _, lidx in valid_pairs], dtype=float)
                    Q_sym = np.array([sym_ca_map[rnum] for rnum, _ in valid_pairs], dtype=float)
                    R_s, cP_s, cQ_s, rmsd_s = kabsch_fit(P_lig, Q_sym)
                    if rmsd_s > 3.0:
                        continue
                    cand_rec = (rec_atoms_raw - cP_s) @ R_s.T + cQ_s
                    # Count contacts across active groove chains; require touching all active chains if >=2
                    per_ch_hits = []
                    for ach in active_tgt_chains:
                        d = np.linalg.norm(cand_rec[:, None, :] - active_pts_by_ch[ach][None, :, :], axis=2)
                        per_ch_hits.append(int((d.min(axis=0) < 5.0).sum()))
                    if len(active_tgt_chains) >= 2 and any(h == 0 for h in per_ch_hits):
                        continue
                    tot_hits = sum(per_ch_hits)
                    if tot_hits > best_contacts:
                        best_contacts = tot_hits
                        receptor_heavy = cand_rec

    return {
        "ref_A_ca": ref_A_ca,
        "tgt_A_len": len(tgt_A_resnums),
        "receptor_heavy": receptor_heavy,
    }


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


def _select_core_anchor_positions(
    binder_res: list[dict[str, Any]],
    concat_target_res: list[dict[str, Any]],
    flat_labels: list[str],
    hotspots: set[str],
    contact_cutoff: float = 4.0,
) -> list[int]:
    """Select clean, deeply buried hydrophobic pocket-anchor positions (1-based) on the binder.

    Unlike `fix_intf` (which freezes 25+ residues including unsatisfied polar sidechains, Met, and
    WW patches), `core_fix` freezes ONLY 4-10 tight (`<=4.0A`) non-oxidizable hydrophobic/aromatic
    anchors (`F, Y, W, L, I, V`, max 2 non-adjacent `W`), freeing SolubleMPNN to optimize all rim
    hydrogen-bond networks, electrostatics, and solvent-exposed surface residues.
    """
    all_target_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in concat_target_res])
    hs_heavy_list = [
        np.asarray(r["heavy_atoms"])
        for idx, r in enumerate(concat_target_res)
        if flat_labels[idx] in hotspots
    ]
    hs_heavy = np.vstack(hs_heavy_list) if hs_heavy_list else None

    candidates: list[tuple[float, float, int, str]] = []
    for pos_idx, bres in enumerate(binder_res):
        aa = bres["aa"]
        if aa not in CORE_ANCHOR_AAS:
            continue
        bh = np.asarray(bres["heavy_atoms"])
        d_tgt = float(np.sqrt(np.min(np.sum((bh[:, None, :] - all_target_heavy[None, :, :]) ** 2, axis=2))))
        if d_tgt > contact_cutoff:
            continue
        if hs_heavy is not None:
            d_hs = float(np.sqrt(np.min(np.sum((bh[:, None, :] - hs_heavy[None, :, :]) ** 2, axis=2))))
        else:
            d_hs = d_tgt
        # Prioritize anchors closest to target groove & functional hotspots
        priority = d_tgt + 0.15 * min(d_hs, 8.0)
        candidates.append((priority, d_tgt, pos_idx + 1, aa))

    candidates.sort(key=lambda x: x[0])
    max_anchors = max(1, min(10, len(binder_res) // 6))
    selected: list[int] = []
    selected_w: set[int] = set()

    for _, _, pos1, aa in candidates:
        if len(selected) >= max_anchors:
            break
        if aa == "W":
            if len(selected_w) >= 2 or (pos1 - 1) in selected_w or (pos1 + 1) in selected_w:
                continue
            selected_w.add(pos1)
        selected.append(pos1)

    if len(selected) >= len(binder_res):
        selected = selected[:-1]
    return sorted(selected)


def audit_and_prepare_backbone(
    relaxed_pdb: pathlib.Path,
    traj_row: dict[str, str],
    target_pdb: pathlib.Path,
    target_meta: dict[str, Any],
    native_trimer_pdb: pathlib.Path | None,
    hotspots: set[str],
    receptor_ref: dict[str, Any] | None = None,
    min_iptm: float = 0.75,
    min_plddt: float = 0.78,
    min_sc: float = 0.54,
    min_second_chain_contacts: int = 3,
    max_trimer_clashes: int = 5,
    max_fix_unsat: float = 4.0,
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
    packstat = float(traj_row.get("PackStat") or 0.0)
    dsasa = float(traj_row.get("dSASA") or 0.0)
    hbonds = int(float(traj_row.get("n_InterfaceHbonds") or 0.0))
    unsat = float(
        traj_row.get("n_InterfaceUnsatHbonds")
        if traj_row.get("n_InterfaceUnsatHbonds") not in (None, "")
        else (traj_row.get("UnsatHbonds") or 0.0)
    )
    surf_hyd = float(traj_row.get("Surface_Hydrophobicity") or 0.0)
    helix_pct = float(traj_row.get("Binder_Helix%") or 0.0)
    strand_pct = float(traj_row.get("Binder_BetaSheet%") or 0.0)
    loop_pct = float(traj_row.get("Binder_Loop%") or 0.0)
    traj_seq = (traj_row.get("Sequence") or "").strip().upper()

    result: dict[str, Any] = {
        "design_id": design_id,
        "pdb_path": relaxed_pdb,
        "i_pTM": iptm,
        "pLDDT": plddt,
        "i_pAE": ipae,
        "dG": dg,
        "SC": sc,
        "PackStat": packstat,
        "dSASA": dsasa,
        "interface_hbonds": hbonds,
        "unsat_hbonds": unsat,
        "surf_hydrophobicity": surf_hyd,
        "binder_helix_pct": helix_pct,
        "binder_strand_pct": strand_pct,
        "binder_loop_pct": loop_pct,
        "is_mixed_ab": strand_pct >= 15.0,
        "traj_seq": traj_seq,
        "passed": False,
        "reject_reasons": [],
    }

    if iptm < min_iptm:
        result["reject_reasons"].append(f"LOW_IPTM_{iptm:.2f}")
    if plddt < min_plddt:
        result["reject_reasons"].append(f"LOW_PLDDT_{plddt:.2f}")
    if sc > 0.0 and sc < min_sc:
        result["reject_reasons"].append(f"LOW_SC_{sc:.2f}")
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

    # Compute heavy-atom contacts (<5.0 A), N/C-terminus clearance (SOP §D-1.5), & distance to masked LYS
    binder_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in binder_res])
    all_target_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in concat_target_res])
    flat_labels: list[str] = target_meta["flat_labels"]
    masked_labels: set[str] = target_meta["masked_labels"]

    n_slice = binder_res[: min(2, len(binder_res))]
    c_slice = binder_res[max(0, len(binder_res) - 2) :]
    nterm_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in n_slice])
    cterm_heavy = np.vstack([np.asarray(r["heavy_atoms"]) for r in c_slice])
    d_nterm = float(np.sqrt(np.min(np.sum((nterm_heavy[:, None, :] - all_target_heavy[None, :, :]) ** 2, axis=2))))
    d_cterm = float(np.sqrt(np.min(np.sum((cterm_heavy[:, None, :] - all_target_heavy[None, :, :]) ** 2, axis=2))))
    result["d_nterm"] = round(d_nterm, 2)
    result["d_cterm"] = round(d_cterm, 2)
    result["min_term_dist"] = round(min(d_nterm, d_cterm), 2)

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
    result["contact_labels"] = contact_labels
    result["contacts_by_chain"] = contacts_by_orig_chain
    result["min_chain_contacts"] = min(contacts_by_orig_chain.values()) if contacts_by_orig_chain else 0
    result["hotspot_hits"] = hotspot_hits
    result["has_b_hotspot"] = any(h.startswith("B") for h in hotspot_hits)
    result["min_mask_dist"] = round(min_mask_dist, 2) if masked_labels else 999.0

    # Align binder into common reference frame for 3D pose clustering & DcR3 receptor competition check
    binder_ca_local = np.array([r["bb"]["CA"] for r in binder_res], dtype=float)
    receptor_clashes = -1
    if receptor_ref is not None and len(concat_target_res) >= receptor_ref["tgt_A_len"]:
        P_tA = np.array([r["bb"]["CA"] for r in concat_target_res[: receptor_ref["tgt_A_len"]]], dtype=float)
        R_ref, cP_ref, cQ_ref, _ = kabsch_fit(P_tA, receptor_ref["ref_A_ca"])
        binder_ca_ref = (binder_ca_local - cP_ref) @ R_ref.T + cQ_ref
        if receptor_ref["receptor_heavy"] is not None:
            b_heavy_ref = (binder_heavy - cP_ref) @ R_ref.T + cQ_ref
            d_rec = np.linalg.norm(b_heavy_ref[:, None, :] - receptor_ref["receptor_heavy"][None, :, :], axis=2)
            receptor_clashes = int((d_rec.min(axis=1) < 3.0).sum())
    else:
        binder_ca_ref = binder_ca_local
    result["com_ref"] = binder_ca_ref.mean(axis=0)
    result["receptor_clashes"] = receptor_clashes

    if masked_labels and min_mask_dist < 8.0:
        result["reject_reasons"].append(f"NEAR_MASKED_RESIDUE_{min_mask_dist:.1f}A")
    if require_groove_bridge and len(target_meta["chain_specs"]) >= 2:
        chains_touched = sum(1 for n in contacts_by_orig_chain.values() if n > 0)
        if chains_touched < 2:
            result["reject_reasons"].append("SINGLE_CHAIN_ONLY_NOT_GROOVE")
        elif result["min_chain_contacts"] < min_second_chain_contacts:
            result["reject_reasons"].append(f"WEAK_GROOVE_BRIDGE_MIN_{result['min_chain_contacts']}")
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
            if 1 <= pos <= len(binder_res) and binder_pdb_seq[pos - 1] not in ("C", "M"):
                fixed_binder_positions.append(pos)
    fixed_binder_positions = sorted(set(fixed_binder_positions))
    if len(fixed_binder_positions) >= len(binder_res):
        fixed_binder_positions = fixed_binder_positions[:-1]

    core_fixed_positions = _select_core_anchor_positions(
        binder_res=binder_res,
        concat_target_res=concat_target_res,
        flat_labels=flat_labels,
        hotspots=hotspots,
    )

    # Guard `fix_intf`: never freeze interfaces with high UnsatHbonds, high Surface_Hydrophobicity, or fixed WW/excess W
    fixed_w_pos = {p for p in fixed_binder_positions if binder_pdb_seq[p - 1] == "W"}
    has_fixed_ww = any((p + 1) in fixed_w_pos for p in fixed_w_pos)
    allow_fix_intf = (
        unsat <= max_fix_unsat
        and (surf_hyd <= 0.42 or surf_hyd == 0.0)
        and (not has_fixed_ww)
        and (len(fixed_w_pos) <= 3)
    )

    result["mpnn_entry"] = mpnn_entry
    result["fixed_chains"] = fixed_chains
    result["binder_mpnn_chain"] = binder_mpnn_ch
    result["fixed_binder_positions"] = fixed_binder_positions
    result["core_fixed_binder_positions"] = core_fixed_positions
    result["allow_fix_intf"] = allow_fix_intf
    result["binder_ca"] = binder_ca_local
    result["passed"] = len(result["reject_reasons"]) == 0
    return result


def compute_backbone_quality_score(bb: dict[str, Any]) -> float:
    """Multi-objective physical quality score for ordering backbones (SOP §5: avoid raw i_pTM over-ranking)."""
    iptm_term = min(bb["i_pTM"], 0.90) * 1.0
    sc_term = bb["SC"] * 1.2
    dg_term = min(max(-bb["dG"], 0.0), 54.0) / 120.0
    hb_term = min(bb["interface_hbonds"], 12) * 0.025
    bridge_term = min(bb.get("min_chain_contacts", 0), 12) * 0.03
    unsat_pen = max(0.0, bb["unsat_hbonds"] - 3.0) * 0.075
    surf_pen = max(0.0, bb["surf_hydrophobicity"] - 0.36) * 0.5

    # SOP §D-1.5: BLI tag immobilization clearance at N/C termini
    d_min_t = bb.get("min_term_dist", 10.0)
    term_adj = 0.08 if d_min_t >= 4.8 else (-0.12 if d_min_t < 3.5 else 0.0)
    if bb.get("d_nterm", 10.0) < 4.0 and bb.get("d_cterm", 10.0) < 4.0:
        term_adj -= 0.12

    # SOP §D-2 diversity & functional anchor bonuses
    div_bonus = 0.0
    if bb.get("is_mixed_ab"):
        div_bonus += 0.22
    if bb.get("has_b_hotspot"):
        div_bonus += 0.18
    if bb.get("receptor_clashes", -1) >= 200:
        div_bonus += 0.05

    return round(iptm_term + sc_term + dg_term + hb_term + bridge_term + term_adj + div_bonus - unsat_pen - surf_pen, 4)


def cluster_and_rank_backbones(passing_bbs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Score backbones by physical quality and group near-duplicate 3D poses into spatial clusters."""
    for bb in passing_bbs:
        bb["bb_quality_score"] = compute_backbone_quality_score(bb)

    ranked = sorted(
        passing_bbs,
        key=lambda x: (x["bb_quality_score"], x["i_pTM"], -x["dG"]),
        reverse=True,
    )

    clusters: list[list[dict[str, Any]]] = []
    for bb in ranked:
        placed = False
        s1 = bb.get("contact_labels", set())
        c1 = bb.get("com_ref")
        for cl_idx, cl in enumerate(clusters, start=1):
            rep = cl[0]
            s2 = rep.get("contact_labels", set())
            c2 = rep.get("com_ref")
            jacc = len(s1 & s2) / max(1, len(s1 | s2)) if (s1 or s2) else 0.0
            com_d = float(np.linalg.norm(c1 - c2)) if (c1 is not None and c2 is not None) else 999.0
            # Keep mixed alpha/beta and B-wall hotspot scaffolds in distinct clusters from plain helices
            same_special = (bb.get("is_mixed_ab") == rep.get("is_mixed_ab")) and (
                bb.get("has_b_hotspot") == rep.get("has_b_hotspot")
            )
            if same_special and jacc >= 0.55 and com_d <= 6.0:
                bb["pose_cluster"] = f"C{cl_idx:02d}"
                bb["cluster_rank"] = len(cl) + 1
                cl.append(bb)
                placed = True
                break
        if not placed:
            new_id = len(clusters) + 1
            bb["pose_cluster"] = f"C{new_id:02d}"
            bb["cluster_rank"] = 1
            clusters.append([bb])

    return ranked


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
    positions_key: str = "fixed_binder_positions",
    require_clean_fix: bool = False,
    omit_aas: str = "CXM",
    seed: int = 37,
    python_exe: str | None = None,
) -> list[dict[str, Any]]:
    """Run official SolubleMPNN (`v_48_020`) across audited groove backbones."""
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
            if require_clean_fix and not bb.get("allow_fix_intf", True):
                continue
            pos_to_fix = bb.get(positions_key, []) if fix_interface else []
            if fix_interface and not pos_to_fix:
                continue
            active_backbones.append(bb)
            name = bb["design_id"]
            pf.write(json.dumps(bb["mpnn_entry"]) + "\n")
            b_ch = bb["binder_mpnn_chain"]
            f_chs = bb["fixed_chains"]
            chain_dict[name] = [[b_ch], f_chs]
            if fix_interface:
                pos_map = {ch: [] for ch in f_chs}
                pos_map[b_ch] = pos_to_fix
                fixed_dict[name] = pos_map

    if not active_backbones:
        return []

    chain_jsonl.write_text(json.dumps(chain_dict) + "\n", encoding="utf-8")
    if fix_interface:
        fixed_jsonl.write_text(json.dumps(fixed_dict) + "\n", encoding="utf-8")

    py_bin = python_exe or sys.executable
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
        omit_aas,
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
    min_iptm: float = 0.75,
    min_plddt: float = 0.78,
    min_sc: float = 0.54,
    min_second_chain_contacts: int = 3,
    max_trimer_clashes: int = 5,
    num_seqs_core_fix: int = 6,
    num_seqs_fix_intf: int = 4,
    num_seqs_full: int = 6,
    temp_core_fix: float = 0.12,
    temp_fix_intf: float = 0.10,
    temp_full: float = 0.15,
    omit_aas: str = "CXM",
    max_ala_ratio: float = 0.16,
    min_ek_ratio: float = 0.26,
    max_trp_count: int = 3,
    max_met_count: int = 3,
    max_abs_charge_ph74: float = 10.5,
    max_per_backbone: int = 2,
    max_sc_rmsd: float = 2.2,
    min_monomer_plddt: float = 78.0,
    skip_esmfold: bool = False,
    af3_batch_size: int = 20,
    mpnn_runner_override=None,
) -> dict[str, Any]:
    """End-to-end harvest of relaxed groove backbones -> SolubleMPNN -> 3D/seq gates -> AF3 Server JSON."""
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

    receptor_complex_pdb: pathlib.Path | None = None
    rec_rel = job.get("receptor_complex_pdb_relpath")
    if rec_rel:
        cand_rec = (REPO_ROOT / rec_rel).resolve()
        if cand_rec.exists():
            receptor_complex_pdb = cand_rec
    elif job.get("target_id") == "target5_TL1A":
        default_3k51 = REPO_ROOT / "targets" / "target5_TL1A" / "evidence" / "3K51.pdb"
        if default_3k51.exists():
            receptor_complex_pdb = default_3k51

    hotspots = {
        h.strip()
        for h in (job.get("hotspots") or "").split(",")
        if h.strip()
    }
    hard_len = job.get("lengths_hard_limit", [40, 100] if job.get("molecule_type") == "minibinder" else [10, 250])
    min_len_hard, max_len_hard = int(hard_len[0]), int(hard_len[1])

    target_meta = build_target_reference_meta(target_pdb, native_trimer_pdb)
    receptor_ref = build_receptor_groove_reference(target_pdb, native_trimer_pdb, receptor_complex_pdb)
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
            receptor_ref=receptor_ref,
            min_iptm=min_iptm,
            min_plddt=min_plddt,
            min_sc=min_sc,
            min_second_chain_contacts=min_second_chain_contacts,
            max_trimer_clashes=max_trimer_clashes,
        )
        audited_all.append(aud)

    passing_bbs = cluster_and_rank_backbones([a for a in audited_all if a["passed"]])
    n_clusters = len({b["pose_cluster"] for b in passing_bbs})
    n_mixed = sum(1 for b in passing_bbs if b.get("is_mixed_ab"))
    n_b_hs = sum(1 for b in passing_bbs if b.get("has_b_hotspot"))
    print(
        f"[Groove Harvest] Audited {len(audited_all)} relaxed trajectories -> "
        f"{len(passing_bbs)} passed 3D groove/SC/0-clash gates across {n_clusters} pose clusters "
        f"(mixed alpha/beta={n_mixed}, B-wall hotspot={n_b_hs}).",
        flush=True,
    )
    if not passing_bbs:
        return {"selected_records": [], "selected_candidates": [], "af3_files": [], "audited_backbones": audited_all}

    harvest_dir = run_dir / "harvest_soluble_mpnn"
    harvest_dir.mkdir(parents=True, exist_ok=True)

    if mpnn_runner_override is not None:
        raw_mpnn = mpnn_runner_override(passing_bbs, harvest_dir)
    else:
        mpnn_repo = ensure_protein_mpnn_repo()
        raw_mpnn = []
        if num_seqs_core_fix > 0:
            print(
                f"[SolubleMPNN] Generating {num_seqs_core_fix} pocket-anchor-fixed (`core_fix`) sequences/backbone "
                f"(T={temp_core_fix}, omit_AAs={omit_aas})...",
                flush=True,
            )
            raw_mpnn.extend(
                run_soluble_mpnn_mode(
                    passing_bbs,
                    harvest_dir,
                    mpnn_repo,
                    mode_tag="core_fix",
                    num_seqs=num_seqs_core_fix,
                    sampling_temp=temp_core_fix,
                    fix_interface=True,
                    positions_key="core_fixed_binder_positions",
                    require_clean_fix=False,
                    omit_aas=omit_aas,
                    seed=101,
                )
            )
        if num_seqs_full > 0:
            print(
                f"[SolubleMPNN] Generating {num_seqs_full} full-redesign (`full`) sequences/backbone "
                f"(T={temp_full}, omit_AAs={omit_aas})...",
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
                    omit_aas=omit_aas,
                    seed=73,
                )
            )
        if num_seqs_fix_intf > 0:
            n_clean_fix = sum(1 for b in passing_bbs if b.get("allow_fix_intf"))
            print(
                f"[SolubleMPNN] Generating {num_seqs_fix_intf} full-interface-fixed (`fix_intf`) sequences on "
                f"{n_clean_fix}/{len(passing_bbs)} low-Unsat clean backbones (T={temp_fix_intf}, omit_AAs={omit_aas})...",
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
                    positions_key="fixed_binder_positions",
                    require_clean_fix=True,
                    omit_aas=omit_aas,
                    seed=37,
                )
            )

    bb_by_id = {b["design_id"]: b for b in passing_bbs}
    backbone_ca_map = {b["design_id"]: b["binder_ca"] for b in passing_bbs}

    # Load known binders for SOP §C-1 novelty check
    known_fasta = REPO_ROOT / "data" / "raw" / "known_binders.fasta"
    known_seqs = check_novelty.read_fasta(str(known_fasta)) if known_fasta.exists() else {}

    # Load registry rows from OTHER jobs so re-harvesting the same job_id is 100% idempotent
    cur_job_id = str(job.get("job_id", ""))
    retained_registry_rows: list[dict[str, str]] = []
    other_job_seq_shas: set[str] = set()
    other_job_ids: set[str] = set()
    if registry_csv.exists():
        with registry_csv.open("r", newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                is_same_harvest = (
                    cur_job_id
                    and row.get("job_id") == cur_job_id
                    and str(row.get("generator", "")).startswith("BindCraft_GrooveRelaxed")
                )
                if is_same_harvest:
                    continue
                retained_registry_rows.append(row)
                if row.get("design_id"):
                    other_job_ids.add(row["design_id"])
                if row.get("sequence_sha256") and row.get("expressibility_verdict") == "PASS":
                    other_job_seq_shas.add(row["sequence_sha256"])

    mode_short_map = {"core_fix": "core", "fix_intf": "fix", "full": "full"}
    seen_batch_seqs: set[str] = set()
    evaluated_candidates: list[dict[str, Any]] = []
    esm_queue: list[dict[str, Any]] = []

    for rec in raw_mpnn:
        seq = rec["sequence"]
        bb = bb_by_id[rec["backbone_id"]]
        short_mode = mode_short_map.get(rec["mode"], rec["mode"])
        design_id = f"{rec['backbone_id']}_smpnn_{short_mode}{rec['sample_idx']}"
        seq_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()

        exp = worker.compute_expressibility_metrics(seq)
        ala_ratio = round(seq.count("A") / max(1, len(seq)), 4)
        ek_sum = round(float(exp.get("E_ratio", 0.0)) + float(exp.get("K_ratio", 0.0)), 4)
        trp_count = seq.count("W")
        met_count = seq.count("M")
        q74 = float(exp.get("net_charge_pH74", 0.0))

        reasons: list[str] = []
        if exp["expressibility_verdict"] != "PASS":
            reasons.append(exp["expressibility_verdict"].replace("REJECT:", ""))
        if len(seq) < min_len_hard or len(seq) > max_len_hard:
            reasons.append(f"HARD_LEN_{len(seq)}_NOT_IN_{min_len_hard}_{max_len_hard}")
        if rec["mode"] == "fix_intf" and not bb.get("allow_fix_intf", True):
            reasons.append(f"UNCLEAN_FIX_INTF_UNSAT_{bb['unsat_hbonds']:.0f}")
        if ala_ratio > max_ala_ratio:
            reasons.append(f"HIGH_ALA_{ala_ratio * 100:.1f}PCT")
        if "AAAA" in seq:
            reasons.append("POLY_ALA_RUN")
        if "WW" in seq:
            reasons.append("DI_TRP_PATCH")
        if trp_count > max_trp_count:
            reasons.append(f"HIGH_TRP_{trp_count}")
        if met_count > max_met_count:
            reasons.append(f"HIGH_MET_{met_count}")
        if ek_sum < min_ek_ratio:
            reasons.append(f"LOW_EK_{ek_sum * 100:.1f}PCT")
        if abs(q74) > max_abs_charge_ph74:
            reasons.append(f"EXTREME_CHARGE_{q74:+.1f}")
        if seq in seen_batch_seqs:
            reasons.append("DUPLICATE_IN_BATCH")
        elif seq_sha in other_job_seq_shas:
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
            "pose_cluster": bb["pose_cluster"],
            "bb_quality_score": bb["bb_quality_score"],
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
            "traj_unsat_hbonds": bb["unsat_hbonds"],
            "traj_hbonds": bb["interface_hbonds"],
            "binder_strand_pct": bb["binder_strand_pct"],
            "trimer_clashes": bb["trimer_clashes"],
            "receptor_clashes": bb["receptor_clashes"],
            "d_nterm": bb["d_nterm"],
            "d_cterm": bb["d_cterm"],
            "hotspot_hits": ";".join(bb["hotspot_hits"]),
            "contacts_by_chain": json.dumps(bb["contacts_by_chain"]),
            "ala_ratio": ala_ratio,
            "trp_count": trp_count,
            "met_count": met_count,
            "max_known_ident": max_known_ident,
            "monomer_esm_plddt": rec.get("monomer_esm_plddt", ""),
            "monomer_sc_rmsd": rec.get("monomer_sc_rmsd", ""),
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
    for cand in esm_queue:
        sc_r = cand["monomer_sc_rmsd"]
        pl_v = cand["monomer_esm_plddt"]
        bb = bb_by_id[cand["backbone_id"]]
        req_plddt = max(74.0, min_monomer_plddt - 4.0) if bb.get("is_mixed_ab") else min_monomer_plddt
        if sc_r != "" and float(sc_r) > max_sc_rmsd:
            cand["gate_reasons"].append(f"MONOMER_SCRMSD_{float(sc_r):.2f}A")
        if pl_v != "" and float(pl_v) < req_plddt:
            cand["gate_reasons"].append(f"MONOMER_PLDDT_{float(pl_v):.1f}")

    for cand in evaluated_candidates:
        cand["final_verdict"] = "PASS" if not cand["gate_reasons"] else ("REJECT:" + ",".join(cand["gate_reasons"]))

    # Plateaued composite rank within each backbone (avoids sacrificing solubility for 0.1A sub-atomic RMSD noise)
    def _within_bb_key(c: dict[str, Any]) -> tuple[float, float]:
        scr = float(c["monomer_sc_rmsd"]) if c["monomer_sc_rmsd"] != "" else 1.2
        pl = float(c["monomer_esm_plddt"]) if c["monomer_esm_plddt"] != "" else 85.0
        ek = float(c.get("E_ratio", 0.0)) + float(c.get("K_ratio", 0.0))
        ala = float(c.get("ala_ratio", 0.0))
        q_abs = abs(float(c.get("net_charge_pH74", 0.0)) + 3.0)
        rmsd_pen = max(0.0, scr - 1.05) * 1.2 + scr * 0.15
        pl_bonus = min(pl, 92.0) * 0.015
        ek_bonus = min(ek, 0.48) * 0.8
        ala_pen = max(0.0, ala - 0.10) * 1.5
        q_pen = max(0.0, q_abs - 4.5) * 0.04
        mode_pref = -0.04 if c["mode"] == "core_fix" else 0.0
        comp = rmsd_pen + float(c["mpnn_score"]) * 0.65 + ala_pen + q_pen + mode_pref - pl_bonus - ek_bonus
        return (round(comp, 5), scr)

    by_bb_pass: dict[str, list[dict[str, Any]]] = {}
    for c in evaluated_candidates:
        if c["final_verdict"] == "PASS":
            by_bb_pass.setdefault(c["backbone_id"], []).append(c)

    per_bb_selected: dict[str, list[dict[str, Any]]] = {}
    for bid, clist in by_bb_pass.items():
        clist.sort(key=_within_bb_key)
        if max_per_backbone >= 2:
            # Prefer 1 best anchor-preserving (`core_fix` or clean `fix_intf`) + 1 best `full` redesign
            best_anchored = next((x for x in clist if x["mode"] in ("core_fix", "fix_intf")), None)
            best_full = next((x for x in clist if x["mode"] == "full"), None)
            picked: list[dict[str, Any]] = []
            if best_anchored and best_full:
                picked = sorted([best_anchored, best_full], key=_within_bb_key)
            for x in clist:
                if len(picked) >= max_per_backbone:
                    break
                if x not in picked:
                    picked.append(x)
            per_bb_selected[bid] = picked[:max_per_backbone]
        else:
            per_bb_selected[bid] = clist[:max_per_backbone]

    # Cluster-first round-robin ordering:
    # Wave 1: Top backbone of EVERY distinct 3D pose cluster (guarantees Part 1 covers all poses, mixed alpha/beta & B317!)
    # Wave 2: Remaining distinct backbones ordered by physical quality score
    cluster_heads = [b["design_id"] for b in passing_bbs if b["design_id"] in per_bb_selected and b.get("cluster_rank", 1) == 1]
    cluster_tails = [b["design_id"] for b in passing_bbs if b["design_id"] in per_bb_selected and b.get("cluster_rank", 1) > 1]
    ordered_bb_ids = cluster_heads + cluster_tails

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
        "design_id", "backbone_id", "pose_cluster", "bb_quality_score", "mode", "sample_idx", "final_verdict",
        "traj_i_pTM", "traj_pLDDT", "traj_i_pAE", "traj_dG", "traj_SC", "traj_dSASA",
        "traj_hbonds", "traj_unsat_hbonds", "binder_strand_pct", "trimer_clashes", "receptor_clashes",
        "d_nterm", "d_cterm", "hotspot_hits", "mpnn_score", "global_score", "seq_recovery",
        "monomer_esm_plddt", "monomer_sc_rmsd", "length_aa", "E_ratio", "K_ratio", "ala_ratio",
        "trp_count", "met_count", "net_charge_pH6", "net_charge_pH74", "gravy", "cys_count",
        "longest_hydrophobic_run", "max_known_ident", "sequence",
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
                f">{c['design_id']} rank={c['af3_priority_rank']} cluster={c['pose_cluster']} mode={c['mode']} "
                f"bbQ={c['bb_quality_score']} traj_ipTM={c['traj_i_pTM']} SC={c['traj_SC']} "
                f"Unsat={c['traj_unsat_hbonds']} scRMSD={c['monomer_sc_rmsd']} esm_pLDDT={c['monomer_esm_plddt']}\n"
                f"{c['sequence']}\n"
            )

    # Rewrite design_registry.csv idempotently (retaining all other jobs + writing updated harvest for cur_job_id)
    lengths = job.get("lengths", [50, 70])
    registry_csv.parent.mkdir(parents=True, exist_ok=True)
    registry_records: list[dict[str, Any]] = []
    for c in selected_candidates:
        if c["design_id"] in other_job_ids:
            continue
        seed_m = re.search(r"_s(\d+)_", c["design_id"])
        rec = {
            "design_id": c["design_id"],
            "job_id": cur_job_id,
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

    with registry_csv.open("w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=worker.REGISTRY_HEADERS, extrasaction="ignore")
        writer.writeheader()
        for old_row in retained_registry_rows:
            writer.writerow(old_row)
        for new_row in registry_records:
            writer.writerow(new_row)

    # Export AlphaFold 3 Server JSON batches in cluster-diverse physical priority order (cleaning stale parts first)
    af3_dir = run_dir / "af3_inputs"
    af3_dir.mkdir(parents=True, exist_ok=True)
    job_id = cur_job_id or "target5_groove"
    for stale_json in af3_dir.glob(f"af3_server_{job_id}_groove*.json"):
        stale_json.unlink()

    af3_tgt_seq = (job.get("af3_target_sequence") or "").strip().upper()
    af3_copies = int(job.get("af3_target_copies", 3))
    if not af3_tgt_seq:
        parts = [p for p in (job.get("target_sequence") or "").strip().upper().split(":") if p]
        if parts:
            af3_tgt_seq = parts[0]
            if len(set(parts)) == 1:
                af3_copies = len(parts)

    af3_files: list[pathlib.Path] = []
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
        description="Harvest golden multi-chain groove backbones via 3-mode SolubleMPNN + ESMFold scRMSD + AF3 JSON export"
    )
    ap.add_argument("--drive-root", default=worker.SHARED_DRIVE_FOLDER_URL, help="Shared Google Drive root or URL")
    ap.add_argument("--job-id", default="job_t5_tl1a_grooveAB_002", help="Run folder name under runs/")
    ap.add_argument("--run-dir", help="Explicit run directory path (overrides --drive-root/--job-id)")
    ap.add_argument(
        "--job-json",
        default=str(REPO_ROOT / "jobs" / "templates" / "job_target5_TL1A_minibinder.json"),
        help="Job configuration JSON file",
    )
    ap.add_argument("--min-iptm", type=float, default=0.75, help="Minimum Stage-4 AF2-Multimer trajectory i_pTM")
    ap.add_argument("--min-plddt", type=float, default=0.78, help="Minimum Stage-4 trajectory pLDDT")
    ap.add_argument("--min-sc", type=float, default=0.54, help="Minimum PyRosetta ShapeComplementarity")
    ap.add_argument("--min-second-chain-contacts", type=int, default=3, help="Min target contacts on secondary groove chain")
    ap.add_argument("--max-trimer-clashes", type=int, default=5, help="Max allowed heavy-atom clashes (<2.5A) vs native trimer")
    ap.add_argument("--num-seqs-core-fix", type=int, default=6, help="SolubleMPNN sequences/backbone with pocket anchors fixed")
    ap.add_argument("--num-seqs-fix-intf", type=int, default=4, help="SolubleMPNN sequences/backbone with full interface fixed (clean Unsat<=4 only)")
    ap.add_argument("--num-seqs-full", type=int, default=6, help="SolubleMPNN sequences/backbone redesigning full binder")
    ap.add_argument("--omit-aas", default="CXM", help="Amino acids omitted by SolubleMPNN (default: CXM)")
    ap.add_argument("--max-ala-ratio", type=float, default=0.16, help="Maximum allowed Alanine fraction (SOP §C-2)")
    ap.add_argument("--min-ek-ratio", type=float, default=0.26, help="Minimum allowed E+K fraction (SOP §C-2)")
    ap.add_argument("--max-per-backbone", type=int, default=2, help="Max diverse sequences kept per groove backbone")
    ap.add_argument("--max-sc-rmsd", type=float, default=2.2, help="Max unbound monomer ESMFold C-alpha scRMSD (A)")
    ap.add_argument("--min-monomer-plddt", type=float, default=78.0, help="Min unbound monomer ESMFold pLDDT (0-100)")
    ap.add_argument("--skip-esmfold", action="store_true", help="Skip ESMFold monomer scRMSD check (fast CPU mode)")
    ap.add_argument("--af3-batch-size", type=int, default=20, help="Max jobs per AlphaFold 3 Server JSON file")
    args = ap.parse_args(argv)

    # Auto-pick matching template if job_id is 003 and user kept default --job-json
    job_json_path = pathlib.Path(args.job_json)
    alt_003 = REPO_ROOT / "jobs" / "templates" / "job_target5_TL1A_grooveAB_003.json"
    if args.job_id.endswith("_003") and job_json_path.name == "job_target5_TL1A_minibinder.json" and alt_003.exists():
        job_json_path = alt_003

    job = json.loads(job_json_path.read_text(encoding="utf-8"))
    if args.run_dir:
        run_dir = pathlib.Path(args.run_dir)
        reg_csv = run_dir.parent.parent / "registry" / "design_registry.csv"
    else:
        dirs = worker.ensure_drive_workspace(pathlib.Path(args.drive_root))
        run_dir = dirs["runs"] / args.job_id
        reg_csv = dirs["registry"] / "design_registry.csv"
        for completed_job_p in sorted(dirs["completed"].glob("*.json")):
            try:
                cjob = json.loads(completed_job_p.read_text(encoding="utf-8"))
                if cjob.get("job_id") == args.job_id:
                    job = cjob
                    break
            except Exception:
                pass
    job["job_id"] = args.job_id

    res = harvest_run(
        run_dir=run_dir,
        job=job,
        registry_csv=reg_csv,
        min_iptm=args.min_iptm,
        min_plddt=args.min_plddt,
        min_sc=args.min_sc,
        min_second_chain_contacts=args.min_second_chain_contacts,
        max_trimer_clashes=args.max_trimer_clashes,
        num_seqs_core_fix=args.num_seqs_core_fix,
        num_seqs_fix_intf=args.num_seqs_fix_intf,
        num_seqs_full=args.num_seqs_full,
        omit_aas=args.omit_aas,
        max_ala_ratio=args.max_ala_ratio,
        min_ek_ratio=args.min_ek_ratio,
        max_per_backbone=args.max_per_backbone,
        max_sc_rmsd=args.max_sc_rmsd,
        min_monomer_plddt=args.min_monomer_plddt,
        skip_esmfold=args.skip_esmfold,
        af3_batch_size=args.af3_batch_size,
    )

    sel = res["selected_candidates"]
    n_clusters = len({b["pose_cluster"] for b in res.get("passing_backbones", [])})
    print(f"\n================ [黃金凹槽骨架採收 v2 完成: {args.job_id}] ================")
    print(f"  通過 3D 凹槽/物理/零碰撞審計骨架數 : {len(res.get('passing_backbones', []))} / {len(res['audited_backbones'])} (共 {n_clusters} 個 3D 姿態簇)")
    print(f"  SolubleMPNN 總生成候選序列數       : {len(res['evaluated_candidates'])}")
    print(f"  入選跨簇零硬傷高多樣性序列數       : {len(sel)} (已冪等更新至 {reg_csv})")
    print(f"  完整審計報表                       : {res['all_csv']}")
    print(f"  FASTA 檔案                         : {res['pass_fasta']}")
    print(f"  AlphaFold 3 Server 批量 JSON       :")
    for fp in res["af3_files"]:
        print(f"    -> {fp}")
    if sel:
        print("\n--- Top 20 優先送檢序列（Part 1 JSON：優先覆蓋所有不同 3D 姿態簇、混合 α/β 與 B317 錨點！）---")
        for c in sel[:20]:
            ek_pct = (float(c["E_ratio"]) + float(c["K_ratio"])) * 100.0
            ala_pct = float(c["ala_ratio"]) * 100.0
            print(
                f"  #{c['af3_priority_rank']:02d} [{c['pose_cluster']}] {c['design_id']:50s} "
                f"bbQ={c['bb_quality_score']:.2f} ipTM={c['traj_i_pTM']:.2f} SC={c['traj_SC']:.2f} "
                f"Unsat={c['traj_unsat_hbonds']:3.0f} β={c['binder_strand_pct']:2.0f}% "
                f"scR={str(c['monomer_sc_rmsd'] or '-'):>5s} E+K={ek_pct:4.1f}% A={ala_pct:4.1f}% "
                f"hs={c['hotspot_hits']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
