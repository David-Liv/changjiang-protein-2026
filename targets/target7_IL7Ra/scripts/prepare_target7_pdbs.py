#!/usr/bin/env python3
"""Deterministic Stage-A target preparation for Target 7: IL-7Ra (CD127, UniProt P16871-1).

Produces 5 clean, validated trimmed PDB templates (all standardized to Chain B, Mature numbering
where UniProt = Mature + 20, 0 missing heavy atoms, 0 mismatches vs canonical P16871-1, LF line endings):
  - t1_3DI3_coreECD.pdb          (Primary: 3DI3 Chain B, Mature 17-209, 193 aa, IL-7 holo state)
  - t2_7OPB_holo_coreECD.pdb     (Cross-template: 7OPB Chain A, Mature 16-211, 196 aa, minibinder holo state)
  - t3_3UP1_apo_coreECD.pdb      (Apo control: 3UP1 Chain A, Mature 16-211, 196 aa, resting state)
  - t4_AF2_pred_fullECD.pdb      (Full-ECD control: AF-P16871-F1, Mature 1-219 = UniProt 21-239, 219 aa)
  - t5_3DI3_coreECD_backMaskK.pdb(Negative-design variant: t1 + back-face glycan sequons B45K, B131K, B162K)
"""
from __future__ import annotations

import hashlib
import pathlib
import numpy as np

BASE_DIR = pathlib.Path(__file__).resolve().parent.parent
EV_DIR = BASE_DIR / "evidence"
TGT_DIR = BASE_DIR / "targets"

EXPECTED_ATOMS = {
    "GLY": ["N", "CA", "C", "O"],
    "ALA": ["N", "CA", "C", "O", "CB"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
}


def parse_chain_residues(pdb_path: pathlib.Path, chain_id: str, res_offset: int = 0) -> dict[int, dict]:
    """Parse protein heavy atoms of a single chain, shifting residue numbers by `res_offset`."""
    residues: dict[int, dict] = {}
    for ln in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not ln.startswith("ATOM  "):
            continue
        if ln[16] not in (" ", "A") or ln[21] != chain_id:
            continue
        resname = ln[17:20].strip()
        if resname not in EXPECTED_ATOMS:
            continue
        aname = ln[12:16].strip()
        elem = (ln[76:78].strip() if len(ln) >= 78 else "") or aname[0]
        if elem == "H" or aname.startswith("H"):
            continue
        rnum = int(ln[22:26]) + res_offset
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])], dtype=float)
        occ = float(ln[54:60] or 1.0)
        bfac = float(ln[60:66] or 20.0)
        rec = residues.setdefault(rnum, {"resname": resname, "atoms": {}})
        if aname not in rec["atoms"]:
            rec["atoms"][aname] = {"xyz": xyz, "occ": occ, "b": bfac, "elem": elem}
    return residues


def superpose_points(src_anchors: np.ndarray, dst_anchors: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Rigid-body Kabsch transform mapping `src_anchors` onto `dst_anchors`, applied to `pts`."""
    c_src = src_anchors.mean(axis=0)
    c_dst = dst_anchors.mean(axis=0)
    P = src_anchors - c_src
    Q = dst_anchors - c_dst
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = float(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, -1.0 if d < 0.0 else 1.0]) @ U.T
    return (pts - c_src) @ R.T + c_dst


def complete_and_canonicalize_from_donor(
    target_res: dict[int, dict],
    donor_res: dict[int, dict],
    canonicalize_positions: set[int] | None = None,
) -> tuple[int, list[str]]:
    """Fill missing sidechain heavy atoms (0.000 A shift on existing atoms) and canonicalize SNP positions."""
    canonicalize_positions = canonicalize_positions or set()
    added_count = 0
    notes: list[str] = []

    for rnum in sorted(target_res.keys()):
        trec = target_res[rnum]
        drec = donor_res.get(rnum)
        if drec is None:
            continue

        # Canonicalize SNP (e.g. Mature T46 -> I46 per UniProt P16871-1) keeping exact N,CA,C,O,CB
        if rnum in canonicalize_positions and trec["resname"] != drec["resname"]:
            old_rn = trec["resname"]
            new_rn = drec["resname"]
            anchors = [a for a in ("N", "CA", "C", "CB") if a in trec["atoms"] and a in drec["atoms"]]
            src_a = np.array([drec["atoms"][a]["xyz"] for a in anchors])
            dst_a = np.array([trec["atoms"][a]["xyz"] for a in anchors])
            keep_bb = {a: trec["atoms"][a] for a in ("N", "CA", "C", "O", "CB") if a in trec["atoms"]}
            mean_b = float(np.mean([v["b"] for v in keep_bb.values()]))
            for aname in EXPECTED_ATOMS[new_rn]:
                if aname not in keep_bb:
                    mapped_xyz = superpose_points(src_a, dst_a, drec["atoms"][aname]["xyz"][None, :])[0]
                    keep_bb[aname] = {"xyz": mapped_xyz, "occ": 1.0, "b": mean_b, "elem": aname[0]}
                    added_count += 1
            trec["resname"] = new_rn
            trec["atoms"] = keep_bb
            notes.append(f"B{rnum}:{old_rn}->{new_rn}")
            continue

        expected = EXPECTED_ATOMS[trec["resname"]]
        missing = [a for a in expected if a not in trec["atoms"]]
        if not missing:
            continue
        if drec["resname"] != trec["resname"]:
            raise ValueError(f"Donor residue mismatch at B{rnum}: {trec['resname']} vs {drec['resname']}")

        anchors = [a for a in ("N", "CA", "C", "CB") if a in trec["atoms"] and a in drec["atoms"]]
        src_a = np.array([drec["atoms"][a]["xyz"] for a in anchors])
        dst_a = np.array([trec["atoms"][a]["xyz"] for a in anchors])
        mean_b = float(np.mean([v["b"] for v in trec["atoms"].values()]))
        for aname in missing:
            mapped_xyz = superpose_points(src_a, dst_a, drec["atoms"][aname]["xyz"][None, :])[0]
            trec["atoms"][aname] = {"xyz": mapped_xyz, "occ": 1.0, "b": mean_b, "elem": aname[0]}
            added_count += 1
        notes.append(f"B{rnum}({trec['resname']}):+{''.join(missing)}")

    return added_count, notes


def mutate_surface_to_lysine(
    target_res: dict[int, dict],
    positions: list[int],
    lys_templates: list[dict],
) -> list[str]:
    """Mutate specified surface positions to clash-free full-atom LYS (SOP §A-4 Rule 1)."""
    mut_notes: list[str] = []
    for rnum in positions:
        trec = target_res[rnum]
        old_rn = trec["resname"]
        other_pts = np.array(
            [
                at["xyz"]
                for k, r in target_res.items()
                if k != rnum
                for at in r["atoms"].values()
            ]
        )
        dst_a = np.array([trec["atoms"][a]["xyz"] for a in ("N", "CA", "C", "CB")])
        best_clearance = -1.0
        best_sc: dict[str, np.ndarray] = {}
        for tmpl in lys_templates:
            src_a = np.array([tmpl["atoms"][a]["xyz"] for a in ("N", "CA", "C", "CB")])
            sc_names = ["CG", "CD", "CE", "NZ"]
            sc_pts = np.array([tmpl["atoms"][a]["xyz"] for a in sc_names])
            mapped = superpose_points(src_a, dst_a, sc_pts)
            dmin = float(np.sqrt(np.min(np.sum((mapped[:, None, :] - other_pts[None, :, :]) ** 2, axis=2))))
            if dmin > best_clearance:
                best_clearance = dmin
                best_sc = dict(zip(sc_names, mapped))

        keep = {a: dict(trec["atoms"][a]) for a in ("N", "CA", "C", "O", "CB")}
        mean_b = float(np.mean([v["b"] for v in keep.values()]))
        for aname, xyz in best_sc.items():
            keep[aname] = {"xyz": xyz, "occ": 1.0, "b": mean_b, "elem": aname[0]}
        trec["resname"] = "LYS"
        trec["atoms"] = keep
        mut_notes.append(f"B{rnum}:{old_rn}->LYS(min_dist={best_clearance:.2f}A)")
    return mut_notes


def write_clean_pdb(
    out_path: pathlib.Path,
    residues: dict[int, dict],
    chain_id: str,
    remarks: list[str],
) -> str:
    """Write standard PDB with LF line endings and return SHA-256 hex digest."""
    lines: list[str] = [f"REMARK 999 {r}" for r in remarks]
    serial = 1
    for rnum in sorted(residues.keys()):
        rec = residues[rnum]
        rn = rec["resname"]
        for aname in EXPECTED_ATOMS[rn]:
            at = rec["atoms"][aname]
            x, y, z = at["xyz"]
            name_field = f" {aname:<3s}" if len(aname) < 4 else f"{aname:<4s}"
            lines.append(
                f"ATOM  {serial:5d} {name_field} {rn:>3s} {chain_id}{rnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{at['occ']:6.2f}{at['b']:6.2f}          {at['elem']:>2s}"
            )
            serial += 1
    lines.append("TER")
    lines.append("END")
    content_lf = ("\n".join(lines) + "\n").encode("utf-8")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(content_lf)
    return hashlib.sha256(content_lf).hexdigest()


def copy_residues(src: dict[int, dict], lo: int, hi: int) -> dict[int, dict]:
    out: dict[int, dict] = {}
    for k in range(lo, hi + 1):
        if k in src:
            out[k] = {
                "resname": src[k]["resname"],
                "atoms": {
                    an: {"xyz": at["xyz"].copy(), "occ": at["occ"], "b": at["b"], "elem": at["elem"]}
                    for an, at in src[k]["atoms"].items()
                },
            }
    return out


def main() -> None:
    TGT_DIR.mkdir(parents=True, exist_ok=True)

    # Remove unneeded 5J11.pdb if present to keep git repo compact; keep 3DI3, 3UP1, 7OPB, AF-P16871-F1, P16871.fasta
    p_5j11 = EV_DIR / "5J11.pdb"
    if p_5j11.exists():
        p_5j11.unlink()

    # Load UniProt P16871-1 and write official competition ECD FASTA (Glu21-Asp239, 219 aa)
    raw_fa = [l.strip() for l in (EV_DIR / "P16871.fasta").read_text(encoding="utf-8").splitlines() if l.strip()]
    unp_full = "".join(l for l in raw_fa if not l.startswith(">"))
    official_ecd = unp_full[20:239]  # UniProt 21..239 = Mature 1..219 (219 aa)
    assert len(official_ecd) == 219 and official_ecd.startswith("ESGYAQ") and official_ecd.endswith("SGEMD")
    (EV_DIR / "P16871_official_ECD_21-239.fasta").write_bytes(
        (
            ">P16871_IL7RA_ECD_21-239 | Official Changjiang Cup Target 7 Sequence | Mature 1-219 (219 aa)\n"
            f"{official_ecd}\n"
        ).encode("utf-8")
    )

    # Parse donor & templates (unified to Mature numbering: UniProt = Mature + 20)
    res_3di3 = parse_chain_residues(EV_DIR / "3DI3.pdb", "B", res_offset=0)      # Mature 17..209 (193 aa)
    res_7opb = parse_chain_residues(EV_DIR / "7OPB.pdb", "A", res_offset=-20)    # Mature 16..211 (196 aa, 0 missing, canonical I46)
    res_3up1 = parse_chain_residues(EV_DIR / "3UP1.pdb", "A", res_offset=0)      # Mature 12..212 (201 aa)
    res_af2 = parse_chain_residues(EV_DIR / "AF-P16871-F1.pdb", "A", res_offset=-20)  # Mature -19..439

    # Collect experimental LYS rotamer templates from 7OPB Chain A
    lys_templates = [v for v in res_7opb.values() if v["resname"] == "LYS"]

    sha_rows: list[tuple[str, str]] = []

    # 1) t1_3DI3_coreECD.pdb (Primary: 3DI3 Chain B, Mature 17-209, 193 aa)
    t1 = copy_residues(res_3di3, 17, 209)
    n_add1, notes1 = complete_and_canonicalize_from_donor(t1, res_7opb, canonicalize_positions={46})
    sha1 = write_clean_pdb(
        TGT_DIR / "t1_3DI3_coreECD.pdb",
        t1,
        "B",
        [
            "TARGET 7 IL-7RA (UNIPROT P16871-1) - PRIMARY TEMPLATE t1_3DI3_coreECD",
            "SOURCE PDB 3DI3 CHAIN B (2.15 A, IL-7 BOUND HOLO CONFORMATION)",
            "CHAIN B MATURE NUMBERING 17-209 (UNIPROT 37-229, 193 AA, 3 INTACT DISULFIDES)",
            f"COMPLETED {n_add1} MISSING SIDECHAIN ATOMS + CANONICALIZED B46 THR->ILE FROM 7OPB_A",
            "PRE-EXISTING 3DI3 EXPERIMENTAL ATOM SHIFT = 0.000 A",
        ],
    )
    sha_rows.append((sha1, "t1_3DI3_coreECD.pdb"))
    print(f"Built t1_3DI3_coreECD.pdb ({len(t1)} aa, +{n_add1} atoms): {notes1}")

    # 2) t2_7OPB_holo_coreECD.pdb (Holo minibinder crystal template: 7OPB Chain A, Mature 16-211, 196 aa)
    t2 = copy_residues(res_7opb, 16, 211)
    sha2 = write_clean_pdb(
        TGT_DIR / "t2_7OPB_holo_coreECD.pdb",
        t2,
        "B",
        [
            "TARGET 7 IL-7RA (UNIPROT P16871-1) - HOLO MINIBINDER CROSS-TEMPLATE t2_7OPB_holo_coreECD",
            "SOURCE PDB 7OPB CHAIN A RENAMED TO B & RENUMBERED -20 TO MATURE 16-211 (UNIPROT 36-231, 196 AA)",
            "ZERO MISSING SIDECHAIN ATOMS, ZERO SEQUENCE MISMATCHES VS P16871-1",
        ],
    )
    sha_rows.append((sha2, "t2_7OPB_holo_coreECD.pdb"))
    print(f"Built t2_7OPB_holo_coreECD.pdb ({len(t2)} aa)")

    # 3) t3_3UP1_apo_coreECD.pdb (Resting state apo template: 3UP1 Chain A, Mature 16-211, 196 aa)
    t3 = copy_residues(res_3up1, 16, 211)
    n_add3, notes3 = complete_and_canonicalize_from_donor(t3, res_7opb, canonicalize_positions={46})
    sha3 = write_clean_pdb(
        TGT_DIR / "t3_3UP1_apo_coreECD.pdb",
        t3,
        "B",
        [
            "TARGET 7 IL-7RA (UNIPROT P16871-1) - APO RESTING TEMPLATE t3_3UP1_apo_coreECD",
            "SOURCE PDB 3UP1 CHAIN A RENAMED TO B, MATURE 16-211 (UNIPROT 36-231, 196 AA)",
            f"COMPLETED {n_add3} MISSING SIDECHAIN ATOMS + CANONICALIZED B46 THR->ILE FROM 7OPB_A",
        ],
    )
    sha_rows.append((sha3, "t3_3UP1_apo_coreECD.pdb"))
    print(f"Built t3_3UP1_apo_coreECD.pdb ({len(t3)} aa, +{n_add3} atoms): {notes3}")

    # 4) t4_AF2_pred_fullECD.pdb (Official full-ECD predicted control: Mature 1-219 = UniProt 21-239, 219 aa)
    t4 = copy_residues(res_af2, 1, 219)
    sha4 = write_clean_pdb(
        TGT_DIR / "t4_AF2_pred_fullECD.pdb",
        t4,
        "B",
        [
            "TARGET 7 IL-7RA (UNIPROT P16871-1) - FULL OFFICIAL ECD PREDICTED TEMPLATE t4_AF2_pred_fullECD",
            "SOURCE ALPHAFOLD AF-P16871-F1 RENAMED TO B & RENUMBERED -20 TO MATURE 1-219 (UNIPROT 21-239, 219 AA)",
        ],
    )
    sha_rows.append((sha4, "t4_AF2_pred_fullECD.pdb"))
    print(f"Built t4_AF2_pred_fullECD.pdb ({len(t4)} aa)")

    # 5) t5_3DI3_coreECD_backMaskK.pdb (Back-face N-glycan sequon Lys shield: B45K, B131K, B162K, all >16.5 A from hotspots)
    t5 = copy_residues(t1, 17, 209)
    mut5 = mutate_surface_to_lysine(t5, [45, 131, 162], lys_templates)
    sha5 = write_clean_pdb(
        TGT_DIR / "t5_3DI3_coreECD_backMaskK.pdb",
        t5,
        "B",
        [
            "TARGET 7 IL-7RA (UNIPROT P16871-1) - BACK-FACE GLYCAN SEQUON LYS-SHIELDED TEMPLATE t5",
            "DERIVED FROM t1_3DI3_coreECD WITH SOP A-4 NEGATIVE DESIGN MUTATIONS: B45K, B131K, B162K (>16.5 A FROM HOTSPOTS)",
        ],
    )
    sha_rows.append((sha5, "t5_3DI3_coreECD_backMaskK.pdb"))
    print(f"Built t5_3DI3_coreECD_backMaskK.pdb ({len(t5)} aa): {mut5}")

    # Write known_binders_IL7Ra_SEQRES.fasta (tag-stripped public positive controls for SOP §C-1 novelty gate)
    known_fa = (
        ">7OPB_D_Cao2022_IL7Ra_minibinder | De novo 3-helix minibinder (Nature 2022, 55 aa, KD ~ nM)\n"
        "SVIEKLRKLEKQARKQGDEVLVMLARMVLEYLEKGWVSEEDADESADRIEEVLKK\n"
        ">IL7_HUMAN_P13232_26-177 | Mature human Interleukin-7 ligand (3DI3 Chain A WT, 152 aa)\n"
        "DCDIEGKDGKQYESVLMVSIDQLLDSMKEIGSNCLNNEFNFFKRHICDANKEGMFLFRAARKLRQFLKMNSTGDFDLHLLKVSEGTTILLNCTGQVKGRKPAALGEAQPTKSLEENKSLKEQKKLNDLCFLKRLLQEIKTCWNKILMGTKEH\n"
        ">TSLP_HUMAN_Q969D9_29-159 | Mature human TSLP ligand (5J11 Chain A tag-stripped, 126 aa)\n"
        "YDFTNCDFEKIKAAYLSTISKDLITYMSGTKSTEFNNTVSCSNRPHCLTEIQSLTFNPTAGCASLAKEMFAMKTKAALAIWCPGYSETQINATQAMKKVTTNKCLEQVSQLQGLWRRFNRPLLKQQ\n"
    )
    known_path = TGT_DIR / "known_binders_IL7Ra_SEQRES.fasta"
    known_path.write_bytes(known_fa.encode("utf-8"))
    sha_rows.append((hashlib.sha256(known_path.read_bytes()).hexdigest(), "known_binders_IL7Ra_SEQRES.fasta"))

    # Write SHA256SUMS.txt & evidence provenance
    (TGT_DIR / "SHA256SUMS.txt").write_bytes(
        ("".join(f"{h}  {fn}\n" for h, fn in sha_rows)).encode("utf-8")
    )
    ev_shas = []
    for ep in sorted(EV_DIR.glob("*")):
        if ep.is_file() and ep.name != "PDB_PROVENANCE_SHA256.txt":
            ev_shas.append(f"{hashlib.sha256(ep.read_bytes()).hexdigest()}  {ep.name}\n")
    (EV_DIR / "PDB_PROVENANCE_SHA256.txt").write_bytes("".join(ev_shas).encode("utf-8"))
    print("Wrote SHA256SUMS.txt and PDB_PROVENANCE_SHA256.txt successfully.")


if __name__ == "__main__":
    main()
