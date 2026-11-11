#!/usr/bin/env python3
"""Export passing candidates from design_registry.csv as AlphaFold 3 Server uploadable JSON files.

Example:
  python3 scripts/export_af3_jobs.py --registry /path/to/design_registry.csv \
      --target-id target5_TL1A --exclude-job job_t5_tl1a_minibinder_001 \
      --out-dir /path/to/af3_batches
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "pipeline") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pipeline"))

import worker  # noqa: E402

WT_TL1A_164 = (
    "LRADGDKPRAHLTVVRQTPTQHFKNQFPALHWEHELGLAFTKNRMNYTNKFLLIPESGDYFIYSQVTFRGMTSECSEIR"
    "QAGRPNKPDSITVVITKVTDSYPEPTQLLMGTKSVCEVGSNWFQPIYLGAMFSLQEGDKLMVNVSDISLVDYTKEDKTFFGAFLL"
)


def main() -> int:
    p = argparse.ArgumentParser(description="Export AlphaFold 3 Server JSON batches from design_registry.csv")
    p.add_argument("--registry", required=True, help="Path to design_registry.csv")
    p.add_argument("--target-id", default="target5_TL1A", help="Filter by target_id")
    p.add_argument("--exclude-job", action="append", default=["job_t5_tl1a_minibinder_001"], help="Job IDs to exclude")
    p.add_argument("--target-seq", default=WT_TL1A_164, help="Wild-type target monomer sequence")
    p.add_argument("--target-copies", type=int, default=3, help="Copies of target chain in AF3 (3 for TL1A homotrimer)")
    p.add_argument("--batch-size", type=int, default=20, help="Max jobs per AF3 Server JSON file (default: 20)")
    p.add_argument("--out-dir", required=True, help="Directory to write AF3 JSON files")
    args = p.parse_args()

    reg_path = pathlib.Path(args.registry)
    rows = list(csv.DictReader(reg_path.open("r", newline="", encoding="utf-8")))
    excl = set(args.exclude_job or [])
    filtered = [
        r for r in rows
        if r.get("target_id") == args.target_id
        and r.get("job_id") not in excl
        and r.get("expressibility_verdict") == "PASS"
    ]
    job_meta = {
        "job_id": f"{args.target_id}_selected",
        "af3_target_sequence": args.target_seq,
        "af3_target_copies": args.target_copies,
    }
    written = worker.export_alphafold3_server_jobs(
        filtered, job_meta, pathlib.Path(args.out_dir), batch_size=args.batch_size
    )
    print(f"Exported {len(filtered)} PASS candidates into {len(written)} AlphaFold 3 Server JSON file(s):")
    for fp in written:
        print(f"  - {fp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
