# Changjiang Cup 2026 — AI Protein Design Pipeline (`AGENTS.md`)

## Project Overview

This repository houses the CS team's automated protein design, validation, and evaluation pipeline for **The Changjiang Cup 2026 Global AI Protein Design Challenge (长江杯 2026 全球 AI 蛋白质设计挑战赛)**.

### Security & Competition Guardrails (CRITICAL)

1. **Never commit generated candidate sequences or experimental scores to Git.**
   - Under Changjiang Cup rules, duplicate sequences across teams are awarded to the **earliest submitter**.
   - All generated candidate sequences (`.fasta`, `candidates.csv`, `design_registry.csv`), model outputs (`.pdb`, `.zip`), and job logs MUST be written strictly to the shared Google Drive workspace (`https://drive.google.com/drive/folders/1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg?usp=drive_link`, resolved on Colab via `/content/drive/.shortcut-targets-by-id/1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg/`).
2. **Enforce Pre-Run and Post-Run Quality Gates.**
   - Before running any generation batch, run `scripts/validate_input.py` to verify target SHA-256 hashes, hotspot chain/residue existence, and atomic completeness.
   - Inverse folding MUST use **SolubleMPNN** (SOP §C-4).
   - Generated sequences MUST pass `scripts/check_expressibility.py` (20 standard amino acids, length limits, even Cys count, charge/hydrophobic checks) and `scripts/check_novelty.py`.
3. **Do not use files under `_superseded/`.**
   - Always use the latest verified target structures (e.g., `targets/target1_EGFR/targets/t1_4KRP_domIII.pdb` for Target 1).

## Agent skills

### Issue tracker

Issues are tracked in GitHub Issues on `David-Liv/changjiang-protein-2026` via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

Uses the five default canonical triage labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context layout (`GLOSSARY.md` at repo root and ADRs in `docs/adr/`). See `docs/agents/domain.md`.
