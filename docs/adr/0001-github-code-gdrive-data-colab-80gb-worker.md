# ADR-0001: Public GitHub for Pipeline Code + Google Drive Job Queue for 80GB Colab Execution

- **Status**: Accepted
- **Date**: 2026-11-06

## Context

1. The team has access to a single **80GB GPU Colab instance** (under David's `cuhktu` account) alongside several CS and Biology teammates.
2. Google Colab binds GPU allocation to the logged-in user executing the notebook; simply sharing a `.ipynb` link does not share David's 80GB GPU with teammates.
3. Under Changjiang Cup 2026 rules, identical or near-identical sequences submitted by multiple teams are credited to the **earliest submitter**, meaning generated candidate sequences must remain private to the team.

## Decision

1. **Public GitHub Repository (`David-Liv/changjiang-protein-2026`) for Code & Public Target Inputs Only**:
   - All pipeline scripts, validation gates, job templates, and trimmed public PDB structures live in the public GitHub repository so any Colab instance or teammate can run `git clone` / `git pull` without managing tokens.
   - `.gitignore` strictly blocks any candidate sequences, output `.fasta`/`.csv` tables, and run artifacts from entering Git.
2. **Shared Google Drive (`https://drive.google.com/drive/folders/1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg?usp=drive_link`) as Private Data Store & Job Queue**:
   - Teammates submit runs by placing JSON job definitions into `jobs/pending/` on the shared Google Drive folder.
   - David's 80GB Colab runs `notebooks/colab_80gb_worker.ipynb`, which pulls the latest code from GitHub, polls `jobs/pending/`, runs pre-validation (`scripts/validate_input.py`), executes generation/evaluation (`BindCraft` + `SolubleMPNN` / `ColabFold`), runs post-validation (`scripts/check_expressibility.py`, `scripts/check_novelty.py`), and writes all outputs and `registry/design_registry.csv` back to Google Drive.
3. **Initial Tooling Scope**:
   - Phase 1 standardizes on **BindCraft (with SolubleMPNN) + ColabFold + SOP Python Validators**, prioritizing **Target 5 (`TL1A`, Minibinder)** and **Target 1 (`EGFR`, VHH)**.
