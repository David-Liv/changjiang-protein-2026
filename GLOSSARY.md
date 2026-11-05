# Glossary

Canonical domain vocabulary for the Changjiang Cup 2026 protein design pipeline.

## Core Terms

- **Target (靶点)**: The biological protein we are designing a binder against. In Changjiang Cup 2026, there are 7 official targets; our team prioritizes **Target 1 (`EGFR`, VHH track, pH 6.0)**, **Target 5 (`TL1A`, Minibinder track)**, and **Target 2 (`TNF-alpha`, VHH track)**.
- **Trimmed PDB (修剪结构)**: A cleaned, cropped `.pdb` coordinate file prepared in Stage A (e.g., `t1_4KRP_domIII.pdb`) with no partial residues, verified against UniProt reference sequences, and pinned by a SHA-256 hash in `SHA256SUMS.txt`.
- **Hotspot (热点残基)**: Specific surface residues on the Target (e.g., `A409,A412,A417,A438,A465,A467`) passed to the generator to steer where the binder docks.
- **Minibinder (迷你结合蛋白)**: A *de novo* single-chain small binding protein of **40–100 amino acids** (20 standard amino acids only, no tags/signal peptides). Used for **Target 5 (`TL1A`)** and **Target 7 (`IL-7Rα`)**.
- **VHH / Nanobody (纳米抗体)**: A single-domain heavy-chain antibody variable region of **100–150 amino acids** containing complete framework regions (`FR1–FR4`) and loops (`CDR1–CDR3`). Used for **Target 1 (`EGFR`)**, **Target 2 (`TNF-α`)**, **Target 3 (`GLP-1R`)**, **Target 4 (`TfR1`)**, and **Target 6 (`TNFR1`)**. *Never submit a plain helical Minibinder to a VHH track.*
- **SolubleMPNN**: Solubility-biased inverse folding weights for ProteinMPNN. Required by SOP §C-4 as the final sequence design step to maximize cell-free/E. coli expression success.
- **Job Spec (任务单)**: A JSON file dropped into `/content/drive/MyDrive/Changjiang_Cup_2026/jobs/pending/` specifying the target PDB, expected SHA-256, hotspots, length bounds, and number of trajectories for the 80GB Colab Worker to execute.
- **Design Registry (设计登录表)**: The single source of truth CSV (`registry/design_registry.csv` on Google Drive) recording every generated candidate, its provenance (`target_sha256`, generator, seed, hotspots, `is_solubleMPNN`), folding metrics (`ipTM`, `pAE`), and expressibility/novelty gate verdicts.
