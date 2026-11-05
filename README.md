# 長江杯 2026 全球 AI 蛋白質設計挑戰賽 — CUHK 團隊 Pipeline

本倉庫為香港中文大學（CUHK）參賽團隊的自動化蛋白質設計、質量門禁與評估流水線（基於 **80GB Colab GPU + GitHub 代碼同步 + Google Drive 任務隊列**）。

---

## 1. 一鍵在 80GB Colab 啟動（兩種入口）

| 入口 Notebook | 適用場景 | 一鍵打開 Colab 連結 |
| :--- | :--- | :--- |
| **`colab_80gb_worker.ipynb`**<br/>（全靶點隊列總控端） | 自動拉取 GitHub 最新代碼、監聽網盤 `Changjiang_Cup_2026/jobs/pending/` 任務單（支持靶點 2 校準、靶點 5 Minibinder 生成、靶點 1 EGFR 生成）並自動串聯 `ipSAE` 評估與寫回總表 | [Open in Colab (`colab_80gb_worker.ipynb`)](https://colab.research.google.com/github/David-Liv/changjiang-protein-2026/blob/main/notebooks/colab_80gb_worker.ipynb) |
| **`TNF_ipSAE_colab.ipynb`**<br/>（靶點 2 TNF-α 專項校準） | 生物組定制的 TNF-α 單體 vs 三聚體、正負對照、72 次變體穩定度與晶體疊合（5M2J/1TNF/2AZ5）測試，跑完自動打包 `tnf_calibration_runs.zip` 存入共享網盤 `Changjiang_Cup_2026/runs/target2_TNFa_calibration/` | [Open in Colab (`TNF_ipSAE_colab.ipynb`)](https://colab.research.google.com/github/David-Liv/changjiang-protein-2026/blob/main/notebooks/TNF_ipSAE_colab.ipynb) |

---

## 2. 主攻賽道與靶點狀態

| 靶點編號 | 靶點名稱 | 賽道分子類型 | 長度要求 | 特殊要求 | 當前狀態與輸入路徑 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Target 1** | **EGFR** | **VHH（奈米抗體）** | 100–150 aa | **pH 6.0 親和力測定** | ✅ 生物組已交付 v2 定版 [`targets/target1_EGFR/`](./targets/target1_EGFR/README.md)（預設使用 `t1_4KRP_domIII.pdb`） |
| **Target 2** | **TNF-α** | **VHH（奈米抗體）** | 100–150 aa | 無特殊附加條件 | ✅ 已接入單體/三聚體 `ESMFold + ipSAE + pDockQ` 校準管線 [`targets/target2_TNFa/`](./targets/target2_TNFa/README_COLAB.md) |
| **Target 5** | **TL1A** | **Minibinder（迷你結合蛋白）** | 40–100 aa | 無特殊附加條件 | ✅ 已準備單體零殘缺結構 [`targets/target5_TL1A/`](./targets/target5_TL1A/README.md)（`t1_2RE9_chainA.pdb`） |

---

## 3. 核心腳本與全自動門禁（Pre-run & Post-run Gates）

1. **跑批前輸入驗證（[`scripts/validate_input.py`](./scripts/validate_input.py)）**：
   驗證靶點 PDB SHA-256 雜湊、熱點殘基存在性、跨鏈殘基索引唯一性及主鏈完整性。
2. **序列生成後複合體折疊與 `ipSAE` / `pDockQ` 評分（[`scripts/tnf_pipeline.py`](./scripts/tnf_pipeline.py)）**：
   - `selftest`：9 項數學與 PDB 分鏈清理單元測試（無需 GPU）。
   - `triage`：長度（≤250 aa）、半胱氨酸（Cys）奇偶性、淨電荷、GRAVY 疏水度與連續疏水段預檢。
   - `fold` & `score`：調用 HuggingFace `ESMFold` 預測複合體，向量化計算 **Dunbrack `ipSAE`（含多聚體跨鏈合併修正 `A-ALL`）** 與 **Bryant `pDockQ` / 最小重原子距離**。
   - `compare`：自動從 RCSB 下載晶體結構並以 Kabsch 算法計算表位區 CA RMSD。
3. **表達性與新穎性硬閘（[`scripts/check_expressibility.py`](./scripts/check_expressibility.py) & [`scripts/check_novelty.py`](./scripts/check_novelty.py)）**：
   計算 E/K 比例、pH 6.0 / pH 7.4 淨電荷、以及與已知陽性結合蛋白的全長/CDR3 序列同一性。
