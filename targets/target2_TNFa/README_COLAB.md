# Colab 執行包（Ultra / A100）

```
colab_bundle/
├── TNF_ipSAE_colab.ipynb     ← 在 Colab 開這個
├── tnf_pipeline.py           ← 整條管線（notebook 內也有一份內嵌副本）
├── spec_stage1_mono.tsv      ← 單體決定性測試 3 項
├── spec_stage1_trimer.tsv    ← 三聚體 2 項
├── spec_stage23.tsv          ← 陽性 4 條 + 陰性 4 條
├── requirements_colab.txt
└── README_COLAB.md
```

## 三種用法，挑一種

**A. 什麼都不上傳（最省事，建議先這樣試）**
把 `TNF_ipSAE_colab.ipynb` 拖到 [colab.research.google.com](https://colab.research.google.com)
（或 `檔案 → 上傳筆記本`）。notebook 第 1 格會自己把 `tnf_pipeline.py` 與三個 spec
寫到 `/content`，不需要上傳任何其他檔案。

**B. 上傳整包（想改 spec 加自己的序列時）**
把 `colab_bundle.zip` 上傳到 Colab 再解開：

```bash
from google.colab import files; files.upload()     # 選 colab_bundle.zip
!unzip -o colab_bundle.zip && cd colab_bundle
```

第 1 格偵測到檔案已存在就會**沿用你的版本**，不會覆蓋你改過的 spec；
`tnf_pipeline.py` 若與內嵌版不同會印出兩邊的 md5 提醒你。

**C. 放 Google Drive（要保留 runs/ 跨 session）**
把這個資料夾複製到 `MyDrive/tnf_pipeline/`，然後跑 notebook 第 0 節的「掛 Drive」格
（把註解解開）。Colab runtime 斷線會清掉 `/content`，但 Drive 不會。

## 執行階段設定

`執行階段 → 變更執行階段類型 → GPU → A100`（Ultra 可選 A100 / H100）。

第 0 格會印出實際拿到的 GPU 與顯存，並把建議的 `--chunk` 寫進 `chunk.txt`
（≥ 38 GB → `0`，代表不分塊、最快；否則 `64`）。後面每一格都從這個檔案讀，
所以**不要跳過第 0 格**；斷線重連後也只需重跑第 0 格。

若拿到的是 T4（16 GB），第 0 格會警告，且第 6 節三聚體很可能 OOM ——
`fold` 偵測到 OOM 會自動用 chunk/4 重試一次，仍失敗就手動加 `--chunk 16` 或跳過。

## 各節在做什麼

| 節 | 內容 | GPU | 粗略時間（A100） |
|---|---|---|---|
| 0 | 環境、GPU 偵測、寫 `chunk.txt` | — | < 1 分 |
| 1 | 寫出腳本與 spec | — | 秒 |
| 1a | 9 項自測（計分、PDB 清理、triage） | — | 秒 |
| 1b | 序列 triage 預檢 | — | 秒 |
| 2 | 冒煙測試（TNF 單體 157 aa） | 需 | 首次含下載 2.6 GB 權重 |
| 3 | 單體決定性比較（最優先的結果） | 需 | 數分鐘 |
| 4 | 陽性／陰性校準（9 條） | 需 | 約十分鐘 |
| 5 | 確定性驗證 + 8 變體穩定度（72 次推理） | 需 | 數十分鐘 |
| 6 | 三聚體（~700 token） | 需 | 每項數分鐘 |
| 7 | 與 5M2J / 1TNF / 2AZ5 疊合 | — | 秒 |
| 8 | 一鍵全跑 | 需 | — |
| 9 | 打包下載 | — | 秒 |

時間是數量級估計，不是實測 —— 這套管線還沒在 GPU 上跑過。

## 判讀三個重點

1. **陽性對照全部低分時，先看表二不是調 cutoff。** 接觸數 0 或 min 重原子 > 5 Å
   表示兩條鏈根本沒對接上（ESMFold 是單序列模型，抗體–抗原界面常常失敗）。
   接觸數正常但 ipSAE 低，才是 cutoff 的問題，去比 `max@15` 欄。

2. **三聚體要看「各鏈對明細」與 `合併A-ALL`。** VHH3 跨 B、C 兩個亞基結合，
   拆成鏈對後每一半的 n 變小、d0 懲罰過重。`score` 已自動取所有 binder–target
   鏈對的最大值；`合併A-ALL` 把 B+C+D 當單一靶點，只作為「低估多少」的診斷，
   **已偏離 Dunbrack 原版定義，不可當成官方分數**。

3. **標準差不是來自 random seed。** ESMFold 在 `eval()` 下是確定性的，
   同序列跑不同 seed 會給出完全相同的 PAE。第 5 節的變異來自
   `--recycles` / `--linkers` / `--chain-orders`，第一格會當場驗證這件事。

## 為什麼是單一檔案

Colab 上單檔比多模組好，三個具體原因：

1. **`!python tnf_pipeline.py` 每次都是新 process**，不會有 `import` 快取問題。
   多模組時你改了 `ipsae.py`，已 import 的舊版還留在 kernel 記憶體裡，
   要 `%autoreload` 或重啟 runtime 才生效 —— 這個坑很安靜，分數會悄悄用舊邏輯算。
2. **斷線重連後要重建環境**，複製一個檔案比六個快，`%cd` 與 `sys.path` 也不用管。
3. **notebook 能把它整個內嵌**，所以用法 A 可以零上傳執行。六個檔案的 heredoc
   也能內嵌，但第 1 格會變成六段，出錯時難定位。

多檔只在「多人用 git 各改一塊」時比較好，那不是 Colab 的工作方式。
整合前的版本留在上層 `legacy/`，僅供對照。
