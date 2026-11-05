# EGFR × VHH 靶點檔案 — 交接包

長江杯 2026 靶點一 ｜ 生物組 → 計算組 ｜ 2026-11-01 ｜ 版本 v2

## 先讀這個

1. **`docs/EGFR_VHH_target_dossier.md`** — 靶點檔案本體。表位決策、熱點字串、禁區、配額建議全在裡面。
2. **`targets/`** — 五個修剪 PDB + `SHA256SUMS.txt`。**預設用 `t1_4KRP_domIII.pdb`。**
3. 跑批前執行 `scripts/validate_input.py`（SOP §B-1 五項斷言）。

## 開跑前的兩個前提

**① 一個送出硬閘尚未解除。** 新穎性判準需主辦方裁示：三條已發表的 EGFR VHH
互相之間全長同一性為 76.6–81.1%，若按「全長同一性 ≤ 75%」判，它們會互相出局。
建議改判 CDR3（不相關結合物實測 35–53%，門檻取 ≤50%）。
**這是送出閘不是開跑閘——現在可以生成，裁示回來前不要送出。**
檢查工具：`scripts/check_novelty.py`，排除集：`targets/known_binders_SEQRES.fasta`。

**② 反向封堵指令有歧義，執行前須確認。** 生物組內部有一條表述為
「反向封堵 Face A1，只用 cetux 專有的 F436/V441/I462/K489/I491」。
這五個 UniProt 殘基換算成熟編號為 412/417/438/465/467，**全部屬於 Pole N 主力熱點組，
沒有一個屬於 Face A1（7D12 面）**。封堵與瞄準是互為反向的操作，照字面執行會把主力熱點封死。
**確認語意前不要執行任何負向設計。**

## 目錄

```
targets/      五個修剪 PDB、SHA256SUMS.txt、新穎性排除集   ← 跑批輸入
docs/         靶點檔案、設計登錄表範本、ChimeraX 檢視腳本
scripts/      輸入驗證器、組成檢查、新穎性檢查、PDB 解析庫
evidence/     六個母結構 PDB、UniProt 序列、全部日誌與中間數值
_superseded/  第一代產物，不得使用（保留供追溯）
```

## 修剪檔速查

| 檔案 | 母結構 | 解析度 | 成熟範圍 | 殘基 | 用途 |
|---|---|---|---|---|---|
| **t1_4KRP_domIII.pdb** ★ | 4KRP | 2.82 Å | 311–514 | 204 | **預設**。零序列錯配 |
| t2_4KRL_domIII.pdb | 4KRL | 2.85 Å | 311–511 | 201 | pH 6.0 晶長條件，交叉模板 |
| t3_4KRO_domIII.pdb | 4KRO | 3.05 Å | 311–514 | 204 | 第三個獨立來源 |
| t4_4KRP_faceB.pdb | 4KRP | 2.82 Å | 200–612 | 405 | Face B 配額專用 |
| t5_4KRP_fullECD.pdb | 4KRP | 2.82 Å | 4–612 | 534 | 最終複合物脈絡檢查 |

全部零切斷雙硫鍵，全部對 P00533-1 零序列錯配，編號為成熟編號（UniProt = 成熟 + 24）。

## 熱點字串（鏈 ID 統一為 A）

```
組 1  Face A2 完整（主力，含 pH 槓桿 His409）
  BindCraft   : A409,A412,A417,A438,A465,A467
  RFdiffusion : ppi.hotspot_res=[A409,A412,A417,A438,A465,A467]

組 2  Face A2 緊緻 β 核心
  BindCraft   : A412,A417,A438,A465,A467
  RFdiffusion : ppi.hotspot_res=[A412,A417,A438,A465,A467]

組 3  Face B（分散押注，配 t4）
  BindCraft   : A310,A400,A401,A403,A405,A431
  RFdiffusion : ppi.hotspot_res=[A310,A400,A401,A403,A405,A431]
```

## 跑批前檢查

```bash
python3 scripts/validate_input.py \
  --target targets/t1_4KRP_domIII.pdb \
  --hotspots A409,A412,A417,A438,A465,A467 \
  --sha 950b0cc8884cd44913edfc50653af91f4cc86b9799fa72d3923f957fa07f5be7
```

產出序列後：
```bash
python3 scripts/check_expressibility.py designs.fasta     # SOP §C-2 組成檢查
python3 scripts/check_novelty.py targets/known_binders_SEQRES.fasta   # 新穎性基準
```

反向摺疊最後一步**必須**是 SolubleMPNN（SOP §C-4）。
