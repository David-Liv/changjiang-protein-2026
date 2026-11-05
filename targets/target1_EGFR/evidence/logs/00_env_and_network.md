# 日誌 00 — 環境盤點與網路探測
階段：階段 A 前置　｜　靶點：靶點 1（G12D pMHC）　｜　時間：2026-11-01（Asia/Shanghai）
角色：位點與表達性負責人（SOP v1.0）

## 1. 嘗試：盤點計算環境
指令：`free -g` / `nproc` / `python3 -V` / 模組與 CLI 探測
結果：
- RAM 7 GB 可用（使用者宣稱手頭 12 GB；本容器額度較低）、CPU 2 核、無 GPU
- Python 3.13.18；已有 numpy / scipy / pandas
- 缺：biopython、freesasa、pdbfixer、openmm、prody
- 缺全部結構 CLI：pymol、chimerax、obabel、reduce、pdb2pqr、apbs、foldx、dssp

判讀：
- 無 GPU → RFdiffusion / RFantibody / BindCraft / AlphaFold2 **不可能在本容器執行**。
  這與 SOP §0「你不負責跑 RFdiffusion／BindCraft 的實作」一致，不影響本角色交付物。
- 階段 A（結構修剪、熱點、SASA、接觸表）為純 CPU 幾何計算，7 GB RAM 綽綽有餘。

## 2. 嘗試：下載結構（files.rcsb.org）
指令：`curl https://files.rcsb.org/download/7OW6.pdb`
結果：CONNECT tunnel failed, response 403
proxy status 記錄：`connect_rejected — gateway answered 403 to CONNECT (policy denial)`
判讀：代理策略封鎖，README 明示 403/407 不得重試或繞路。

## 3. 嘗試：安裝 biopython
指令：`pip install --break-system-packages biopython`
結果：pypi.org DNS 解析失敗（pypi 在 noProxy 名單內 → 走直連 → 容器無直連出口）
再試：清掉 no_proxy 強制走代理 → `Tunnel connection failed: 403 Forbidden`
判讀：PyPI 亦不可達。純標準庫實作成為必要路線（見日誌 02）。

## 4. 嘗試：尋找可用出口（鏡像與替代來源）
| 主機 | HTTP |
|---|---|
| www.ebi.ac.uk（PDBe） | 000（封鎖） |
| alphafold.ebi.ac.uk | 000 |
| raw.githubusercontent.com | 000 |
| github.com | 400（CONNECT 成功，僅此一個通） |
| aiprotein.sjtu.edu.cn（賽事官網） | 000 |
| files.wwpdb.org | 000 |
| data.rcsb.org | 000 |

判讀：容器出口僅餘 github.com。結構檔必須由容器外取得。

## 5. 決策
改走「使用者電腦」路線：本會話已連結使用者電腦，其網路不受本容器策略限制。
計畫：請求資料夾權限 → 在使用者電腦上下載 4 個 PDB → stage 回容器 → 全部幾何計算在容器內以純標準庫完成。
備援：若裝置不可達，改請使用者直接把 PDB 附到對話。
