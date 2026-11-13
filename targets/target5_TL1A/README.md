# Target 5: TL1A (TNFSF15) — Minibinder 靶点档案（v2 双链凹槽裁剪加速版）

长江杯 2026 · 靶点五 · Minibinder 赛道（40–100 aa）

## 1. 官方赛题规格与复合物晶体依据

| 项目 | 内容 |
| :--- | :--- |
| **靶点名称** | TL1A / TNFSF15（肿瘤坏死因子超家族成员 15） |
| **设计分子类型** | **Minibinder（迷你结合蛋白，单链 40–100 aa）** |
| **UniProt ID** | `O95150` / `Q8NFE9`（胞外域 `Leu72–Leu251`） |
| **天然寡聚状态** | **紧密非共价同源三聚体（Homotrimer, Chains A+B+C）** |
| **参考晶体结构** | `evidence/2RE9.pdb`（2.10 Å 天然三聚体）+ `evidence/3K51.pdb`（TL1A–DcR3 受体复合物晶体结构） |

## 2. 编号映射、结构修剪与负向设计（SOP §A-3 / §A-4）

- **编号映射**：`2RE9` 中 `PDB 1–180` $\leftrightarrow$ `UniProt 72–251`（`UniProt = PDB + 71`）；`3K51` 中 `Chain A 5–184` $\leftrightarrow$ `2RE9 1–180`（`3K51 = 2RE9 + 4`）。
- **为什么停用裸单体 `t1_2RE9_chainA.pdb`**：
  首轮测试（`job_t5_tl1a_minibinder_001`）证明，单独切出 Chain A 会将原本掩埋在三聚体核心的 65 个疏水 $\beta$-折叠残基暴露于溶剂，导致生成模型 100% 忽略外侧热点并钻入天然 Chain B/C 所在的空间（与天然三聚体重度空间碰撞）。
- **`t2_2RE9_grooveAB_cropped.pdb`（双链结合槽裁剪加速 + Lys 负向封锁版，314 aa）**：
  1. **保留完整天然受体凹槽**：保留构成天然 DcR3/DR3 中和结合槽的相邻双链 `Chain A` 与 `Chain B`。
  2. **柔性末端修剪**：在脯氨酸铰链 `Pro24`（B-factor 由 `67–104 Å²` 骤降至 `32.5 Å²`）处切除 `17–23` 柔性 N 端尾段，保留零断链核心域 `24–180`（每条链 157 aa，双链共 **314 aa**，加 50–70 aa Minibinder 总长 $\le 384$ aa，A100 轨迹速度比完整三聚体快约 **2.8 倍**）。
  3. **跨链唯一编号与晶体缺原子补齐**：`Chain A` 编号为 `24–180`，`Chain B` 偏移 `+200` 编号为 `224–380`；用对称链 `A39` 精确补齐 `B239 (PHE)` 缺失苯环侧链，严格通过 `validate_input.py --strict-sidechains` 全部五项断言。
  4. **背侧三聚体切面负向设计（Lysine Block，SOP §A-4 规则 1）**：
     对背侧因移除 `Chain C` 而暴露（$\Delta\text{SASA} = +26.5 \sim +96.5\text{ Å}^2$）且远离 `A:B` 受体凹槽（$>18.9\text{ Å}$）的 8 个朝外大疏水残基突变为全原子无碰撞 **赖氨酸 `LYS (K)`**（构建 `+8` 正电荷亲水屏障，彻底驱离模型不再打向缺失的 `Chain C` 侧）：
     - `Chain A`：`M125K`、`M146K`、`F147K`、`L180K`
     - `Chain B`：`L254K`（原 B54）、`F256K`（原 B56）、`F338K`（原 B138）、`F374K`（原 B174）

## 3. 修剪文件与 SHA-256

| 文件 | 范围 | 残基数 | SHA-256 (LF) | 状态说明 |
| :--- | :--- | :---: | :--- | :--- |
| `targets/t2_2RE9_grooveAB_cropped.pdb` | `A24–180` + `B224–380` | **314** | `d9120f7a63ba10dfd52f916db2ac5a2e4c16ceb49edf579e4a7b33ae9de9a35b` | ✅ **当前主力**：双链受体凹槽 + 8×Lys 背侧封锁，零残缺侧链 |
| `targets/t1_2RE9_chainA.pdb` | `A17–180` | 164 | `69e31cbc88fe26fa407a20b3959eff14b74716b08f9988b9bfb5ade748ced9eb` | ⚠️ **仅供历史溯源**：裸单体暴露三聚体核心，禁止用于新批次生成 |

## 4. 精确中和表位热点（基于 `3K51` DcR3–TL1A 晶体复合物叠合）

- **主力中和凹槽热点斑块**：`A52,A87,A137,A138,A167,B317`
  - 左槽壁（`Chain A`）：AA' loop `A52 (LEU)`、侧面 `A87 (MET)`、`A137 (TRP)`、`A138 (PHE)`、GH loop `A167 (TYR)`
  - 右槽壁（`Chain B`）：DE loop 核心锚点 `B317 (TYR)`（即文献 `Tyr121` / `2RE9 B117`，突变后直接废除 DR3/DcR3 结合）
- **推荐生成长度**：`50,70`（总 token $\le 384$，符合官方 Minibinder 40–100 aa 要求且表达率更高）。
- **后处理双保险**：
  1. 自动叠合 `evidence/2RE9.pdb` 天然三聚体检查与 `Chain C` 零空间碰撞，并在野生型 `1:3` 同源三聚体上跑 `ESMFold + Dunbrack merged ipSAE`。
  2. 针对 BindCraft 内部单序列单体 PTM 对双链浅凹槽的系统性盲筛（SOP §D-1），使用 `scripts/harvest_groove_backbones.py`（或 Colab Notebook **第 5 格**）对 `Trajectory/Relaxed/` 中 `i_pTM >= 0.80` 且 `Chain C` 零碰撞的凹槽骨架执行双模式 **SolubleMPNN（`fix_intf` + `full`）+ ESMFold 单体 `scRMSD` 自洽门控 + 可表达性/新颖性硬闸**。
  3. 自动按跨骨架轮询优先级导出 **AlphaFold 3 Server (`1×Minibinder + 3×野生型 TL1A`)** 批量验证 JSON 文件供最终正交复核。
