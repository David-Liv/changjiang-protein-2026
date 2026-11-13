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

## 4. 精确中和表位热点与双批次互补策略（基于 `3K51` DcR3–TL1A 晶体复合物叠合）

- **主力批次（`jobs/templates/job_target5_TL1A_minibinder.json`, `002`）**：
  - **热点斑块**：`A52,A87,A137,A138,A167,B317`（左槽壁 `Chain A` AA'/GH 环 + 右槽壁 `Chain B` DE loop `B317/Tyr121`），长度 `50–70 aa`，35 条轨迹。
- **互补增强批次（`jobs/templates/job_target5_TL1A_grooveAB_003.json`, `003` v2.1）**：
  - **对称跨槽 2:2 热点斑块**：`A52,A167,B317,B325`（跨槽 `Cα` 距离收敛于 `12.7–19.9 Å`，剔除偏外侧的 `A87`，将左槽壁 `A52 (LEU)` + `A167 (TYR)` 与右槽壁核心功能热点 `B317 (TYR, Tyr121)` + 槽底刚性锚点 `B325 (MET, B=29 Å²)` 对称配对，显著提升右槽壁 `B317` 命中率）。
  - **拓扑多样性增强（SOP §D-2）**：长度放宽至 `55–78 aa`（总 token $\le 392$），并开启 `advanced_overrides: {"random_helicity": true, "omit_AAs": "CM"}`，突破单一全 $\alpha$-螺旋束局限，补充混合 $\alpha/\beta$ 跨槽拓扑。
- **v2.1 黄金凹槽骨架采收与正交门控（`scripts/harvest_groove_backbones.py` / Colab 第 5 格）**：
  1. **三维空间与竞争阻断审计**：自动叠合 `evidence/2RE9.pdb` 验证与缺失 `Chain C` 零重原子碰撞（`< 2.5 Å`），同时将 `evidence/3K51.pdb`（`TL1A–DcR3` 复合物）叠合至凹槽，定量计算候选分子与天然受体的重原子空间占位重叠数（`receptor_clashes`）及 N/C 端溶剂暴露距离（`d_nterm`, `d_cterm`，规避 BLI His-tag 位阻）。
  2. **三模式 SolubleMPNN 逆折叠（`--omit_AAs CXM`）**：
     - `core_fix`：仅固定最深处 `4–10` 个 `<= 4.0 Å` 高置信度疏水锚点（`F/Y/W(<=2)/L/I/V`），释放边缘极性残基以消除掩埋未满足氢键洞（`UnsatHbonds`）；
     - `full`：全界面 + 表面重设计，搜索全新极性氢键网络；
     - `fix_intf`（受控开启）：仅当原轨迹 `UnsatHbonds <= 4`、`Surface_Hydrophobicity <= 0.42` 且无 `WW`/`W>3` 时启用。
  3. **SOP 序列物理门控与三维姿态去冗余**：强制拦截聚丙氨酸假象（`Ala <= 16%`、禁 `AAAA`）、色氨酸/甲硫氨酸风险（禁 `WW`、`W <= 3`、`M <= 3`）、表面带电亲水度（`E+K >= 26%`、`|net_charge_pH74| <= 10.5`）及 ESMFold 单体自洽门控（`scRMSD <= 2.2 Å`、`pLDDT >= 78.0`，混合 $\alpha/\beta$ 放宽至 `74.0`）；在 `scRMSD <= 1.05 Å` 平台期内按界面物理质量与 `B317` 命中综合排序，并通过三维接触斑块聚类（`Jaccard >= 0.55` 且质心距 `<= 6.0 Å`）按簇轮询导出幂等更新的 **AlphaFold 3 Server (`1×Minibinder + 3×野生型 TL1A`)** 批量验证 JSON。
