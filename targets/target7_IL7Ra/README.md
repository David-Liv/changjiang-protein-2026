# Target 7: IL-7Rα (CD127) — Minibinder 靶点档案（全晶体刚性肘部表位与糖链避让版）

长江杯 2026 · 靶点七 · Minibinder 赛道（40–100 aa）

## 1. 官方赛题规格与复合物晶体依据

| 项目 | 内容 |
| :--- | :--- |
| **靶点名称** | Human `IL-7Rα` / `CD127`（白细胞介素-7 受体 $\alpha$ 亚基） |
| **设计分子类型** | **Minibinder（迷你结合蛋白，单链 `40–100 aa`，20 种标准氨基酸，无标签）** |
| **UniProt ID** | `P16871`（经典异构体 `P16871-1`） |
| **官方靶点区段** | **`Glu21–Asp239`（成熟胞外域 `Mature 1–219`，共 `219 aa`）** |
| **天然寡聚状态** | **单体 І 型细胞因子受体胞外域（双纤连蛋白 III 型结构域 `FnIII D1 + D2`）** |
| **官方/晶体参考结构** | `evidence/3DI3.pdb`（2.15 Å，`IL-7:IL-7Rα` 复合物，官方指定参考结构）<br>`evidence/7OPB.pdb`（2.49 Å，Baker 实验室从头设计 Minibinder `b7` 结合人源 `IL-7Rα` 晶体结构，Cao et al., *Nature* 2022）<br>`evidence/3UP1.pdb`（2.70 Å，未结合配体的静息态 `IL-7Rα` 晶体结构）<br>`evidence/AF-P16871-F1.pdb`（AlphaFold DB v6 全长模型，提取完整官方 `21–239` 胞外域） |
| **新颖性硬闸（SOP §C-1）** | 与已知阳性结合分子比对覆盖度 $\ge 80\%$ 时，**序列同一性必须 $\le 70\%$**（排除集 `targets/known_binders_IL7Ra_SEQRES.fasta` 已收录去标签后的 `7OPB_D` 55 aa Minibinder、人源 `IL-7` 152 aa 及 `TSLP` 126 aa） |

### 官方靶点序列（`P16871` `Glu21–Asp239`，219 aa）

```text
ESGYAQNGDLEDAELDDYSFSCYSQLEVNGSQHSLTCAFEDPDVNITNLEFEICGALVEVKCLNFRKLQEIYFIETKKFLLIGKSNICVKVGEKSLTCKKIDLTTIVKPEAPFDLSVVYREGANDFVVTFNTSHLQKKYVKVLMHDVAYRQEKDENKWTHVNLSSTKLTLLQRKLQPAAMYEIKVRSIPDHYFKGFWSEWSPSYYFRTPEINNSSGEMD
```

---

## 2. 编号映射、多态性校正与结构修剪（SOP §A-1 ~ §A-4）

1. **统一编号规则（`UniProt = Mature + 20`，链 ID 统一为 `B`）**：
   - `3DI3`（Chain B）与 `3UP1` 采用去信号肽后的 **成熟编号（`Mature 1–219`）**，而 `7OPB`、`5J11` 与 `AF-P16871-F1` 采用 **UniProt 全长编号（`21–239`）**。
   - 本目录下全部 **5 个修剪靶点 PDB（`t1`–`t5`）** 均已标准化为 **`Chain B` + 成熟编号（`B1–B219`）**，确保同一套热点字符串（如 `B77,B79,B81,B138,B139,B192`）在全部 5 个模板上零偏移通用。
2. **编码区多态性校正（`Mature T46I` / `UniProt T66I`）**：
   - `3DI3` 与 `3UP1` 晶体构建体携带常见 SNP `Thr46`（`UniProt Thr66`），而官方 `P16871-1` 标准序列（以及 `7OPB`、`AF2`）在该位点为 **`Ile46`（位于 D1 背侧，距 `IL-7` 结合肘部 $>30\text{ Å}$）**。
   - `scripts/prepare_target7_pdbs.py` 已将 `t1`、`t3`、`t5` 的 `B46` 确定性标准化为无碰撞 `ILE`，确保与官方 `P16871-1` **零序列错配**。
3. **残缺侧链确定性补齐与二硫键完整性**：
   - `3DI3` 表面 14 个末端缺原子的溶剂暴露侧链已通过 `7OPB_A` 局部主链 `N-CA-C` 刚体叠合确定性补齐（`3DI3` 原有重原子位移 `0.000 Å`），5 个模板全部通过 `scripts/validate_input.py --strict-sidechains` 五项断言。
   - 全部 3 对天然二硫键（`C22–C37`、`C54–C62`、`C88–C98`）完整保留（`Sγ–Sγ = 2.02–2.06 Å`），核心修剪范围 `B17–B209`（`193 aa`）或 `B16–B211`（`196 aa`）避开 N 端无序尾段 `1–16` 与 C 端近膜柔性尾段 `210–219`（含 `WSXWS` 结构基序 `W196-S197-E198-W199-S200` 完整保留），使复合物总长度控制在 `< 285 tokens`，A100 轨迹速度较三聚体靶点提升约 **2 倍**。
   - **正交评估对齐官方全长**：进入阶段 C（`ESMFold + ipSAE`）及导出 `AlphaFold 3 Server` JSON 时，统一采用官方完整 **`219 aa`（`Glu21–Asp239`）** 靶点序列。

---

## 3. 修剪靶点文件与 SHA-256 校验表

| 文件 | 母结构 | 成熟编号范围 | 残基数 | SHA-256 (LF) | 用途与配额建议 |
| :--- | :--- | :---: | :---: | :--- | :--- |
| **`targets/t1_3DI3_coreECD.pdb`** ★ | `3DI3` (2.15 Å) | `B17–B209` | **193** | `d59f9f2848fbf4ec66e0887cb9f9830c28dcffd0bd8009d03355bece559f6ca0` | **默认主力（40% 配额）**：官方 `IL-7` 结合态核心域，零序列错配、零缺原子 |
| `targets/t2_7OPB_holo_coreECD.pdb` | `7OPB` (2.49 Å) | `B16–B211` | **196** | `34af7e4fb3a3657e33ba3b0dfb5ff406537055796c8dfd42262ae12d2fb5c924` | **跨晶体交叉模板（25% 配额）**：Minibinder 结合态实测构象，适配混合 $\alpha/\beta$ 采样 |
| `targets/t3_3UP1_apo_coreECD.pdb` | `3UP1` (2.70 Å) | `B16–B211` | **196** | `a45b173849d11b03f61147eccb7272ae6450f95a6150c70cc5d61c2dad366c54` | **静息态交叉模板（20% 配额）**：无配体 Apo 构象，确保识别天然静息态受体 |
| `targets/t4_AF2_pred_fullECD.pdb` | `AF-P16871-F1` | `B1–B219` | **219** | `f5a64944e50da22642aec56401ddff457c57bf5422cb5d0d4078dc6076b37d93` | **全长上下文对照**：完整 `Glu21–Asp239` 219 aa，用于末端位阻复核 |
| `targets/t5_3DI3_coreECD_backMaskK.pdb` | `3DI3` + 背侧封锁 | `B17–B209` | **193** | `f87f1110750bc969031c91be3cac137ab6b406b83566b00479f663dc3748516e` | **背侧糖基化位点负向封锁版（15% 配额）**：`N45K, N131K, N162K`（距热点 $>16.5\text{ Å}$）阻止模型漂移至非功能背面 |
| `targets/known_binders_IL7Ra_SEQRES.fasta` | `7OPB_D` / `IL7` / `TSLP` | — | 3 条 | `6081373fce4b7412558fa3447306c71aeef1afe2ce89df8dd92961ba42a16dbf` | **新颖性硬闸排除集**：强制拦截与 `7OPB_D`、`IL-7`、`TSLP` 同一性 $\ge 70\%$ 的序列 |

---

## 4. 跨晶体表位刚性、糖基化位阻审计与四组热点字符串

### 4.1 表位刚性与 N-糖基化避让结论

- **极高跨构象刚性**：`D1–D2` 肘部核心表位 11 个残基（`Y59, K77, F79, L80, L81, I82, K138, Y139, H191, Y192, F193`）在 `3DI3`（`IL-7` 结合态）、`7OPB`（Minibinder 结合态）、`3UP1`（Apo 静息态）与 `AF2` 四套结构间的两两 `Cα RMSD` 仅为 **`0.42–0.55 Å`**，无需诱导契合即可直接结合。
- **晶体 N-糖链精确避让**：`IL-7Rα` 胞外域含 6 个 `N-X-S/T` 糖基化序列子（成熟编号 `N29, N45, N131, N162, N212, N213`），其中 `3DI3` 晶体直接解析出 `N29`（`Chain C`）、`N45`（`Chain D`）、`N131`（`Chain E`）三处壳二糖（`NAG-NAG`）：
  - 注意 `B31 (Ser31)` 本身属于 `N29-G30-S31` 糖基化序列子（距 `Asn29` 晶体 `NAG` 仅 `7.67 Å`），故**不将 `B31` 放入默认主力热点组 `HS1`**；
  - 默认主力组 **`HS1` (`B77,B79,B81,B138,B139,B192`)** 各残基距 `3DI3` 全部晶体糖链重原子均保持 **`12.0–22.6 Å` 安全净空**，同时与 `IL-7`（`3DI3`）、`TSLP`（`5J11`）及 `7OPB_D` 正交阻断位点 100% 重叠；
  - `pipeline/worker_target7.py` 在采收阶段自动将每个候选复合物叠合回 `evidence/3DI3.pdb`，硬性拦截与 `Chains C/D/E` 晶体糖链发生重原子碰撞（`< 2.5 Å` 碰撞数 `> 2`）的构象，并验证与 `IL-7`（`Chain A`）的直接空间竞争重叠。

### 4.2 热点字符串速查（链 ID 统一为 `B`，成熟编号）

```text
HS1  主力肘部疏水/盐桥双锚核心（配 t1_3DI3_coreECD / t5_3DI3_coreECD_backMaskK，默认推荐）
  BindCraft   : B77,B79,B81,B138,B139,B192
  UniProt 对应: K97, F99, L101, K158, Y159, Y212

HS2  上沿扩展组（含 D1 顶部 B31/Ser51，仅限小尺寸 48–68 aa 探索）
  BindCraft   : B31,B77,B79,B81,B139,B192
  UniProt 对应: S51, K97, F99, L101, Y159, Y212

HS3  下肘部疏水微移组（配 t2_7OPB_holo_coreECD + random_helicity 混合拓扑）
  BindCraft   : B80,B82,B138,B139,B192,B193
  UniProt 对应: L100, I102, K158, Y159, Y212, F213

HS4  Apo 静息态上鞍部多样性组（配 t3_3UP1_apo_coreECD）
  BindCraft   : B59,B77,B138,B139,B192
  UniProt 对应: Y79, K97, K158, Y159, Y212
```

---

## 5. 跑批前与产出后质量门控命令

```bash
# 1. 跑批前输入结构、热点与全原子侧链完整性校验（SOP §B-1）
python scripts/validate_input.py \
  --target targets/target7_IL7Ra/targets/t1_3DI3_coreECD.pdb \
  --hotspots B77,B79,B81,B138,B139,B192 \
  --sha d59f9f2848fbf4ec66e0887cb9f9830c28dcffd0bd8009d03355bece559f6ca0 \
  --strict-sidechains

# 2. 独立队列干跑自检（零触碰 main 分支 T5 的 jobs/ 队列与注册表）
python pipeline/worker_target7.py \
  --drive-root ./_local_t7_dryrun \
  --job-file jobs/templates/job_target7_IL7Ra_minibinder.json \
  --dry-run

# 3. 产出序列后组成与新颖性复核（SOP §C-1 / §C-2）
python scripts/check_expressibility.py designs.fasta
python scripts/check_novelty.py targets/target7_IL7Ra/targets/known_binders_IL7Ra_SEQRES.fasta
```
