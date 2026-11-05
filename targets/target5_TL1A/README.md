# Target 5: TL1A (TNFSF15) — Minibinder 靶点档案

长江杯 2026 · 靶点五 · Minibinder 赛道（40–100 aa）

## 1. 官方赛题规格

| 项目 | 内容 |
| :--- | :--- |
| **靶点名称** | TL1A / TNFSF15（肿瘤坏死因子超家族成员 15） |
| **设计分子类型** | **Minibinder（迷你结合蛋白，单链 40–100 aa）** |
| **UniProt ID** | `O95150` / `Q8NFE9`（胞外域 `Leu72–Leu251`） |
| **官方参考结构** | `2RE9`（2.10 Å X-ray 晶体结构，同源三聚体 Chains A/B/C） |

## 2. 编号映射与结构核查

- `2RE9` 的 `DBREF` 为：`PDB 1–180` $\leftrightarrow$ `UniProt 72–251`（即 `UniProt = PDB + 71`）。
- `PDB 0`（GLY）为表达标签，且 `PDB 0–16` 在晶体中无电子密度（REMARK 465）。
- 有完整三维坐标的残基范围为 **PDB `17–180`**（对应 `UniProt 88–251`，共 **164 个残基**，全主链与侧链原子 100% 完整）。

## 3. 修剪文件与 SHA-256

| 文件 | 范围 (PDB) | 残基数 | SHA-256 | 说明 |
| :--- | :--- | :--- | :--- | :--- |
| `targets/t1_2RE9_chainA.pdb` | Chain A `17–180` | 164 | `69e31cbc88fe26fa407a20b3959eff14b74716b08f9988b9bfb5ade748ced9eb` | 单体结构，已通过 `validate_input.py` 五项断言，零残缺原子 |

## 4. 候选热点（Hotspots）建议

- **探索组 1（受体结合槽侧表面残基）**：`A33,A86,A88,A148,A150`（待生物组最终复核确认后也可留空 `"null"` 让 AlphaFold2 自动探索高分表面作为对照）。
- **推荐生成长度**：`55,85`（符合官方 Minibinder 40–100 aa 要求）。
