# 长江杯 2026 — Target 7 (`IL-7Rα` / `CD127`) 独立分支双开操作指南 (`target7-il7ra`)

[![Open In Colab (Target 7 Isolated Worker)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/David-Liv/changjiang-protein-2026/blob/target7-il7ra/notebooks/colab_target7_IL7Ra_worker.ipynb)

> **核心结论（为什么采用独立分支 + 独立队列双开）**：
> `main` 分支现有的 `pipeline/worker.py` 会无差别扫描共享网盘的 `jobs/running/*.json` 与 `jobs/pending/*.json`，并并发追加写入同一个 `registry/design_registry.csv`。如果直接在第二个 Colab 上跑 `main` 分支，会与当前正在跑的 **Target 5 (`TL1A`)** 互相抢单或引发 Google Drive FUSE 文件覆盖冲突。
> 因此，本分支 **`target7-il7ra`** 将 **Target 7 (`IL-7Rα` Minibinder, 40–100 aa)** 完整落地在独立的模块与网盘队列命名空间中——**你可以直接打开一个新的 Colab 同时开跑 T7，对现有的 T5 Colab 100% 零干扰，且仍能唯读复用网盘已打包的 BindCraft 45 秒极速镜像缓存！**

---

## 1. 双 Colab 并发隔离对照表（T5 vs T7）

| 资源维度 | 现有 Colab #1（Target 5 `TL1A`） | 新开 Colab #2（Target 7 `IL-7Rα`） | 是否会冲突？ |
| :--- | :--- | :--- | :---: |
| **GitHub 分支** | `main` | **`target7-il7ra`** | **0 冲突** |
| **Colab Notebook** | `notebooks/colab_80gb_worker.ipynb` | **`notebooks/colab_target7_IL7Ra_worker.ipynb`** | **0 冲突** |
| **Colab 本地代码路径** | `/content/changjiang-protein-2026` | **`/content/changjiang-protein-2026-t7`** | **0 冲突** |
| **网盘任务队列目录** | `<共享网盘>/jobs/{pending,running,...}` | **`<共享网盘>/jobs_target7_IL7Ra/{pending,running,...}`** | **0 冲突（代码级硬拦截）** |
| **网盘序列注册表** | `<共享网盘>/registry/design_registry.csv` | **`<共享网盘>/registry/design_registry_T7_IL7Ra.csv`** | **0 冲突** |
| **网盘运行输出目录** | `<共享网盘>/runs/job_t5_tl1a_*` | **`<共享网盘>/runs/job_t7_il7ra_*`** | **0 冲突** |
| **ESMFold 权重锁 (`HF_HOME`)** | `<共享网盘>/cache/huggingface` | **本地 NVMe `/content/hf_cache`（从网盘唯读预载）** | **0 锁竞争** |
| **BindCraft 镜像缓存** | `<共享网盘>/cache/bindcraft_colab_cache.tar.zst` | **唯读解压同一份 `.tar.zst`（45 秒秒开，绝不覆写）** | **安全唯读复用** |

---

## 2. 新开 Colab 跑 Target 7 四步极速流程

1. **一键打开 T7 专属 Notebook**：
   点击上方徽章或直接访问：
   https://colab.research.google.com/github/David-Liv/changjiang-protein-2026/blob/target7-il7ra/notebooks/colab_target7_IL7Ra_worker.ipynb
   在 Colab 顶部菜单选择 **`执行阶段 -> 变更执行阶段类型 -> A100 GPU`**。
2. **运行『第 1 格』**：
   挂载团队共享 Google Drive（ID: `1NPtj9CpPoq1JJX42EZbdp7wfdQ8THsQg`），自动拉取 `target7-il7ra` 分支，运行 9 项数学与结构自检，并打印 `jobs_target7_IL7Ra/` 独立队列状态。
3. **运行『第 2 格』**：
   自动检测共享网盘已有的 `cache/bindcraft_colab_cache.tar.zst`，**唯读多核解压至当前 Colab 本地 SSD（仅需 ~45 秒，无需重新编译安装 12 分钟）**。
4. **运行『第 3 格』 + 『第 4 格』开跑**：
   - 在 **第 3 格** 选择批次模板（首次直接选默认的 `target7_IL7Ra_core_3DI3_HS1`），投递至 `jobs_target7_IL7Ra/pending/`；
   - 运行 **第 4 格** 启动 `pipeline/worker_target7.py`（由于 `IL-7Rα` 修剪后核心域仅 `193 aa`，复合物总长 `< 280 tokens`，单条轨迹比三聚体靶点快约 **2 倍**！）。
   - 跑完后自动在 `<共享网盘>/runs/<job_id>/af3_inputs/` 生成 **`1×Minibinder + 1×219aa 官方全长 IL-7Rα`** 的 AlphaFold 3 Server 批量验证 JSON，并在 `submission_preview/` 生成符合官方提交格式（`name,sequence`）的 Top-20 预览表（也可随时运行 **第 5 格** 进行 1 分钟 `SolubleMPNN` 增强采收并触发浏览器自动下载 AF3 JSON）。

---

## 3. Target 7 (`IL-7Rα`) 四套预设批次与多模板配额建议（SOP §D-2）

详细结构生物学审计见 [`targets/target7_IL7Ra/README.md`](targets/target7_IL7Ra/README.md)。

| 第 3 格下拉选项 | 对应模板 JSON | 靶点结构 (`Chain B` 成熟编号) | 热点斑块 | 长度范围 | 设计目的与推荐配额 |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **`target7_IL7Ra_core_3DI3_HS1`** ★ | `job_target7_IL7Ra_minibinder.json` | `t1_3DI3_coreECD.pdb` (`193 aa`, `T46I` 校正) | `B77,B79,B81,B138,B139,B192` | `55–85 aa` | **主力批次 (40%)**：官方 `3DI3` 刚性肘部核心，距全部晶体糖链 `12.0–22.6 Å` 净空 |
| `target7_IL7Ra_holo_7OPB_HS3_mixed` | `job_target7_IL7Ra_holo7OPB.json` | `t2_7OPB_holo_coreECD.pdb` (`196 aa`) | `B80,B82,B138,B139,B192,B193` | `52–78 aa` | **跨晶体混合拓扑 (25%)**：`7OPB` 实测结合构象 + `random_helicity` 采样混合 $\alpha/\beta$ |
| `target7_IL7Ra_apo_3UP1_HS4` | `job_target7_IL7Ra_apo3UP1.json` | `t3_3UP1_apo_coreECD.pdb` (`196 aa`) | `B59,B77,B138,B139,B192` | `55–85 aa` | **静息态交叉验证 (20%)**：无配体 `3UP1` 构象，确保结合天然静息态细胞表面受体 |
| `target7_IL7Ra_backMaskK_HS1` | `job_target7_IL7Ra_backMaskK.json` | `t5_3DI3_coreECD_backMaskK.pdb` (`193 aa`) | `B77,B79,B81,B138,B139,B192` | `50–75 aa` | **背侧赖氨酸封锁紧凑版 (15%)**：`N45K/N131K/N162K` 正电荷屏障防漂移 |

---

## 4. 内建赛题合规与生物物理硬闸一览

1. **SOP §B-1 输入完整性硬闸**：每次跑批前自动调用 `scripts/validate_input.py --strict-sidechains` 校验 SHA-256、链 `B` 热点存在性与零残缺重原子。
2. **SOP §C-4 逆折叠硬闸**：全程强制 `SolubleMPNN`（`mpnn_weights="soluble"`, 默认禁用易错配/易氧化的 `Cys/Met`）。
3. **SOP §C-1 新颖性硬闸（`<= 70.0%`）**：同时采用全局 Needleman–Wunsch 比对与局部滑窗扫描比对 `targets/target7_IL7Ra/targets/known_binders_IL7Ra_SEQRES.fasta`（含已发表 `7OPB_D` 55 aa Minibinder、人源 `IL-7` 152 aa、`TSLP` 126 aa），凡覆盖度 $\ge 80\%$ 且序列同一性 $> 70\%$ 的序列一律自动拦截（`REJECT:NOVELTY_FAIL_*`）。
4. **`3DI3` 晶体糖链零碰撞与 `IL-7` 竞争位点硬闸**：将每个候选复合物刚体叠合回 `3DI3`（`2.15 Å`），自动拦截与 `Asn29/Asn45/Asn131` 晶体糖链（`Chains C/D/E`）发生 `< 2.5 Å` 重原子碰撞的序列，同时确认与天然配体 `IL-7`（`Chain A`）存在 `< 3.5 Å` 直接正交空间竞争重叠。
5. **官方 219 aa 全长正交筛选**：`ESMFold + ipSAE` 打分与 `AlphaFold 3 Server` JSON 导出统一采用官方完整 **`Glu21–Asp239`（219 aa）** 序列，与赛题初筛口径 100% 对齐。
