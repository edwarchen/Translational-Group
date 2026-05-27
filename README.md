# TCR Pipeline Optimization

**MiXCR 预设参数对比与标准化分析流程** — 评估 generic-amplicon 和 rna-seq 两个 MiXCR 预设对 gDNA TCR 扩增子数据的影响，建立统一的标准化分析管线。

## 核心结论

同为 MiXCR 4.6.0 时，**generic-amplicon 和 rna-seq 的免疫组库分析结果在定量上不可区分**：

| 证据线 | 指标 | 数值 |
|---|---|---|
| CDR3 AA 序列重叠 | Jaccard | **97.5–98.8%** |
| 克隆频率相关性 | Pearson r | **> 0.999** |
| 多样性 (Shannon) | Pielou J' 差异 | **< 0.0001** |
| Morisita-Horn 交叉验证 | VDJTools vs Custom | **差异 4.2×10⁻⁶** |
| V-J 配对频率 | Raw Pearson r | **> 0.9999** |
| 聚类分析 | 同一样本 amp/rna 1−MH | **~10⁻⁵**（样本间 0.65–0.96） |

**版本升级（3.0.4 → 4.6.0）的影响远超预设差异。** 对于 gDNA 扩增子数据，推荐 generic-amplicon，但 rna-seq 也可接受（99%+ 一致）。

详细分析见 [`comparison/preset_comparison.md`](comparison/preset_comparison.md)。

---

## 管线架构

三层调用结构，单命令切换预设：

```
tcr_process_v3.sh [preset]
  │  PRESET="${1:-generic-amplicon}"    # 默认 generic-amplicon
  │
  └─ tcr_get_shell_fixed_primers_v3.py  # 生成批量任务脚本
       │  读取样本信息表 → 为每个样本生成命令行
       │
       └─ TCR_analysis_pipeline.v7.sh   # 单样本全流程
            │  R1 R2 outdir key threads chain PRESET
            │
            ├─ [1/9] fastp         → 质量修剪 + 接头去除
            ├─ [2/9] FLASH         → 双端合并
            ├─ [3/9] seqkit + V引物 → V 引物筛选
            ├─ [4/9] seqkit + J引物 → J 引物筛选
            ├─ [5/9] MiXCR align   → V/D/J 基因比对 (--preset $PRESET)
            ├─ [6/9] MiXCR assemble→ 克隆型聚类
            ├─ [7/9] export + 下游 → CDR3 统计 + 格式转换
            ├─ [8/9] VDJTools      → 多样性/片段使用/光谱分型
            └─ [9/9] total_stat    → 总体统计
```

**两个 preset 的唯一差异在 [5/9] align 步骤**：

| Preset | 参数 |
|---|---|
| `generic-amplicon` | `--preset generic-amplicon --species hs --dna --floating-left-alignment-boundary --floating-right-alignment-boundary J` |
| `rna-seq` | `--preset rna-seq --species hs` |

其余 8 个步骤完全相同。

---

## 快速开始

### 1. 运行分析管线

```bash
# 默认 generic-amplicon（gDNA 扩增子推荐）
bash tcr_process_v3.sh

# rna-seq 预设
bash tcr_process_v3.sh rna-seq

# 显式指定
bash tcr_process_v3.sh generic-amplicon
```

**前置条件**：
- 样本信息表 `fq_matched.tsv`：包含 `Sample_ID`, `chain`, `Raw_Path_R1`, `Raw_Path_R2` 列
- 服务器环境：Java 11, MiXCR 4.6.0, VDJTools 1.2.1, fastp, FLASH, seqkit
- 配置文件中的路径需根据实际环境修改

### 2. 生成预设对比报告

```bash
cd comparison/

# P1–P10: 主对比分析
python compare_presets.py

# P3b/P6b: 过滤后补充分析
python compare_mixcr_filtered.py

# P11: 多样性指数
python diversity_metrics.py

# P12: VDJTools 交叉验证
python compare_vdjtools_vs_custom.py
python compare_vdjtools_full.py

# P13: V-J 配对比对
python plot_vjusage_heatmaps.py
```

### 3. VDJTools 批量分析（独立运行）

```bash
# 默认过滤 singletons (readCount<2)
bash batch_vdjtools_presets.sh

# 不过滤
MIN_READCOUNT=1 bash batch_vdjtools_presets.sh
```

**VDJTools 分析模块**：
- Phase 1: CalcDiversityStats / CalcSegmentUsage / CalcSpectratype（每个样本）
- Phase 2: OverlapPair（同一样本 amp vs rna）
- Phase 3: RarefactionPlot（全部 10 个样本）
- Phase 4: CalcPairwiseDistances + ClusterSamples
- Phase 5: PlotFancyVJUsage（V-J 配对矩阵）

---

## 预设对比分析维度（P1–P13）

| 分析 | 内容 | 方法 |
|---|---|---|
| P1 | Alignment 效率 | 比对率、未命中率、J 缺失率 |
| P2 | Assembly 效率 | 克隆数、读段使用率、PCR 纠错 |
| P3 | 克隆检测灵敏度 | 克隆数、Unique AA（过滤前后） |
| P4 | CDR3 长度分布 | NT/AA 长度密度图 |
| P5 | V/J 基因使用 | Spearman 秩相关 |
| P6 | CDR3 AA 重叠 | Jaccard 指数、Venn 图 |
| P7 | 频率相关性 | AA 聚合 + NT 合并双重验证 |
| P9 | Top-N 一致性 | Top 10/50/100 重叠率 |
| P10 | 生物学验证 | CDR3 长度、C-end F/W、Cys 起始 |
| P11 | 多样性与一致性 | Shannon/Pielou/Morisita-Horn/Bhattacharyya |
| P12 | VDJTools 交叉验证 | 第三方工具正交验证 |
| P13 | V-J 配对一致性 | 二维 V×J 频率矩阵对比 |

**关键过滤条件**：CDR3 AA 长度 5–45 aa，readCount > 1（去 singleton）。

---

## 目录结构

```
TCR_Pipeline_Optim/
├── scripts/
│   ├── TCR_analysis_pipeline.v7.sh          # 单样本全流程（支持 preset 参数）
│   ├── tcr_process_v3.sh                    # 批量入口脚本
│   ├── tcr_get_shell_fixed_primers_v3.py    # Python 编排器（生成批量任务）
│   └── batch_vdjtools_presets.sh            # VDJTools 批量对比分析
├── comparison/
│   ├── preset_comparison.md                 # 完整分析报告
│   ├── compare_presets.py                   # P1–P10 主对比脚本
│   ├── compare_mixcr_filtered.py            # 过滤后补充分析
│   ├── diversity_metrics.py                 # P11 多样性指数计算
│   ├── compare_vdjtools_vs_custom.py        # P12.1 VDJTools vs 自定义交叉验证
│   ├── compare_vdjtools_full.py             # P12.2–12.5 VDJTools 综合分析
│   ├── plot_vjusage_heatmaps.py             # P13 V-J 配对热图
│   ├── figures_preset/                      # 全部报告图（25 张 PNG）
│   └── *.tsv                                # 中间数据表
├── vjusage/                                 # V-J 配对频率矩阵（TXT）
└── vdjtools_comparison_minCount2/           # VDJTools 分析结果
    ├── generic-amplicon/                    #   每样本 diversity/segment/spectratype
    ├── rna-seq/
    ├── overlap/                             #   OverlapPair 结果
    ├── rarefaction/                         #   RarefactionPlot 数据
    └── cluster/                             #   Pairwise distance + 聚类
```

## 依赖

| 工具 | 版本 | 用途 |
|---|---|---|
| MiXCR | 4.6.0 | TCR 比对与克隆组装 |
| VDJTools | 1.2.1 | 多样性/片段使用/重叠/聚类 |
| Java | 11+ | MiXCR / VDJTools 运行时 |
| fastp | 0.19.7 | FASTQ 质量修剪 |
| FLASH | 1.2.11 | 双端合并 |
| seqkit | — | 引物序列筛选 |
| Python | 3.8+ | 下游分析与可视化 |
| pandas / numpy / scipy | — | 数据处理与统计 |
| matplotlib | — | 绑图 |

---

## 协作指南

### 分支策略

```
main (受保护，PR + 审批强制)
├── edwarchen/dev           # 你的开发分支
├── zhangsan/feature        # 同事A的分支
├── lisi/fix                # 同事B的分支
└── ...
```

- **main** 只读保护，任何人不能直接 push
- 每个人在 `{username}/{suffix}` 分支上独立开发
- 合并到 main 必须提 PR 且至少 1 人审批

### 新成员加入

1. 仓库管理员在 **Settings → Collaborators** 邀请新成员
2. 新成员接受邀请后，在 **Actions → Create User Branch → Run workflow** 输入自己的 GitHub 用户名，自动创建 `{username}/dev` 分支
3. 克隆仓库，切换到自己的分支开始工作：

```bash
git clone https://github.com/edwarchen/TCR_Pipeline_Optim.git
cd TCR_Pipeline_Optim
git fetch origin
git checkout your-username/dev
```

### 提交代码

```bash
# 在自己的分支上开发
git add <files>
git commit -m "feat: description"
git push

# 通过 GitHub 网页创建 Pull Request → 选择 main 作为 target
# 至少 1 人审批通过后合并
```
