# TCR 免疫组库自动化报告生成系统 — 技术文档

## 1. 项目概览

将 TCR 免疫组库测序结果自动生成为独立的 HTML 交互式报告。图表使用 Plotly（交互）和 Logomaker（序列 logo）。

```
云平台 TSV + clonotype 文件
        │
        ▼
  convert_from_tsv.py
        │
        ▼
    13 张 CSV ──→ TCRDataLoader ──→ Chart Functions ──→ ReportAssembler ──→ Jinja2 ──→ report.html
```

**技术栈**：Python 3 + Plotly + Pandas + Jinja2 + Logomaker + Matplotlib

---

## 2. 项目目录结构

```
tcr_report/
├── main.py                      # CLI 入口
├── convert_from_tsv.py          # 数据入口：云平台 TSV + clonotype → 13 CSVs
├── convert_clonotype.py         # 备用：纯 clonotype 文件 → 13 CSVs
├── config/
│   └── default_config.yaml      # YAML 配置模板
├── src/
│   ├── data/
│   │   └── loader.py            # 13 张 CSV Schema 校验 + 加载
│   ├── charts/
│   │   ├── utils.py             # 莫兰迪配色、Plotly 主题、工具函数
│   │   ├── single_sample.py     # 单样本图表函数（17 个，10 个在报告中使用）
│   │   └── multi_sample.py      # 多样本图表函数（7 个）
│   └── report/
│       ├── assembler.py         # 报告组装器（核心编排）
│       └── renderer.py          # Jinja2 渲染引擎
├── templates/
│   └── base.html                # 单文件 HTML 模板（内联 CSS + Jinja2）
├── pngs/                        # 可选静态 PNG 替代图
├── notebooks/                   # 原始 Jupyter Notebook 参考
└── output/                      # 生成的 HTML 报告
```

---

## 3. 数据层

### 3.1 数据入口

日常分析使用 `convert_from_tsv.py`，输入云平台产出的两个文件：

```
python convert_from_tsv.py Results_vdj_tcrmatch.tsv sample.clonotypes.TRB.txt -o data_dir
```

| 输入文件 | 内容 | 提供的数据 |
|----------|------|-----------|
| `Results_vdj_tcrmatch.tsv` | 云平台汇总表（每行一个样本，326 列） | QC、多样性、V/J 基因频率、癌症风险、抗原匹配、CDR3 性质 |
| `.clonotypes.TRB.txt` | 单样本逐克隆表（每行一个克隆） | CDR3 序列、V/D/J 基因、read 计数 → 逐克隆图表 |

TSV 提供汇总级数据（图表 4.2-4.3 的大部分），clonotype 提供逐克隆数据（图表 4.1 及 V-J 配对）。两者互补，缺一不可。

**数据清洗流程**（内置在转换脚本中）：

```
原始克隆 → 去除非功能性 CDR3 (*, _, 非标准AA)
         → 去除低丰度克隆 (--min-reads, 默认 ≥2)
         → 有效 V+J 基因
         → 重归一化频率 (sum=1.0)
```

### 3.2 13 张数据表

| CSV 文件 | 核心列 | 用途 |
|----------|--------|------|
| `sample_info.csv` | SampleID, PatientID, Group, Type, SampleType | 样本信息 |
| `qc_table.csv` | Sample, rawYield, cleanYield, rawReads, cleanReads, cleanQ30, cleanRatio, clonalReads, clonalRatio | 测序质控（9 列） |
| `cdr3_length.csv` | Sample, Length, OverallFrequency, Clone_ID, CloneFrequency | CDR3 长度分布 |
| `v_gene_freq.csv` | Sample, V_gene, Frequency | V 基因频率 |
| `j_gene_freq.csv` | Sample, J_gene, Frequency | J 基因频率 |
| `vj_pairing.csv` | Sample, V_gene, J_gene, Frequency | V-J 配对 |
| `clone_diversity.csv` | Sample, CloneReads, SamplingReads, CloneCount, TopFreq, Clonality, ShannonIndex, SimpsonIndex, Evenness | 克隆多样性 |
| `clone_composition.csv` | Sample, FrequencyRange, Proportion | 克隆组成 |
| `functional_annotation.csv` | Sample, count, freq, cdr3aa, v, d, j, antigen.species | TCR 功能注释 |
| `cancer_risk.csv` | sample_id, KRAS_score, EGFR_score, REF_score, LUNG_CANCER_GDNA_score, LUNG_CANCER_TISSUE_score | 癌症风险 |
| `tcr_motif.csv` | Sample, Clone_ID, CDR3_AA, Position, AminoAcid, Frequency | CDR3 基序 |
| `antigen_summary.csv` | Sample, antigen, freq, count | 抗原匹配汇总 |
| `significant_expansion.csv` | name, sample_Ref, sample_Amp, FC, p.adj | 显著扩增 |

### 3.3 TCRDataLoader

```python
from src.data.loader import TCRDataLoader

loader = TCRDataLoader("data_dir")
data = loader.load_all()              # 一次性加载 13 张表
df_v = loader.load_v_gene_freq()      # 按需加载单表
df_ag = loader.load_antigen_summary() # 抗原汇总
```

---

## 4. 图表层

### 4.1 通用工具 (`src/charts/utils.py`)

| 对象 | 说明 |
|------|------|
| `MORANDI_COLORS` | 20 色莫兰迪调色板 |
| `OTHER_GREY` | `#E0E0E0`，用于"其他"类别 |
| `apply_tcr_theme(fig)` | Plotly 统一主题（字体 #333、标题 16px、白底、网格线） |
| `get_color_palette(n)` | 取 n 种莫兰迪色 |
| `hex_to_rgba(hex, alpha)` | 十六进制 → rgba |
| `figure_to_html_div(fig)` | Figure → `<div>` + `<script>` HTML 片段 |

图表函数签名统一：
```python
def plot_xxx(df: pd.DataFrame, *, sample_id: str = None, ...) -> go.Figure
```

### 4.2 单样本图表（报告中使用 10 张）

| 图号 | 函数 | 类型 | 数据源 |
|------|------|------|--------|
| 4.1.1 | `plot_cdr3_length_distribution` | 堆叠柱状图 | `functional_annotation.csv` |
| 4.1.2 | `plot_v_gene_cdr3_length` | 堆叠柱状图 | `functional_annotation.csv` |
| 4.1.3 | `plot_aa_composition` | 水平堆叠条图 | `functional_annotation.csv` |
| 4.1.4 | `plot_cdr3_motif` | 序列 Logo | Logomaker → `tcr_motif.csv` → base64 PNG |
| 4.1.5 | `plot_cdr3_physicochemical` | 散点图 (Scattergl) | `functional_annotation.csv` |
| 4.2.1.1 | `plot_v_gene_freq_bar` | 水平条形图 | `v_gene_freq.csv` |
| 4.2.1.2 | `plot_j_gene_freq_bar` | 水平条形图 | `j_gene_freq.csv` |
| 4.2.2.1 | `plot_vj_sankey` | 桑基图 | `vj_pairing.csv` |
| 4.3.1 | `plot_diversity_bar` | 水平条形图 | `clone_diversity.csv` |
| 4.3.2 | `plot_antigen_bar` | 棒棒糖图 | `antigen_summary.csv` |
| 4.3.3 | `plot_clone_sunburst` | 旭日图 | `functional_annotation.csv` |
| 4.3.4 | `plot_cancer_risk_violin` | 柱状图 (log10) | `cancer_risk.csv` |

### 4.3 多样本图表（10 张）

| 图号 | 函数 | 类型 | 数据源 |
|------|------|------|--------|
| 5.0.1 | `plot_clone_frequency_lines` | 折线图 | `clone_diversity.csv` |
| 5.0.2 | `plot_clone_composition_stacked` | 堆叠百分比条图 | `clone_composition.csv` |
| 5.0.3 | `plot_antigen_bubble` | 气泡图 | `functional_annotation.csv` |
| 5.1 | `plot_pca_scatter` | PCA 散点图 | `v_gene_freq.csv` |
| 5.2.1 | `plot_gene_clustermap` (V) | 聚类热图 | `v_gene_freq.csv` |
| 5.2.2 | `plot_gene_clustermap` (J) | 聚类热图 | `j_gene_freq.csv` |
| 5.3 | `plot_clone_similarity_heatmap` | Morisita 热图 | `v_gene_freq.csv` |
| 5.4 | `plot_group_comparison_box` | 箱线图 | `clone_diversity.csv` |
| 5.5 | `plot_venn_diagram` | 韦恩图 | `functional_annotation.csv` |
| 5.6 | `plot_significant_expansion` | 散点+折线图 | `significant_expansion.csv` |

---

## 5. 报告组装器

### 5.1 ReportAssembler

```python
assembler = ReportAssembler(
    data_dir="data_dir",
    template_dir="templates",
    report_title="TCR 免疫组库测序报告",
    company="深圳海普洛斯医学检验实验室",
    static_images_dir=None,  # 默认不启用静态 PNG
)
assembler.assemble(output_path="output/report.html", sample_id="1B")
```

### 5.2 执行流程

```
assemble()
  ├── loader.load_all()                        # 加载 13 张 CSV
  ├── _assemble_single_sample_charts()          # 遍历 SINGLE_CHART_SPECS
  │     ├── 检查 static_images_dir 是否有匹配 PNG
  │     │     ├── 有 → base64 <img>
  │     │     └── 无 → _make_single_chart() → Plotly → HTML div
  │     └── 返回 {"sec41": [...], "sec42": [...], "sec43": [...]}
  ├── _assemble_multi_sample_charts()           # 遍历 MULTI_CHART_SPECS
  ├── _df_to_html_table()                       # 样本信息表、质控表
  └── renderer.render(...) → renderer.write()   # Jinja2 渲染 → .html
```

### 5.3 图表分区

图表按标题前缀自动分入三个区块，模板中对应渲染：

```html
<h3>4.1 CDR3 分析</h3>
{% for card in single_sample_charts.sec41 %} ... {% endfor %}

<h3>4.2 基因分析</h3>
{% for card in single_sample_charts.sec42 %} ... {% endfor %}

<h3>4.3 克隆分析</h3>
{% for card in single_sample_charts.sec43 %} ... {% endfor %}
```

### 5.4 静态 PNG 替代（默认关闭）

通过 `--static-dir` 显式开启：

```bash
python main.py run single -d data_dir --static-dir pngs -o report.html
```

### 5.5 容错

每个图表包裹在 `try/except` 中，失败时输出 `[SKIP]` 并插入占位符，不阻断其他图表。

---

## 6. HTML 模板

`templates/base.html` — 单文件，内联 CSS + Jinja2。

**结构**：

```
封面（logo + 标题 + 单位 + 日期，100vh 全屏）
sticky header（滚动后固定在顶部）
一、背景介绍（含 1.3 SOP 流程图）
二、样本信息（{{ sample_info_table }}）
三、数据质控（{{ qc_table }}）
四、单样本分析
    4.1 CDR3 分析 → sec41 图表
    4.2 基因分析  → sec42 图表
    4.3 克隆分析  → sec43 图表
五、多样本分析
    {% if multi_sample_charts %} ...详细分析... {% else %} ...概述说明... {% endif %}
六、软件列表
七、参考文献
页脚（联系方式 + 版权）
```

模板变量：`report_title`, `company`, `report_date`, `logo_base64`, `sample_info_table`, `qc_table`, `single_sample_charts` (dict), `multi_sample_charts` (list), `version`.

---

## 7. CLI

```bash
# 转换数据
python convert_from_tsv.py Results.tsv sample.clonotypes.TRB.txt -o data_dir

# 生成单样本报告
python main.py run single -d data_dir -o output/report.html

# 生成多样本报告（含单样本分析）
python main.py run multi -d data_dir -o output/report.html -s SAMPLE_ID

# 启用静态 PNG 替代
python main.py run single -d data_dir --static-dir pngs -o report.html
```

| 命令 | 单样本图 | 多样本图 | 用途 |
|------|:---:|:---:|------|
| `run single` | 10 | 0 | 单个样本深度分析 |
| `run multi` / `run full` | 10 | 10 | 选中样本 + 多样本比较 |

---

## 8. 新增图表步骤

1. 在 `single_sample.py` 或 `multi_sample.py` 中编写函数：
```python
def plot_my_chart(df: pd.DataFrame, *, sample_id: str = None) -> go.Figure:
    fig = go.Figure(...)
    fig = apply_tcr_theme(fig)
    return fig
```

2. 在 `assembler.py` 的 `SINGLE_CHART_SPECS` 或 `MULTI_CHART_SPECS` 添加：
```python
(None, "my_chart", {"title": "图X.X.X ...", "caption": "..."}),
```

3. 在 `_make_single_chart()` 或 `_make_multi_chart()` 添加分发分支。

4. 如需新数据表，在 `loader.py` 添加 Schema + load 方法，在转换脚本中生成 CSV。

---

## 9. 依赖

```
plotly>=5.0
pandas>=1.3
numpy>=1.20
jinja2>=3.0
pyyaml>=5.4
matplotlib>=3.4
matplotlib-venn>=0.11
scipy>=1.7
logomaker>=0.8
```

Plotly.js 通过 CDN 加载：`<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>`

---

## 10. 快速开始

```bash
pip install -r requirements.txt
python convert_from_tsv.py Results_vdj_tcrmatch.tsv sample.clonotypes.TRB.txt -o data_dir
python main.py run single -d data_dir -o output/report.html
```
