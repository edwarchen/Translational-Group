# TCR 免疫组库自动化报告生成系统

> 从 CSV 数据自动生成交互式 HTML TCR 分析报告，部署于云平台。

## 速查

```bash
python convert_from_tsv.py Results.tsv sample.clonotypes.TRB.txt -o data_dir
python main.py run single -d data_dir -o output/report.html
python main.py run multi -d data_dir -o output/report_multi.html
```

## 架构

```
CSV → TCRDataLoader → Chart Functions (Plotly) → ReportAssembler → Jinja2 → HTML
```

| 层 | 文件 | 职责 |
|---|------|------|
| CLI | `main.py` | argparse 入口，3 种子命令 |
| 数据 | `src/data/loader.py` | 13 张 CSV schema 校验 + 加载（含 antigen_summary） |
| 转换 | `convert_from_tsv.py` | 云平台 TSV + clonotype → 12 CSVs（生产入口） |
| 转换 | `convert_clonotype.py` | 纯 clonotype 文件 → 12 CSVs（备用入口） |
| 图表 | `src/charts/single_sample.py` | 17 个 Plotly 函数（10 个用于单样本报告） |
| 图表 | `src/charts/multi_sample.py` | 7 个多样本 Plotly 函数 |
| 图表 | `src/charts/utils.py` | 莫兰迪配色、Plotly 主题、`figure_to_html_div()` |
| 组装 | `src/report/assembler.py` | 核心编排器：数据→图表→卡片→渲染；图表按标题前缀自动分入 4.1/4.2/4.3 |
| 渲染 | `src/report/renderer.py` | Jinja2 引擎，logo base64 内嵌 |
| 模板 | `templates/base.html` | 单文件 HTML 模板，sticky header |

## 图表函数约定

- 签名：`def plot_xxx(df: pd.DataFrame, *, sample_id: str = None) -> go.Figure`
- 所有图使用 `apply_tcr_theme(fig)` 统一风格
- 莫兰迪色板：`MORANDI_COLORS`（20 色），取色用 `get_color_palette(n)`
- 静态 PNG 替代机制（默认关闭）：使用 `--static-dir pngs` 显式开启，匹配 PNG 以 base64 嵌入替代 Plotly

## 新增图表的步骤

1. 在 `single_sample.py` 或 `multi_sample.py` 中编写函数
2. 在 `assembler.py` 的 `SINGLE_CHART_SPECS` 或 `MULTI_CHART_SPECS` 添加条目
3. 在 `_make_single_chart()` 或 `_make_multi_chart()` 添加分发分支
4. 如需静态替代，将 PNG 放入 `pngs/` 并在条目中指定文件名

## 报告模式

| 命令 | 单样本图 | 多样本图 | 用途 |
|------|:---:|:---:|------|
| `run single` | 10 张 | 0 | 单个样本深度分析 |
| `run multi` / `run full` | 10 张 | 10 张 | 选中样本分析 + 多样本比较 |

`-s SAMPLE_ID` 指定单样本分析的目标样本。

## 红线

- 报告必须自包含：所有图片 base64 内嵌，无外部文件引用（Plotly.js CDN 除外）
- CSV schema 在 `loader.py` 中定义，修改列名必须同步更新 schema
- `figure_to_html_div()` 返回 `<div>` + `<script>` 片段，不要用 `fig.write_html()`
- 模板中图表编号（图X.X.X）在 `assembler.py` 的 specs 中定义，不在模板硬编码
- `run multi` 和 `run full` 行为必须一致（均生成单样本 + 多样本）
- 图表容错：每个图包裹在 try/except 中，失败不阻断其他图

## 深入文档

- 完整技术文档：`DOCUMENTATION.md`
- 原始参考报告：`TCR_report_all.pdf`（20 页）
- 测试数据：`tests/test_data/`（12 张 CSV）
- 依赖：`logomaker`（图4.1.4 序列 logo）
