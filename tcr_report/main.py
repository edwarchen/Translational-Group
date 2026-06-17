#!/usr/bin/env python3
"""TCR 免疫组库自动化报告生成系统 — CLI 入口。

Usage::

    # 生成完整报告
    python main.py run full --data-dir tests/test_data -o output/report.html

    # 只生成单样本分析
    python main.py run single --data-dir tests/test_data -o output/single_report.html

    # 只生成多样本分析
    python main.py run multi --data-dir tests/test_data -o output/multi_report.html

    # 使用配置文件
    python main.py run full --config config/custom.yaml

    # 生成模拟数据
    python main.py mock --output-dir tests/test_data

    # 列出可用图表
    python main.py charts
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure src/ is importable
sys.path.insert(0, str(Path(__file__).parent))


def cmd_run(args):
    """Generate the report."""
    from src.report.assembler import ReportAssembler

    # Load config if provided
    config = {}
    if args.config:
        import yaml
        with open(args.config) as f:
            config = yaml.safe_load(f)

    data_dir = args.data_dir or config.get("data", {}).get("qc_table", "tests/test_data")
    # Handle case where data_dir is a dict
    if isinstance(data_dir, dict):
        data_dir = args.data_dir or "tests/test_data"
    if isinstance(data_dir, str) and not Path(data_dir).is_dir():
        # Try the config's base data dir
        data_dir = args.data_dir or args.data_dir or "tests/test_data"

    # Actually, just use args.data_dir directly
    data_dir = args.data_dir or "tests/test_data"

    report_title = config.get("report", {}).get("title", "TCR 免疫组库测序报告")
    company = config.get("report", {}).get("company", "深圳海普洛斯医学检验实验室")

    assembler = ReportAssembler(
        data_dir=data_dir,
        template_dir="templates",
        report_title=report_title,
        company=company,
        static_images_dir=args.static_dir,
    )

    # single: only single-sample analysis
    # multi / full: single-sample analysis + multi-sample comparison
    skip_single = False
    skip_multi = args.mode == "single"

    assembler.assemble(
        output_path=args.output,
        sample_id=args.sample,
        skip_single=skip_single,
        skip_multi=skip_multi,
    )


def cmd_mock(args):
    """Generate mock data for testing."""
    from generate_mock_data import main as mock_main
    mock_main(args.output_dir)


def cmd_charts(args):
    """List available chart functions."""
    charts = [
        # Single sample
        ("单样本", "plot_cdr3_length_distribution", "CDR3氨基酸长度分布堆叠柱状图"),
        ("单样本", "plot_v_gene_lollipop", "V基因频率棒棒糖图"),
        ("单样本", "plot_j_gene_lollipop", "J基因频率棒棒糖图"),
        ("单样本", "plot_vj_heatmap", "V-J配对热图"),
        ("单样本", "plot_vj_sankey", "V-J配对桑基图"),
        ("单样本", "plot_clone_sunburst", "TCR克隆旭日图"),
        ("单样本", "plot_antigen_bubble", "功能注释抗原气泡图"),
        ("单样本", "plot_cancer_risk_violin", "癌症风险评估小提琴图"),
        ("单样本", "plot_clone_composition_stacked", "克隆组成分布堆叠条图"),
        ("单样本", "plot_clone_frequency_lines", "克隆频率分布线图"),
        ("单样本", "plot_aa_composition", "氨基酸组成分布堆叠条图"),
        # Multi sample
        ("多样本", "plot_pca_scatter", "PCA聚类散点图"),
        ("多样本", "plot_gene_clustermap", "V/J基因聚类热图"),
        ("多样本", "plot_clone_similarity_heatmap", "克隆相似性热图 (Morisita)"),
        ("多样本", "plot_group_comparison_box", "组间多样性指数比较箱线图"),
        ("多样本", "plot_venn_diagram", "克隆重叠韦恩图"),
        ("多样本", "plot_significant_expansion", "显著扩增TCR克隆折线图"),
    ]

    print(f"{'分类':<8} {'函数名':<40} {'说明'}")
    print("-" * 90)
    for cat, func, desc in charts:
        print(f"{cat:<8} {func:<40} {desc}")


def main():
    parser = argparse.ArgumentParser(
        description="TCR 免疫组库自动化报告生成系统",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python main.py run full -d tests/test_data -o output/report.html
  python main.py mock -o tests/test_data
  python main.py charts
        """,
    )

    sub = parser.add_subparsers(dest="command", help="Available commands")

    # "run" command
    p_run = sub.add_parser("run", help="Generate report")
    p_run.add_argument("mode", choices=["full", "single", "multi"],
                       help="Report mode: full (all), single (single-sample only), multi (multi-sample only)")
    p_run.add_argument("--data-dir", "-d", default="tests/test_data",
                       help="Path to CSV data directory")
    p_run.add_argument("--output", "-o", default="output/report.html",
                       help="Output HTML path")
    p_run.add_argument("--config", "-c", default=None,
                       help="YAML config file path")
    p_run.add_argument("--sample", "-s", default=None,
                       help="Sample ID for single-sample analysis (default: first sample)")
    p_run.add_argument("--static-dir", default=None,
                       help="Directory for static PNG overrides (optional; when set, matching PNGs replace Plotly charts)")
    p_run.set_defaults(func=cmd_run)

    # "mock" command
    p_mock = sub.add_parser("mock", help="Generate mock data")
    p_mock.add_argument("--output-dir", "-o", default="tests/test_data",
                        help="Output directory for CSV files")
    p_mock.set_defaults(func=cmd_mock)

    # "charts" command
    p_charts = sub.add_parser("charts", help="List available chart functions")
    p_charts.set_defaults(func=cmd_charts)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
