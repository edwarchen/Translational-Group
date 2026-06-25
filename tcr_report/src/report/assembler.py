"""Report assembler: orchestrates data loading, chart generation, and HTML rendering.

This is the central module that connects all components:
    CSV data → Charts (Plotly) → HTML snippets → Jinja2 template → report.html

Supports both interactive Plotly charts AND static PNG images (base64-embedded).
When a matching static image exists, it takes priority over Plotly generation.
"""

from __future__ import annotations

import base64
from pathlib import Path
from typing import Optional

import pandas as pd
import plotly.graph_objects as go

from src.data.loader import TCRDataLoader
from src.charts.single_sample import (
    plot_cdr3_length_distribution,
    plot_v_gene_cdr3_length,
    plot_j_gene_freq_bar,
    plot_v_gene_freq_bar,
    plot_vj_sankey,
    plot_clone_sunburst,
    plot_clone_frequency_lines,
    plot_clone_composition_stacked,
    plot_antigen_bubble,
    plot_cancer_risk_violin,
    plot_aa_composition,
    plot_cdr3_motif,
    plot_cdr3_physicochemical,
    plot_diversity_bar,
    plot_antigen_bar,
)
from src.charts.multi_sample import (
    plot_pca_scatter,
    plot_gene_clustermap,
    plot_clone_similarity_heatmap,
    plot_group_comparison_box,
    plot_venn_diagram,
    plot_significant_expansion,
)
from src.charts.utils import figure_to_html_div
from src.report.renderer import ReportRenderer


class ReportAssembler:
    """Assembles the complete TCR report from CSV data.

    Usage::

        assembler = ReportAssembler(
            data_dir="tests/test_data",
            static_images_dir="pngs",
        )
        assembler.assemble(output_path="output/report.html")
    """

    def __init__(
        self,
        data_dir: str | Path = "tests/test_data",
        template_dir: str | Path = "templates",
        report_title: str = "TCR 免疫组库测序报告",
        company: str = "深圳海普洛斯医学检验实验室",
        version: str = "1.0.0",
        static_images_dir: str | None = None,
    ):
        self.data_dir = Path(data_dir)
        self.loader = TCRDataLoader(data_dir)
        self.renderer = ReportRenderer(
            template_dir=template_dir,
            logo_path="海普洛斯-蓝橙.png",
        )
        self.report_title = report_title
        self.company = company
        self.version = version
        self.static_dir = Path(static_images_dir) if static_images_dir else None

    # ── Static image helper ──────────────────────────────────────────────

    def _png_to_base64_img(self, filename: str) -> str | None:
        """Convert a PNG file to a base64 data URI <img> tag.

        Returns None if the file doesn't exist.
        """
        if self.static_dir is None:
            return None
        filepath = self.static_dir / filename
        if not filepath.exists():
            return None
        b64 = base64.b64encode(filepath.read_bytes()).decode("utf-8")
        return (
            f'<img src="data:image/png;base64,{b64}" '
            f'class="static-img" alt="{filename}" style="max-width:100%;">'
        )

    def _static_or_fig(
        self,
        static_filename: str,
        fig_factory,
        title: str,
        caption: str | None = None,
        *fig_args,
        **fig_kwargs,
    ) -> dict:
        """Return a card dict — static image if available, otherwise Plotly figure.

        Parameters
        ----------
        static_filename : str
            PNG filename to check in static_images_dir.
        fig_factory : callable
            Function that returns a plotly figure (called if static not found).
        title : str
            Chart title for the card.
        caption : str or None
            Figure caption.
        fig_args, fig_kwargs
            Passed to fig_factory.
        """
        img_html = self._png_to_base64_img(static_filename)
        if img_html is not None:
            return {"title": title, "html": img_html, "caption": caption or ""}

        try:
            fig = fig_factory(*fig_args, **fig_kwargs)
            return self._fig_to_card(fig, title=title, caption=caption)
        except Exception as e:
            print(f"  [WARN] {title}: {e}")
            return {"title": title, "html": f'<p style="color:#999;">Chart unavailable: {e}</p>', "caption": ""}

    def _fig_to_card(
        self,
        fig: go.Figure,
        title: str,
        caption: str | None = None,
    ) -> dict:
        """Convert a Plotly figure to a report card dict."""
        return {
            "title": title,
            "html": figure_to_html_div(fig),
            "caption": caption or "",
        }

    def _df_to_html_table(self, df: pd.DataFrame, max_rows: int = 20) -> str:
        """Convert a DataFrame to an HTML table string."""
        if len(df) > max_rows:
            df = df.head(max_rows)
            note = f'<p style="font-size:12px;color:#888;">Showing {max_rows} of many rows…</p>'
        else:
            note = ""

        return note + df.to_html(
            classes="data-table",
            index=False,
            border=0,
            float_format="%.4f",
            na_rep="—",
        )

    # ── Chart assembly ───────────────────────────────────────────────────

    # Maps chart → (static_filename, title, caption)
    SINGLE_CHART_SPECS = [
        # (static_png_filename, chart_key, {title, caption})
        # ── 4.1 CDR3 Analysis ──
        ("optimized_cdr3_length_conda_dl.png", "cdr3_length", {
            "title": "图4.1.1 Top10克隆在不同长度CDR3中的分布",
            "caption": "横坐标表示CDR3序列长度，纵坐标表示克隆的频率。灰色为其他克隆，彩色为频率Top10的CDR3克隆，不同颜色代表不同克隆序列。",
        }),
        ("optimized_vgene_cdr3_length_conda_dl.png", "v_gene_cdr3", {
            "title": "图4.1.2 Top12克隆V基因在不同长度CDR3中的分布",
            "caption": "图展示在不同CDR3长度分布下，频率top12克隆的V基因分布情况。",
        }),
        (None, "aa_composition", {
            "title": "图4.1.3 氨基酸组成分布",
            "caption": "横坐标表示不同氨基酸的组成比例，不同颜色代表不同的氨基酸类型。统计样本中所有CDR3序列的氨基酸组成。",
        }),
        (None, "cdr3_motif", {
            "title": "图4.1.4 高丰度TCR克隆CDR3β基序",
            "caption": "T细胞受体（TCR）的互补决定区3（CDR3）是直接参与抗原识别的核心区域。对单个样本中丰度最高的前100个TCR克隆进行统计，展示CDR3氨基酸序列在各个位置上出现的残基及其频率，可直观反映序列的保守性和功能特征。",
        }),
        (None, "cdr3_physicochemical", {
            "title": "图4.1.5 CDR3 理化性质分析",
            "caption": "散点图展示样本中每条CDR3序列的平均疏水性（Kyte-Doolittle标度）与净电荷分布。每个点代表一条CDR3氨基酸序列，点大小按克隆频率加权。虚线为中性参考线，四个象限分别标识疏水/亲水与正/负电荷的组合分布。",
        }),
        # ── 4.2 Gene Analysis ──
        (None, "v_gene_freq_bar", {
            "title": "图4.2.1.1 V基因频率分布",
            "caption": "横坐标表示V基因种类，纵坐标为使用频率。V基因使用频率以水平条形图展示，按频率降序排列。",
        }),
        (None, "j_gene_freq_bar", {
            "title": "图4.2.1.2 J基因频率分布",
            "caption": "横坐标表示J基因种类，纵坐标为使用频率。J基因使用频率以水平条形图展示，按频率降序排列。",
        }),
        ("VJ_pairing.png", "vj_sankey", {
            "title": "图4.2.2.1 V-J基因配对Circos图",
            "caption": "每个颜色块代表一种V/J基因，颜色块越宽频率越高；连线表示V-J组合。",
        }),
        # ── 4.3 Clone Analysis ──
        (None, "diversity_bar", {
            "title": "图4.3.1 克隆多样性指标",
            "caption": "Shannon指数、Simpson指数、Evenness均匀度、Clonality克隆性。Shannon越高多样性越好，Clonality越低说明无优势克隆。",
        }),
        (None, "clone_sunburst", {
            "title": "图4.3.3 TCR克隆分布旭日图",
            "caption": "按圈展示不同克隆类型的频率分布及对应的氨基酸序列。",
        }),
        (None, "antigen_bar", {
            "title": "图4.3.2 功能注释 — 抗原匹配",
            "caption": "展示该样本中匹配到已知抗原特异性TCR的频率和克隆数。匹配基于IEDB/VDJdb等公开数据库。",
        }),
        (None, "cancer_risk", {
            "title": "图4.3.4 癌症风险评估",
            "caption": "展示不同样本组的打分结果，包含KRAS、EGFR、健康对照、肺癌血液和肺癌组织的TCR信号评分。",
        }),
    ]

    MULTI_CHART_SPECS = [
        # ── Moved from single-sample (require multiple samples) ──
        (None, "clone_frequency", {
            "title": "图5.0.1 克隆频率分布",
            "caption": "纵坐标表示样本的克隆累积数量，横坐标表示克隆的频率，不同颜色的曲线表示不同样品。",
        }),
        (None, "clone_composition", {
            "title": "图5.0.2 克隆组成分布",
            "caption": "横坐标代表样本，纵坐标代表不同频率范围组的克隆种类数的占比。",
        }),
        (None, "antigen_bubble", {
            "title": "图5.0.3 功能注释气泡图",
            "caption": "气泡图展示每个样本中的TCR注释结果。横坐标代表不同样本，纵坐标代表注释到的抗原物种。",
        }),
        # ── 5.1 PCA ──
        (None, "pca", {
            "title": "图5.1 样本间PCA聚类分析",
            "caption": "每个点代表一个样本。彼此间越靠近的样品，相似性越高。",
        }),
        (None, "v_clustermap", {
            "title": "图5.2.1 V基因聚类热图",
            "caption": "颜色表示对应样本V基因使用频率，红为高频，蓝为低频。上方树形图表示样本聚类，左侧表示V基因聚类。",
        }),
        (None, "j_clustermap", {
            "title": "图5.2.2 J基因聚类热图",
            "caption": "颜色表示对应样本J基因使用频率，红为高频，蓝为低频。",
        }),
        (None, "clone_similarity", {
            "title": "图5.3 组间样本克隆相似性热图",
            "caption": "颜色代表Morisita指数大小，红色表示相似性高，蓝色相似性低。",
        }),
        (None, "group_box", {
            "title": "图5.4 组间多样性指数比较",
            "caption": "箱线图展示分组间ShannonIndex、SimpsonIndex、Evenness比较。",
        }),
        (None, "venn", {
            "title": "图5.5 样本间TCR克隆重叠关系（韦恩图）",
            "caption": "每个圈代表一例样本的TCR克隆数量，交叉部分代表共享克隆。",
        }),
        (None, "expansion", {
            "title": "图5.6 发生显著扩增的TCR克隆",
            "caption": "展示不同条件下鉴定到的T细胞克隆频率变化。FC > 1, p.adj < 0.05。",
        }),
    ]

    def _make_single_chart(self, chart_key: str, sample_id: str) -> go.Figure | None:
        """Generate a single-sample chart by key. Returns None if it fails."""
        try:
            if chart_key == "cdr3_length":
                df_fa = self.loader.load_functional_annotation()
                return plot_cdr3_length_distribution(df_fa, sample_id=sample_id)
            elif chart_key == "v_gene_cdr3":
                df_fa = self.loader.load_functional_annotation()
                return plot_v_gene_cdr3_length(df_fa, sample_id=sample_id)
            elif chart_key == "aa_composition":
                df_fa = self.loader.load_functional_annotation()
                return plot_aa_composition(df_fa, sample_id=sample_id)
            elif chart_key == "cdr3_motif":
                df_m = self.loader.load_tcr_motif()
                return plot_cdr3_motif(df_m, sample_id=sample_id)
            elif chart_key == "cdr3_physicochemical":
                df_fa = self.loader.load_functional_annotation()
                return plot_cdr3_physicochemical(df_fa, sample_id=sample_id)
            elif chart_key == "j_gene_freq_bar":
                df_j = self.loader.load_j_gene_freq()
                return plot_j_gene_freq_bar(df_j, sample_id=sample_id)
            elif chart_key == "v_gene_freq_bar":
                df_v = self.loader.load_v_gene_freq()
                return plot_v_gene_freq_bar(df_v, sample_id=sample_id)
            elif chart_key == "vj_sankey":
                df_vj = self.loader.load_vj_pairing()
                return plot_vj_sankey(df_vj, sample_id=sample_id)
            elif chart_key == "diversity_bar":
                df_cd = self.loader.load_clone_diversity()
                return plot_diversity_bar(df_cd, sample_id=sample_id)
            elif chart_key == "clone_sunburst":
                df_fa = self.loader.load_functional_annotation()
                return plot_clone_sunburst(df_fa, sample_id=sample_id)
            elif chart_key == "antigen_bar":
                df_ag = self.loader.load_antigen_summary()
                return plot_antigen_bar(df_ag, sample_id=sample_id)
            elif chart_key == "cancer_risk":
                df_cr = self.loader.load_cancer_risk()
                return plot_cancer_risk_violin(df_cr)
        except Exception as e:
            print(f"  [SKIP] {chart_key}: {e}")
        return None

    def _make_multi_chart(self, chart_key: str) -> go.Figure | None:
        """Generate a multi-sample chart by key. Returns None if it fails."""
        try:
            if chart_key == "clone_frequency":
                df_cd = self.loader.load_clone_diversity()
                return plot_clone_frequency_lines(df_cd)
            elif chart_key == "clone_composition":
                df_cc = self.loader.load_clone_composition()
                return plot_clone_composition_stacked(df_cc)
            elif chart_key == "antigen_bubble":
                df_fa = self.loader.load_functional_annotation()
                return plot_antigen_bubble(df_fa)
            elif chart_key == "pca":
                df_v = self.loader.load_v_gene_freq()
                df_si = self.loader.load_sample_info()
                return plot_pca_scatter(df_v, group_df=df_si)
            elif chart_key == "v_clustermap":
                df_v = self.loader.load_v_gene_freq()
                return plot_gene_clustermap(df_v, gene_type="V")
            elif chart_key == "j_clustermap":
                df_j = self.loader.load_j_gene_freq()
                return plot_gene_clustermap(df_j, gene_type="J", gene_col="J_gene")
            elif chart_key == "clone_similarity":
                df_v = self.loader.load_v_gene_freq()
                return plot_clone_similarity_heatmap(df_v)
            elif chart_key == "group_box":
                df_cd = self.loader.load_clone_diversity()
                df_si = self.loader.load_sample_info()
                return plot_group_comparison_box(df_cd, group_df=df_si)
            elif chart_key == "venn":
                df_fa = self.loader.load_functional_annotation()
                sample_ids = df_fa["Sample"].unique()[:3].tolist()
                if len(sample_ids) >= 2:
                    return plot_venn_diagram(df_fa, sample_ids=sample_ids)
            elif chart_key == "expansion":
                df_se = self.loader.load_significant_expansion()
                return plot_significant_expansion(df_se)
        except Exception as e:
            print(f"  [SKIP] {chart_key}: {e}")
        return None

    def _assemble_single_sample_charts(self, sample_id: str) -> dict[str, list[dict]]:
        """Generate all single-sample charts, split by section.

        Returns dict with keys 'sec41', 'sec42', 'sec43'."""
        all_cards = {"sec41": [], "sec42": [], "sec43": []}
        for static_name, chart_key, meta in self.SINGLE_CHART_SPECS:
            title = meta["title"]
            caption = meta["caption"]
            # Determine section from title prefix
            if "4.1" in title:
                sec = "sec41"
            elif "4.2" in title:
                sec = "sec42"
            else:
                sec = "sec43"

            if static_name:
                img_html = self._png_to_base64_img(static_name)
                if img_html is not None:
                    all_cards[sec].append({"title": title, "html": img_html, "caption": caption})
                    continue

            fig = self._make_single_chart(chart_key, sample_id)
            if fig is not None:
                all_cards[sec].append(self._fig_to_card(fig, title=title, caption=caption))
            else:
                all_cards[sec].append({
                    "title": title,
                    "html": f'<p style="color:#999;font-style:italic;">图表生成失败，请检查输入数据</p>',
                    "caption": caption,
                })

        return all_cards

    def _assemble_multi_sample_charts(self) -> list[dict]:
        """Generate all multi-sample charts, preferring static images."""
        cards = []
        for static_name, chart_key, meta in self.MULTI_CHART_SPECS:
            title = meta["title"]
            caption = meta["caption"]

            if static_name:
                img_html = self._png_to_base64_img(static_name)
                if img_html is not None:
                    cards.append({"title": title, "html": img_html, "caption": caption})
                    continue

            fig = self._make_multi_chart(chart_key)
            if fig is not None:
                cards.append(self._fig_to_card(fig, title=title, caption=caption))
            else:
                cards.append({
                    "title": title,
                    "html": f'<p style="color:#999;font-style:italic;">图表生成失败，请检查输入数据</p>',
                    "caption": caption,
                })

        return cards

    def assemble(
        self,
        output_path: str | Path = "output/report.html",
        sample_id: str | None = None,
        skip_single: bool = False,
        skip_multi: bool = False,
    ) -> str:
        """Assemble the complete report.

        Parameters
        ----------
        output_path : str or Path
            Output HTML file path.
        sample_id : str or None
            Specific sample for single-sample analysis. If None, uses first sample.
        skip_single : bool
            If True, skip single-sample analysis section.
        skip_multi : bool
            If True, skip multi-sample analysis section.

        Returns
        -------
        str
            Path to the generated report.
        """
        print("Loading data...")
        data = self.loader.load_all()

        # Pick sample for single-sample analysis
        if sample_id is None:
            sample_id = data["sample_info"]["SampleID"].iloc[0]
        print(f"Single-sample analysis: {sample_id}")

        # Generate single-sample charts
        single_cards = {"sec41": [], "sec42": [], "sec43": []}
        if not skip_single:
            print("Generating single-sample charts...")
            single_cards = self._assemble_single_sample_charts(sample_id)
            total = sum(len(v) for v in single_cards.values())
            print(f"  Generated {total} charts (4.1:{len(single_cards['sec41'])} 4.2:{len(single_cards['sec42'])} 4.3:{len(single_cards['sec43'])})")

        # Generate multi-sample charts
        multi_cards = []
        if not skip_multi:
            print("Generating multi-sample charts...")
            multi_cards = self._assemble_multi_sample_charts()
            print(f"  Generated {len(multi_cards)} charts")

        # Data tables
        sample_info_html = self._df_to_html_table(
            data["sample_info"][["SampleID", "PatientID", "Group", "Type", "SampleType"]],
            max_rows=30,
        )
        qc_html = self._df_to_html_table(
            data["qc_table"][["Sample", "rawYield", "cleanYield", "rawReads", "cleanReads",
                               "cleanQ30", "cleanRatio", "clonalReads", "clonalRatio"]],
            max_rows=30,
        )

        # Render
        print("Rendering HTML report...")
        html = self.renderer.render(
            report_title=self.report_title,
            company=self.company,
            sample_info_table=sample_info_html,
            qc_table=qc_html,
            single_sample_charts=single_cards,
            multi_sample_charts=multi_cards,
            version=self.version,
        )

        # Write
        self.renderer.write(html, output_path)
        print(f"Report saved to: {Path(output_path).resolve()}")
        print(f"  Size: {Path(output_path).stat().st_size / 1024 / 1024:.1f} MB")

        return str(Path(output_path).resolve())
