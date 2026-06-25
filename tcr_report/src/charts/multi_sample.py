"""Multi-sample comparison charts for TCR report.

Charts covering sections 5.1–5.6 of the reference report:
- PCA clustering
- V/J gene clustering heatmaps
- Clone similarity heatmap (Morisita index)
- Group comparison box plots
- Venn diagram (matplotlib fallback)
- Significantly expanded TCR analysis
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

from src.charts.utils import (
    MORANDI_COLORS,
    OTHER_GREY,
    apply_tcr_theme,
    get_color_palette,
    hex_to_rgba,
)


# ══════════════════════════════════════════════════════════════════════════
# 5.1 PCA Analysis
# ══════════════════════════════════════════════════════════════════════════

def plot_pca_scatter(
    df: pd.DataFrame,
    *,
    sample_col: str = "Sample",
    value_col: str = "Frequency",
    gene_col: str = "V_gene",
    group_df: pd.DataFrame | None = None,
    group_col: str = "Group",
    sample_id_col: str = "SampleID",
    output_path: str | None = None,
    title: str = "PCA Clustering Analysis",
    n_components: int = 2,
) -> go.Figure:
    """PCA scatter plot from V gene frequency matrix.

    Performs PCA on samples × V-genes matrix and plots PC1 vs PC2.

    Parameters
    ----------
    df : pd.DataFrame
        V gene frequency data with Sample, V_gene, Frequency columns.
    group_df : pd.DataFrame or None
        Sample info with Group column for coloring. If None, all same color.
    """
    # Pivot to samples × genes matrix
    pivot = df.pivot_table(
        index=sample_col, columns=gene_col, values=value_col, fill_value=0
    )
    # PCA via SVD
    from numpy.linalg import svd
    X = pivot.values
    X_centered = X - X.mean(axis=0)
    U, S, Vt = svd(X_centered, full_matrices=False)
    pc_scores = U * S  # n_samples × n_components

    pc1 = pc_scores[:, 0]
    pc2 = pc_scores[:, 1] if pc_scores.shape[1] > 1 else np.zeros_like(pc1)
    var_explained = (S ** 2) / (S ** 2).sum()
    pc1_var = var_explained[0] * 100
    pc2_var = var_explained[1] * 100 if len(var_explained) > 1 else 0

    # Group coloring
    if group_df is not None:
        group_map = dict(zip(group_df[sample_id_col], group_df[group_col]))
        sample_groups = [group_map.get(s, "Unknown") for s in pivot.index]
    else:
        sample_groups = ["Sample"] * len(pivot.index)

    unique_groups = sorted(set(sample_groups))
    colors = get_color_palette(len(unique_groups))
    group_color_map = dict(zip(unique_groups, colors))

    fig = go.Figure()

    for group in unique_groups:
        mask = [g == group for g in sample_groups]
        fig.add_trace(go.Scatter(
            x=pc1[mask],
            y=pc2[mask],
            mode="markers",
            name=group,
            marker=dict(
                color=group_color_map[group],
                size=10,
                line=dict(color="white", width=1),
            ),
            text=[s for s, m in zip(pivot.index, mask) if m],
            hovertemplate="%{text}<br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<extra></extra>",
        ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Based on V gene frequencies</sub>",
            x=0, xref="paper",
        ),
        xaxis_title=f"PC1 ({pc1_var:.1f}%)",
        yaxis_title=f"PC2 ({pc2_var:.1f}%)",
        legend=dict(x=1.02, y=1),
        height=600,
        width=800,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 5.2 V/J Gene Clustering Heatmaps
# ══════════════════════════════════════════════════════════════════════════

def plot_gene_clustermap(
    df: pd.DataFrame,
    *,
    sample_col: str = "Sample",
    gene_col: str = "V_gene",
    freq_col: str = "Frequency",
    output_path: str | None = None,
    gene_type: str = "V",
    top_genes: int = 30,
    title: str | None = None,
    zscore: bool = True,
) -> go.Figure:
    """Clustered heatmap: gene usage across samples with dendrograms.

    Parameters
    ----------
    df : pd.DataFrame
        Gene frequency data.
    gene_type : str
        "V" or "J" for axis labels.
    top_genes : int
        Top N genes by total frequency to include.
    zscore : bool
        If True, Z-score normalize rows before plotting.
    """
    # Pivot and select top genes
    pivot = df.pivot_table(
        index=gene_col, columns=sample_col, values=freq_col, fill_value=0
    )
    gene_totals = pivot.sum(axis=1).sort_values(ascending=False)
    top = gene_totals.head(top_genes).index
    pivot = pivot.loc[top]

    # Z-score normalize rows
    if zscore:
        data = pivot.values
        row_means = data.mean(axis=1, keepdims=True)
        row_stds = data.std(axis=1, keepdims=True) + 1e-10
        data_z = (data - row_means) / row_stds
    else:
        data_z = pivot.values

    # Hierarchical clustering on columns (samples)
    col_linkage = linkage(pdist(data_z.T), method="ward")
    col_dendro = dendrogram(col_linkage, no_plot=True)
    col_order = col_dendro["leaves"]

    # Clustering on rows (genes)
    row_linkage = linkage(pdist(data_z), method="ward")
    row_dendro = dendrogram(row_linkage, no_plot=True)
    row_order = row_dendro["leaves"]

    # Reorder data
    data_ordered = data_z[np.array(row_order)][:, np.array(col_order)]
    x_labels = [pivot.columns[i] for i in col_order]
    y_labels = [pivot.index[i] for i in row_order]

    if title is None:
        title = f"{gene_type} Gene Clustering"

    fig = go.Figure(data=go.Heatmap(
        z=data_ordered,
        x=x_labels,
        y=y_labels,
        colorscale="RdBu_r",
        zmid=0,
        colorbar=dict(title="Z-score", len=0.6),
        hovertemplate=(
            f"{gene_type}: %{{y}}<br>"
            "Sample: %{x}<br>"
            "Z-score: %{z:.2f}<extra></extra>"
        ),
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Z-score normalized, hierarchical clustering</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="Sample",
        yaxis_title=f"TRB{gene_type} Genes",
        height=max(600, len(y_labels) * 25),
        width=max(800, len(x_labels) * 12),
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(tickangle=90, tickfont_size=8, showticklabels=False)
    fig.update_yaxes(tickfont_size=10)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 5.3 Clone Similarity Analysis (Morisita Index)
# ══════════════════════════════════════════════════════════════════════════

def plot_clone_similarity_heatmap(
    df: pd.DataFrame,
    *,
    sample_col: str = "Sample",
    gene_col: str = "V_gene",
    freq_col: str = "Frequency",
    output_path: str | None = None,
    title: str = "Clone Similarity Heatmap (Morisita Index)",
) -> go.Figure:
    """Heatmap: pairwise clone similarity between samples using Morisita index.

    Parameters
    ----------
    df : pd.DataFrame
        V gene frequency data across samples.
    """
    pivot = df.pivot_table(
        index=sample_col, columns=gene_col, values=freq_col, fill_value=0
    )

    # Morisita-Horn similarity index
    n = len(pivot)
    sim_matrix = np.zeros((n, n))
    X = pivot.values

    for i in range(n):
        for j in range(i, n):
            xi, xj = X[i], X[j]
            # Morisita-Horn
            numerator = 2 * np.sum(xi * xj)
            denom_a = np.sum(xi ** 2) / (np.sum(xi) ** 2) if np.sum(xi) > 0 else 1
            denom_b = np.sum(xj ** 2) / (np.sum(xj) ** 2) if np.sum(xj) > 0 else 1
            denom = (denom_a + denom_b) * np.sum(xi) * np.sum(xj)
            sim = numerator / denom if denom > 0 else 0
            sim_matrix[i, j] = sim_matrix[j, i] = min(sim, 1.0)

    labels = pivot.index.tolist()

    fig = go.Figure(data=go.Heatmap(
        z=sim_matrix,
        x=labels,
        y=labels,
        colorscale="Reds",
        zmin=0,
        zmax=1,
        colorbar=dict(title="Similarity", len=0.6),
        hovertemplate=(
            "Sample A: %{x}<br>"
            "Sample B: %{y}<br>"
            "Morisita: %{z:.4f}<extra></extra>"
        ),
    ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        height=max(600, len(labels) * 15),
        width=max(800, len(labels) * 15),
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(tickangle=90, tickfont_size=7, showticklabels=False)
    fig.update_yaxes(tickfont_size=7)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 5.4 Group Comparison Box Plots
# ══════════════════════════════════════════════════════════════════════════

def plot_group_comparison_box(
    df: pd.DataFrame,
    *,
    group_df: pd.DataFrame | None = None,
    sample_col: str = "Sample",
    group_col: str = "Group",
    sample_id_col: str = "SampleID",
    metrics: list[str] | None = None,
    output_path: str | None = None,
    title: str = "Group Comparison of Diversity Metrics",
) -> go.Figure:
    """Box plots: Shannon, Simpson, Evenness comparison between sample groups.

    Parameters
    ----------
    df : pd.DataFrame
        Clone diversity data with Sample column.
    group_df : pd.DataFrame
        Sample info with Group column.
    metrics : list or None
        List of metric columns to plot. Default: ShannonIndex, SimpsonIndex, Evenness.
    """
    if metrics is None:
        metrics = ["ShannonIndex", "SimpsonIndex", "Evenness"]

    if group_df is not None:
        group_map = dict(zip(group_df[sample_id_col], group_df[group_col]))
        df = df.copy()
        df["Group"] = df[sample_col].map(group_map)
    else:
        df = df.copy()
        df["Group"] = "All"

    groups = sorted(df["Group"].dropna().unique())
    n_metrics = len(metrics)
    colors = get_color_palette(len(groups))

    fig = go.Figure()

    for i, metric in enumerate(metrics):
        if metric not in df.columns:
            continue
        for j, group in enumerate(groups):
            values = df[df["Group"] == group][metric].dropna()
            offset = i * (len(groups) + 1) + j  # position groups side by side

            fig.add_trace(go.Box(
                y=values,
                name=f"{group}",
                x=[metric] * len(values),
                marker_color=colors[j],
                legendgroup=group,
                showlegend=(i == 0),  # Show legend only once per group
                offsetgroup=group,
                hovertemplate=f"{metric}: %{{y:.4f}}<extra>{group}</extra>",
            ))

    fig.update_layout(
        boxmode="group",
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        yaxis_title="Metric Value",
        xaxis_title="",
        height=500,
        legend=dict(x=1.02, y=1),
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 5.5 Venn Diagram (Matplotlib fallback)
# ══════════════════════════════════════════════════════════════════════════

def plot_venn_diagram(
    df: pd.DataFrame,
    *,
    sample_ids: list[str],
    output_path: str | None = None,
    title: str = "Clone Overlap Analysis",
    clone_id_col: str = "cdr3aa",
    sample_col: str = "Sample",
) -> go.Figure:
    """Venn diagram showing clone overlap between 2–5 samples.

    Uses matplotlib_venn for rendering, then returns as Plotly Figure
    with an embedded static image.

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data with clone CDR3 sequences.
    sample_ids : list
        2–5 sample IDs to compare.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2, venn3
    import io
    import base64

    if len(sample_ids) < 2 or len(sample_ids) > 5:
        raise ValueError("sample_ids must have 2–5 samples")

    # Extract clone sets per sample
    sets = {}
    for sid in sample_ids:
        df_s = df[df[sample_col] == sid]
        sets[sid] = set(df_s[clone_id_col].dropna().unique())

    n = len(sample_ids)
    fig_mpl, ax = plt.subplots(figsize=(8, 6))

    if n == 2:
        venn2([sets[s] for s in sample_ids], set_labels=sample_ids, ax=ax)
    elif n == 3:
        venn3([sets[s] for s in sample_ids], set_labels=sample_ids, ax=ax)
    else:
        # For 4-5 samples, use upset-style bar chart as alternative
        plt.close(fig_mpl)
        return _plot_upset_alternative(sets, title=title, output_path=output_path)

    ax.set_title(title, fontweight="bold", fontsize=14)
    plt.tight_layout()

    # Save to buffer and encode
    buf = io.BytesIO()
    fig_mpl.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig_mpl)

    # Embed in Plotly figure
    fig = go.Figure()
    fig.add_layout_image(
        source=f"data:image/png;base64,{img_base64}",
        x=0.5, y=0.5,
        xref="paper", yref="paper",
        sizex=1, sizey=1,
        xanchor="center", yanchor="middle",
    )
    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        height=600,
        width=700,
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def _plot_upset_alternative(
    sets: dict[str, set],
    *,
    title: str = "Clone Overlap Analysis",
    output_path: str | None = None,
) -> go.Figure:
    """Alternative to Venn for 4–5 samples: bar chart with overlap info."""
    sample_names = list(sets.keys())
    counts = [len(s) for s in sets.values()]

    colors = get_color_palette(len(sample_names))

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=sample_names,
        y=counts,
        marker_color=colors,
        text=counts,
        textposition="outside",
        hovertemplate="%{x}: %{y:,} unique clones<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Unique clone counts per sample (Venn not available for >3 groups)</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="Sample",
        yaxis_title="Unique Clone Count",
        height=500,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 5.6 Significantly Expanded TCR
# ══════════════════════════════════════════════════════════════════════════

def plot_significant_expansion(
    df: pd.DataFrame,
    *,
    output_path: str | None = None,
    fc_col: str = "FC",
    p_adj_col: str = "p.adj",
    freq_ref_col: str = "freq_Ref",
    freq_amp_col: str = "freq_Amp",
    name_col: str = "name",
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
    title: str = "Significantly Expanded TCR Clones",
) -> go.Figure:
    """Scatter + line chart showing TCR clone expansion after antigen stimulation.

    Parameters
    ----------
    df : pd.DataFrame
        Significant expansion data.
    """
    df_sig = df[(df[fc_col].abs() >= fc_threshold) & (df[p_adj_col] <= p_threshold)].copy()
    if df_sig.empty:
        df_sig = df.copy()  # Show all if none pass threshold

    df_sig = df_sig.sort_values(fc_col, ascending=False)

    # Aggregate by name across samples
    agg = df_sig.groupby(name_col).agg(
        avg_freq_ref=(freq_ref_col, "mean"),
        avg_freq_amp=(freq_amp_col, "mean"),
        avg_fc=(fc_col, "mean"),
    ).reset_index()
    agg = agg.sort_values("avg_fc", ascending=False).head(30)

    fig = go.Figure()

    for _, row in agg.iterrows():
        name = row[name_col]
        fig.add_trace(go.Scatter(
            x=["Baseline", "Expanded"],
            y=[row["avg_freq_ref"], row["avg_freq_amp"]],
            mode="lines+markers",
            name=name[:15],
            marker=dict(size=8),
            hovertemplate=(
                f"{name}<br>"
                "%{x}: %{y:.6f}<br>"
                f"FC: {row['avg_fc']:.2f}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>FC > {fc_threshold}, p.adj < {p_threshold}</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="Condition",
        yaxis_title="Clone Frequency",
        yaxis_type="log",
        legend=dict(x=1.02, y=1, font_size=8),
        height=550,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig
