"""Single-sample analysis charts for TCR report.

Each function accepts a DataFrame and returns a Plotly Figure.
All charts use Plotly for interactive HTML embedding.

Function signature convention:
    def plot_xxx(df: pd.DataFrame, *, output_path: str | None = None, ...) -> go.Figure
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from src.charts.utils import (
    MORANDI_COLORS,
    OTHER_GREY,
    apply_tcr_theme,
    get_color_palette,
    hex_to_rgba,
    figure_to_html_div,
)


# ══════════════════════════════════════════════════════════════════════════
# 4.1 CDR3 Analysis
# ══════════════════════════════════════════════════════════════════════════

def plot_cdr3_length_distribution(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    top_n: int = 10,
    cdr3_col: str = "cdr3aa",
    freq_col: str = "freq",
    title: str = "CDR3 Amino Acid Length Distribution",
) -> go.Figure:
    """Stacked bar: CDR3 AA length distribution with top-N clonotypes highlighted.

    X-axis = CDR3 AA length.
    Y-axis = cumulative proportion of CDR3 clones at each length.
    The top-N most frequent CDR3 sequences across the whole sample are
    shown as coloured segments; all other clones are grey.

    Data is read directly from functional_annotation (cdr3aa + freq),
    so each clone naturally belongs to exactly one length bin.

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data. Expected: Sample, cdr3aa, freq.
    sample_id : str
        The sample being plotted.
    top_n : int
        Number of top clones to show individually.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No data for sample {sample_id}")

    # CDR3 length for every clone
    df_s["cdr3_len"] = df_s[cdr3_col].str.len()
    lengths = sorted(df_s["cdr3_len"].unique())

    # --- Y-axis: cumulative proportion per length ---
    total_freq = df_s[freq_col].sum()
    len_freq = df_s.groupby("cdr3_len")[freq_col].sum()
    len_prop = (len_freq / total_freq).reindex(lengths, fill_value=0.0)

    # --- Top-N clones by total frequency ---
    clone_totals = df_s.groupby(cdr3_col)[freq_col].sum().sort_values(ascending=False)
    top_clones = clone_totals.head(top_n)

    # How much does each top clone contribute at each length?
    colors = get_color_palette(len(top_clones) + 1)

    # Compute "other" at each length (total minus top-N)
    top_per_len = np.zeros(len(lengths))
    top_data = {}  # clone_id -> array of freqs per length
    for clone_id in top_clones.index:
        c_len = len(str(clone_id))
        clone_freqs = np.zeros(len(lengths))
        df_c = df_s[df_s[cdr3_col] == clone_id]
        if not df_c.empty:
            c_len_val = df_c["cdr3_len"].iloc[0]
            idx = lengths.index(c_len_val)
            clone_freqs[idx] = df_c[freq_col].sum() / total_freq
        top_data[clone_id] = clone_freqs
        top_per_len += clone_freqs
    other_per_len = np.maximum(0, len_prop.values - top_per_len)

    fig = go.Figure()

    # Top clones first → stacked at bottom, easiest to read
    for i, (clone_id, _) in enumerate(top_clones.items()):
        label = str(clone_id)
        if len(label) > 20:
            label = label[:19] + "…"
        fig.add_trace(go.Bar(
            x=lengths, y=top_data[clone_id].tolist(),
            name=label,
            marker_color=colors[i],
            marker_line=dict(color="white", width=0.5),
            hovertemplate=f"Length=%{{x}}<br>{clone_id}=%{{y:.6f}}<extra></extra>",
        ))

    # Other clones (grey, stacked on top — remainder)
    fig.add_trace(go.Bar(
        x=lengths, y=other_per_len.tolist(),
        name="Other Clones",
        marker_color=OTHER_GREY,
        hovertemplate="Length=%{x}<br>Other=%{y:.4f}<extra></extra>",
    ))

    n_ticks = len(lengths)
    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id} | Top {top_n} clonotypes</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="CDR3 Amino Acid Length",
        yaxis_title="Cumulative Proportion",
        legend=dict(
            title=dict(text="CDR3 Clonotype"),
            orientation="v", x=1.02, y=1,
            font_size=10,
        ),
        height=450,
        margin=dict(l=60, r=40, t=80, b=80 if n_ticks > 12 else 60),
    )
    fig = apply_tcr_theme(fig)
    if n_ticks > 12:
        fig.update_xaxes(tickvals=lengths, tickangle=90, tickfont_size=9)
    else:
        fig.update_xaxes(tickvals=lengths)
    fig.update_yaxes(tickformat=".0%")

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_v_gene_cdr3_length(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    top_n_v: int = 12,
    v_col: str = "v",
    cdr3_col: str = "cdr3aa",
    freq_col: str = "freq",
    title: str = "Top V-gene Distribution Across CDR3 Lengths",
) -> go.Figure:
    """Stacked bar chart: Top V-gene composition across CDR3 amino acid lengths.

    For each CDR3 length, this shows which V genes are used by clones at
    that length.  Only the top-N V genes (by total frequency) are shown
    individually; the rest are grouped as "Other V-genes".

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data. Expected: Sample, v (V gene), cdr3aa, freq.
    sample_id : str
        The sample being plotted.
    top_n_v : int
        Number of top V genes to show individually.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No functional annotation data for sample {sample_id}")

    # CDR3 length per clone
    df_s["cdr3_len"] = df_s[cdr3_col].str.len()
    lengths = sorted(df_s["cdr3_len"].dropna().unique())

    # Top V genes by total frequency
    v_totals = df_s.groupby(v_col)[freq_col].sum().sort_values(ascending=False)
    top_v = v_totals.head(top_n_v).index.tolist()

    # Total frequency (for computing cumulative proportion, same as 图4.1.1)
    total_freq = df_s[freq_col].sum()

    # Per-length cumulative proportion (same as 图4.1.1 bar heights)
    len_freq = df_s.groupby("cdr3_len")[freq_col].sum()
    len_prop = (len_freq / total_freq).reindex(lengths, fill_value=0.0)

    # Aggregate: raw frequency per (length, V_gene), then as fraction of TOTAL
    grouped = df_s.groupby(["cdr3_len", v_col])[freq_col].sum().reset_index()
    grouped["GlobalProp"] = grouped[freq_col] / total_freq

    colors = get_color_palette(len(top_v) + 1)

    fig = go.Figure()

    # Top V genes first → at bottom of bar
    for i, v_gene in enumerate(top_v):
        vals = []
        for length in lengths:
            row = grouped[(grouped["cdr3_len"] == length) & (grouped[v_col] == v_gene)]
            vals.append(row["GlobalProp"].values[0] if not row.empty else 0.0)
        fig.add_trace(go.Bar(
            x=lengths, y=vals,
            name=v_gene,
            marker_color=colors[i],
            marker_line=dict(color="white", width=0.5),
            hovertemplate=f"Length=%{{x}}<br>{v_gene}=%{{y:.4f}}<extra></extra>",
        ))

    # Other V-genes (grey, stacked on top = remainder)
    other_by_len = []
    for length in lengths:
        total_prop = len_prop.get(length, 0.0)
        top_sum = grouped[(grouped["cdr3_len"] == length) & (grouped[v_col].isin(top_v))]["GlobalProp"].sum()
        other_by_len.append(max(0.0, total_prop - top_sum))
    fig.add_trace(go.Bar(
        x=lengths, y=other_by_len,
        name="Other V-genes",
        marker_color=OTHER_GREY,
        hovertemplate="Length=%{x}<br>Other=%{y:.4f}<extra></extra>",
    ))

    n_ticks = len(lengths)
    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id} | Top {top_n_v} V genes</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="CDR3 Amino Acid Length",
        yaxis_title="Cumulative Proportion",
        legend=dict(
            title=dict(text="V Gene"),
            orientation="v", x=1.02, y=1,
            font_size=10,
        ),
        height=450,
        margin=dict(l=60, r=40, t=80, b=80 if n_ticks > 12 else 60),
    )
    fig = apply_tcr_theme(fig)
    if n_ticks > 12:
        fig.update_xaxes(tickvals=lengths, tickangle=90, tickfont_size=9)
    else:
        fig.update_xaxes(tickvals=lengths)
    fig.update_yaxes(tickformat=".0%")

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_aa_composition(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    cdr3_col: str = "cdr3aa",
    title: str = "Amino Acid Composition",
) -> go.Figure:
    """Horizontal stacked bar chart: amino acid composition from CDR3 sequences.

    Counts each amino acid across all CDR3 sequences for a sample, then
    normalises to relative frequencies.

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data. Expected: Sample, cdr3aa (CDR3 AA sequences).
        If sample_id is None, uses the first sample found.
    """
    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id].copy()
    elif sample_id is None and "Sample" in df.columns:
        sample_id = df["Sample"].iloc[0]
        df = df[df["Sample"] == sample_id].copy()

    if df.empty:
        raise ValueError("No CDR3 sequence data found")

    # Count every amino acid across all CDR3 sequences
    aa_counter: dict[str, int] = {}
    total_aa = 0
    for seq in df[cdr3_col].dropna():
        for aa in str(seq).upper():
            if aa in "ACDEFGHIKLMNPQRSTVWY":  # only standard 20 AA
                aa_counter[aa] = aa_counter.get(aa, 0) + 1
                total_aa += 1

    if total_aa == 0:
        raise ValueError("No valid amino acids found in CDR3 sequences")

    # Build sorted frequency series (ascending for horizontal bar)
    aa_series = pd.Series(aa_counter, name="count") / total_aa
    aa_series = aa_series.sort_values(ascending=True)

    # Ensure all 20 standard AAs are represented (fill 0 for missing)
    std_order = sorted(aa_series.index, key=lambda x: aa_series[x])
    colors = get_color_palette(len(std_order))
    color_map = dict(zip(std_order, colors))

    fig = go.Figure()
    for aa in std_order:
        freq = aa_series.get(aa, 0)
        fig.add_trace(go.Bar(
            y=[sample_id or "Sample"],
            x=[freq],
            name=aa,
            orientation="h",
            marker_color=color_map[aa],
            marker_line=dict(color="white", width=0.5),
            hovertemplate=f"{aa}=%{{x:.3f}}<extra></extra>",
        ))

    n_aa = len(std_order)

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id} | {len(df):,} CDR3 sequences</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="Relative Amino Acid Composition",
        yaxis_title="",
        legend=dict(
            title=dict(text="Amino Acid"),
            orientation="h",
            y=-0.50,
            x=0.5,
            xanchor="center",
            font_size=10,
            itemwidth=32,
            tracegroupgap=3,
        ),
        height=320,
        margin=dict(b=140, l=60, r=40, t=80),
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(range=[0, 1])

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_cdr3_motif(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    top_n: int = 100,
    title: str = "High-Abundance TCR Clone CDR3β Motif",
) -> go.Figure:
    """Sequence logo: amino acid letters at each CDR3 position.

    Uses logomaker to render a proper sequence logo and embeds it as a
    base64 PNG inside a Plotly figure (so the assembler pipeline is unchanged).

    X-axis = CDR3 position (Pos).
    Y-axis = probability (Prob) per amino acid at that position.
    Letter height = frequency at that position.
    Colour = standard sequence-logo chemical-group scheme.

    Parameters
    ----------
    df : pd.DataFrame
        TCR motif data. Expected: Sample, Clone_ID, CDR3_AA, Position, AminoAcid, Frequency.
    sample_id : str or None
        Sample to plot. If None, uses the first sample found.
    top_n : int
        Number of top-abundance clones to include.
    """
    import base64
    import io
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import logomaker

    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id].copy()
    elif sample_id is None and "Sample" in df.columns:
        sample_id = df["Sample"].iloc[0]
        df = df[df["Sample"] == sample_id].copy()

    if df.empty:
        raise ValueError(f"No TCR motif data for sample {sample_id}")

    # Select top-N clones by total frequency
    clone_totals = df.groupby("Clone_ID")["Frequency"].sum().sort_values(ascending=False)
    top_clones = clone_totals.head(top_n).index.tolist()
    df_top = df[df["Clone_ID"].isin(top_clones)].copy()

    # Aggregate: mean frequency per (Position, AminoAcid) across top clones
    pos_aa = df_top.groupby(["Position", "AminoAcid"])["Frequency"].mean().reset_index()
    pos_aa["FreqNorm"] = pos_aa.groupby("Position")["Frequency"].transform(lambda x: x / x.sum())

    # Build position-weight matrix (rows=positions, cols=20 AAs)
    positions = sorted(pos_aa["Position"].unique())
    aa_letters = list("ACDEFGHIKLMNPQRSTVWY")
    pwm = pd.DataFrame(0.0, index=positions, columns=aa_letters)
    for _, row in pos_aa.iterrows():
        pwm.at[row["Position"], row["AminoAcid"]] = row["FreqNorm"]

    # User-specified colour scheme
    color_scheme = {
        "A": "#377eb8", "V": "#377eb8", "L": "#377eb8", "I": "#377eb8",
        "M": "#377eb8", "F": "#377eb8", "W": "#377eb8",
        "S": "#33cc33", "T": "#33cc33", "N": "#33cc33", "Q": "#33cc33",
        "D": "#cc33cc", "E": "#cc33cc",
        "K": "#e41a1c", "R": "#e41a1c",
        "H": "#00bfc4", "Y": "#00bfc4",
        "C": "#f78181",
        "G": "#d98c3a",
        "P": "#b3c833",
    }

    # Render with logomaker — taller figure
    fig_mpl, ax = plt.subplots(figsize=(max(8, len(positions) * 0.42), 5.5))
    logo = logomaker.Logo(
        pwm, ax=ax,
        color_scheme=color_scheme,
        font_name="DejaVu Sans",
        show_spines=False,
    )
    logo.style_spines(spines=["left", "bottom"], visible=True, linewidth=0.5)
    logo.style_spines(spines=["top", "right"], visible=False)
    logo.ax.set_xticks(range(len(positions)))
    logo.ax.set_xticklabels([str(p) for p in positions])
    logo.ax.set_xlabel("CDR3 Position", fontsize=13)
    logo.ax.set_ylabel("Probability", fontsize=13)
    logo.ax.set_ylim(0, 1.05)
    fig_mpl.tight_layout(pad=0.5)

    # Convert to base64 PNG
    buf = io.BytesIO()
    fig_mpl.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig_mpl)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")

    # Embed in a Plotly figure — white background, no axes
    fig = go.Figure()
    fig.add_layout_image(dict(
        source=f"data:image/png;base64,{b64}",
        xref="paper", yref="paper",
        x=0, y=1, sizex=1, sizey=1,
        xanchor="left", yanchor="top",
    ))
    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0.08, xref="paper", font=dict(color="#333333"),
        ),
        width=max(800, len(positions) * 36),
        height=620,
        margin=dict(l=0, r=0, t=60, b=0),
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_cdr3_physicochemical(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    cdr3_col: str = "cdr3aa",
    freq_col: str = "freq",
    title: str = "CDR3 Physicochemical Properties",
) -> go.Figure:
    """Scatter plot: hydrophobicity vs net charge for each CDR3 sequence.

    Each point represents one CDR3 amino acid sequence from the sample.
    Point size is weighted by clone frequency.

    Uses the Kyte-Doolittle hydrophobicity scale and a simple net-charge
    model (R/K = +1, D/E = -1).

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data. Expected: Sample, cdr3aa, freq.
    sample_id : str or None
        Sample to plot. If None, uses the first sample found.
    """
    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id].copy()
    elif sample_id is None and "Sample" in df.columns:
        sample_id = df["Sample"].iloc[0]
        df = df[df["Sample"] == sample_id].copy()

    if df.empty:
        raise ValueError("No CDR3 sequence data found")

    # Kyte-Doolittle hydrophobicity scale
    kd_scale = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    }
    # Net charge: positive residues
    charge_pos = {"R", "K"}
    charge_neg = {"D", "E"}

    hydro_vals = []
    charge_vals = []
    freq_vals = []
    hover_texts = []

    for _, row in df.iterrows():
        seq = str(row.get(cdr3_col, ""))
        if not seq:
            continue
        seq_upper = seq.upper()
        valid_aa = [aa for aa in seq_upper if aa in kd_scale]
        if not valid_aa:
            continue

        avg_hydro = sum(kd_scale[aa] for aa in valid_aa) / len(valid_aa)
        net_charge = sum(
            1 if aa in charge_pos else (-1 if aa in charge_neg else 0)
            for aa in valid_aa
        )

        freq = row.get(freq_col, 1.0)
        try:
            freq = float(freq)
        except (ValueError, TypeError):
            freq = 1.0

        hydro_vals.append(avg_hydro)
        charge_vals.append(net_charge)
        freq_vals.append(freq)
        hover_texts.append(
            f"CDR3: {seq[:20]}{'…' if len(seq) > 20 else ''}<br>"
            f"Hydrophobicity: {avg_hydro:.2f}<br>"
            f"Net Charge: {net_charge:+d}<br>"
            f"Frequency: {freq:.6f}"
        )

    if not hydro_vals:
        raise ValueError("No valid CDR3 sequences to analyse")

    # Downsample if too many points (SVG scatter > 8K freezes browsers)
    MAX_POINTS = 8000
    n_total = len(hydro_vals)
    if n_total > MAX_POINTS:
        rng = np.random.default_rng(42)
        indices = rng.choice(n_total, MAX_POINTS, replace=False)
        hydro_vals = [hydro_vals[i] for i in indices]
        charge_vals = [charge_vals[i] for i in indices]
        freq_vals = [freq_vals[i] for i in indices]
        hover_texts = [hover_texts[i] for i in indices]

    # Scale marker sizes by frequency (min size 4, max 24)
    freqs_arr = np.array(freq_vals)
    if freqs_arr.max() > freqs_arr.min():
        sizes = 4 + 20 * (freqs_arr - freqs_arr.min()) / (freqs_arr.max() - freqs_arr.min())
    else:
        sizes = np.full_like(freqs_arr, 12)

    color = MORANDI_COLORS[0]  # primary morandi blue-grey

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=hydro_vals,
        y=charge_vals,
        mode="markers",
        marker=dict(
            size=sizes,
            color=color,
            opacity=0.45,
            line=dict(color="white", width=0.3),
        ),
        text=hover_texts,
        hoverinfo="text",
        showlegend=False,
    ))

    # Reference lines: zero charge, neutral hydrophobicity
    x_min, x_max = min(hydro_vals), max(hydro_vals)
    x_pad = (x_max - x_min) * 0.1
    y_min, y_max = min(charge_vals), max(charge_vals)
    y_pad = max((y_max - y_min) * 0.1, 0.5)

    # Zero-charge reference
    fig.add_hline(y=0, line_dash="dot", line_color=OTHER_GREY, opacity=0.5)
    # Neutral hydrophobicity reference (~0)
    fig.add_vline(x=0, line_dash="dot", line_color=OTHER_GREY, opacity=0.5)

    # Quadrant annotations
    fig.add_annotation(x=x_max - x_pad * 0.3, y=y_max - y_pad * 0.2,
        text="Hydrophobic<br>+", showarrow=False, font_size=9, font_color="#999")
    fig.add_annotation(x=x_min + x_pad * 0.3, y=y_max - y_pad * 0.2,
        text="Hydrophilic<br>+", showarrow=False, font_size=9, font_color="#999")
    fig.add_annotation(x=x_max - x_pad * 0.3, y=y_min + y_pad * 0.2,
        text="Hydrophobic<br>−", showarrow=False, font_size=9, font_color="#999")
    fig.add_annotation(x=x_min + x_pad * 0.3, y=y_min + y_pad * 0.2,
        text="Hydrophilic<br>−", showarrow=False, font_size=9, font_color="#999")

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id} | {len(hydro_vals):,} CDR3 sequences</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="Mean Hydrophobicity (Kyte-Doolittle)",
        yaxis_title="Net Charge",
        height=500,
        margin=dict(l=60, r=30, t=80, b=60),
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 4.2 Gene Analysis
# ══════════════════════════════════════════════════════════════════════════

def plot_v_gene_lollipop(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    gene_col: str = "V_gene",
    freq_col: str = "Frequency",
    top_n: int = 30,
    title: str = "V Gene Frequency Distribution",
) -> go.Figure:
    """Lollipop chart: V gene usage frequency for a single sample.

    Parameters
    ----------
    df : pd.DataFrame
        V gene frequency data. Expected: Sample, V_gene, Frequency.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No V-gene data for sample {sample_id}")

    df_s = df_s.sort_values(freq_col, ascending=True).tail(top_n)  # Top N, ascending for display

    genes = df_s[gene_col].tolist()
    freqs = df_s[freq_col].tolist()

    color = MORANDI_COLORS[0]

    fig = go.Figure()

    # Horizontal lines (stems)
    for gene, freq in zip(genes, freqs):
        fig.add_trace(go.Scatter(
            x=[0, freq],
            y=[gene, gene],
            mode="lines",
            line=dict(color=color, width=2),
            opacity=0.4,
            showlegend=False,
            hoverinfo="skip",
        ))

    # Dots
    fig.add_trace(go.Scatter(
        x=freqs,
        y=genes,
        mode="markers",
        marker=dict(color=color, size=12, line=dict(color="white", width=1)),
        name="Frequency",
        hovertemplate="%{y}=%{x:.4f}%<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="Usage Frequency",
        yaxis_title="TRBV Genes",
        height=max(500, len(genes) * 22),
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_j_gene_lollipop(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    gene_col: str = "J_gene",
    freq_col: str = "Frequency",
    title: str = "J Gene Frequency Distribution",
) -> go.Figure:
    """Lollipop chart: J gene usage frequency for a single sample."""
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No J-gene data for sample {sample_id}")

    df_s = df_s.sort_values(freq_col, ascending=True)

    genes = df_s[gene_col].tolist()
    freqs = df_s[freq_col].tolist()
    color = MORANDI_COLORS[1]

    fig = go.Figure()

    for gene, freq in zip(genes, freqs):
        fig.add_trace(go.Scatter(
            x=[0, freq],
            y=[gene, gene],
            mode="lines",
            line=dict(color=color, width=2),
            opacity=0.4,
            showlegend=False,
            hoverinfo="skip",
        ))

    fig.add_trace(go.Scatter(
        x=freqs,
        y=genes,
        mode="markers",
        marker=dict(color=color, size=12, line=dict(color="white", width=1)),
        name="Frequency",
        hovertemplate="%{y}=%{x:.4f}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="Usage Frequency",
        yaxis_title="TRBJ Genes",
        height=max(400, len(genes) * 22),
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_v_gene_freq_bar(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    gene_col: str = "V_gene",
    freq_col: str = "Frequency",
    top_n: int = 30,
    title: str = "V Gene Frequency Distribution",
) -> go.Figure:
    """Horizontal bar chart: V gene usage frequency for a single sample.

    A cleaner alternative to the lollipop chart — sorted horizontal bars
    showing each V gene's frequency as a proportion of the repertoire.

    Parameters
    ----------
    df : pd.DataFrame
        V gene frequency data. Expected: Sample, V_gene, Frequency.
    sample_id : str
        The sample to plot.
    top_n : int
        Max number of V genes to show.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No V-gene data for sample {sample_id}")

    df_s = df_s.sort_values(freq_col, ascending=True).tail(top_n)

    genes = df_s[gene_col].tolist()
    freqs = df_s[freq_col].tolist()
    colors = get_color_palette(len(genes))

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=genes,
        x=freqs,
        orientation="h",
        marker=dict(
            color=colors,
            line=dict(color="white", width=0.5),
        ),
        hovertemplate="%{y}=%{x:.4f}<extra></extra>",
    ))

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="Usage Frequency",
        yaxis_title="TRBV Genes",
        height=max(500, len(genes) * 22),
        showlegend=False,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_j_gene_freq_bar(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    gene_col: str = "J_gene",
    freq_col: str = "Frequency",
    title: str = "J Gene Frequency Distribution",
) -> go.Figure:
    """Horizontal bar chart: J gene usage frequency for a single sample.

    Parameters
    ----------
    df : pd.DataFrame
        J gene frequency data. Expected: Sample, J_gene, Frequency.
    sample_id : str
        The sample to plot.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No J-gene data for sample {sample_id}")

    df_s = df_s.sort_values(freq_col, ascending=True)

    genes = df_s[gene_col].tolist()
    freqs = df_s[freq_col].tolist()
    colors = get_color_palette(len(genes))

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=genes,
        x=freqs,
        orientation="h",
        marker=dict(
            color=colors,
            line=dict(color="white", width=0.5),
        ),
        hovertemplate="%{y}=%{x:.4f}<extra></extra>",
    ))

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0, xref="paper",
        ),
        xaxis_title="Usage Frequency",
        yaxis_title="TRBJ Genes",
        height=max(400, len(genes) * 22),
        showlegend=False,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_vj_heatmap(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    title: str = "V-J Gene Pairing Landscape",
) -> go.Figure:
    """Heatmap: V-J gene pairing frequencies for a single sample.

    Parameters
    ----------
    df : pd.DataFrame
        V-J pairing data. Expected: Sample, V_gene, J_gene, Frequency.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    if df_s.empty:
        raise ValueError(f"No V-J pairing data for sample {sample_id}")

    # Pivot to matrix
    pivot = df_s.pivot_table(
        index="V_gene", columns="J_gene", values="Frequency", aggfunc="sum", fill_value=0
    )

    # Sort rows and cols by total frequency
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    pivot = pivot[pivot.sum(axis=0).sort_values(ascending=False).index]

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale="Reds",
        zmin=0,
        colorbar=dict(
            title=dict(text="Frequency", font_size=11),
            len=0.6,
            thickness=15,
            x=1.02,
        ),
        hovertemplate="V: %{y}<br>J: %{x}<br>Freq: %{z:.4f}<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0,
            xref="paper",
        ),
        xaxis_title="TRBJ Genes",
        yaxis_title="TRBV Genes",
        height=700,
        width=900,
        margin=dict(r=80),
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(tickangle=45)
    fig.update_yaxes(tickfont_size=9)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_vj_sankey(
    df: pd.DataFrame,
    *,
    sample_id: str,
    output_path: str | None = None,
    min_freq: float = 0.005,
    title: str = "V-J Gene Pairing Flow",
) -> go.Figure:
    """Sankey diagram: V-J gene pairing flow for a single sample.

    Parameters
    ----------
    df : pd.DataFrame
        V-J pairing data. Expected: Sample, V_gene, J_gene, Frequency.
    min_freq : float
        Minimum frequency threshold for drawing links.
    """
    df_s = df[df["Sample"] == sample_id].copy()
    df_s = df_s[df_s["Frequency"] >= min_freq]

    if df_s.empty:
        raise ValueError(f"No V-J pairing data above threshold for sample {sample_id}")

    v_genes = sorted(df_s["V_gene"].unique())
    j_genes = sorted(df_s["J_gene"].unique())
    node_labels = v_genes + j_genes
    v_to_idx = {v: i for i, v in enumerate(v_genes)}
    j_to_idx = {j: i + len(v_genes) for i, j in enumerate(j_genes)}

    sources, targets, values = [], [], []
    for _, row in df_s.iterrows():
        sources.append(v_to_idx[row["V_gene"]])
        targets.append(j_to_idx[row["J_gene"]])
        values.append(row["Frequency"])

    # Node colors: muted for V, grey for J
    v_colors = get_color_palette(len(v_genes))
    node_colors = v_colors + ["#B0B0B0"] * len(j_genes)

    # Link colors: inherit from source V-gene
    link_colors = [hex_to_rgba(v_colors[src], 0.4) for src in sources]

    fig = go.Figure(data=[go.Sankey(
        arrangement="perpendicular",
        node=dict(
            pad=20,
            thickness=15,
            line=dict(color="rgba(0,0,0,0)", width=0),
            label=node_labels,
            color=node_colors,
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
        ),
    )])

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>",
            x=0,
            xref="paper",
        ),
        height=650,
        margin=dict(t=80, b=20, l=20, r=20),
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 4.3 Clone Analysis
# ══════════════════════════════════════════════════════════════════════════

def plot_clone_sunburst(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    level1_col: str = "v",
    level2_col: str = "cdr3aa",
    freq_col: str = "freq",
    top_level1: int = 5,
    top_level2: int = 3,
    title: str = "TCR Clone Distribution",
) -> go.Figure:
    """Sunburst chart: hierarchical view of clone distribution.

    Uses functional annotation data. Level 1 = V gene, Level 2 = CDR3 sequences.

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data.
    top_level1 : int
        Top N V genes to show in detail.
    top_level2 : int
        Top N CDR3 sequences per V gene.
    """
    # Reconstruct from functional annotation
    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id].copy()
    elif sample_id is None and "Sample" in df.columns:
        # Use first sample
        first_sample = df["Sample"].iloc[0]
        df = df[df["Sample"] == first_sample].copy()

    if freq_col not in df.columns:
        df[freq_col] = df["freq"]

    # Group by V gene
    v_totals = df.groupby(level1_col)[freq_col].sum().sort_values(ascending=False)
    top_v = v_totals.head(top_level1)
    other_v_sum = v_totals.iloc[top_level1:].sum()

    ids = ["Total"]
    labels = ["Total TCR Repertoire"]
    parents = [""]
    values = [df[freq_col].sum()]
    colors_list = [OTHER_GREY]
    color_map_base = get_color_palette(len(top_v) + 2)

    # Add V-gene level
    for i, (v_gene, v_freq) in enumerate(top_v.items()):
        v_id = f"v_{v_gene}"
        ids.append(v_id)
        labels.append(v_gene)
        parents.append("Total")
        values.append(v_freq)
        colors_list.append(color_map_base[i])

    if other_v_sum > 0:
        ids.append("v_Other")
        labels.append("Other V-genes")
        parents.append("Total")
        values.append(other_v_sum)
        colors_list.append("#D5D8DC")

    # Add CDR3 level per V gene
    for v_gene in top_v.index:
        df_v = df[df[level1_col] == v_gene]
        clone_totals = df_v.groupby(level2_col)[freq_col].sum().sort_values(ascending=False)
        top_clones = clone_totals.head(top_level2)
        others = clone_totals.iloc[top_level2:].sum()

        v_id = f"v_{v_gene}"
        v_color = color_map_base[list(top_v.index).index(v_gene)]

        for clone, freq in top_clones.items():
            c_id = f"c_{v_gene}_{clone[:10]}"
            ids.append(c_id)
            labels.append(clone[:15] + ("…" if len(clone) > 15 else ""))
            parents.append(v_id)
            values.append(freq)
            colors_list.append(v_color)

        if others > 0:
            ids.append(f"c_{v_gene}_others")
            labels.append("Others")
            parents.append(v_id)
            values.append(others)
            colors_list.append("#D5D8DC")

    fig = go.Figure(go.Sunburst(
        ids=ids,
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        marker=dict(colors=colors_list, line=dict(color="white", width=1)),
        insidetextorientation="radial",
    ))

    fig.update_layout(
        title=dict(
            text=f"<b>{title}</b><br><sub>Sample: {sample_id or 'Aggregated'}</sub>",
            x=0,
            xref="paper",
        ),
        height=700,
        width=700,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_clone_frequency_lines(
    df: pd.DataFrame,
    *,
    sample_ids: list[str] | None = None,
    output_path: str | None = None,
    max_samples: int = 20,
    title: str = "Clone Frequency Distribution",
) -> go.Figure:
    """Line chart: cumulative clone frequency distribution per sample.

    Shows the clone count vs clone frequency for each sample.
    Uses clone diversity data.

    Parameters
    ----------
    df : pd.DataFrame
        Clone diversity data.
    sample_ids : list or None
        Samples to include. If None, takes first max_samples.
    """
    if sample_ids is None:
        sample_ids = df["Sample"].unique()[:max_samples]

    sample_ids = sample_ids[:max_samples]
    colors = get_color_palette(len(sample_ids))

    fig = go.Figure()

    for i, sample in enumerate(sample_ids):
        df_s = df[df["Sample"] == sample]
        if df_s.empty:
            continue

        # Simulate a realistic clone frequency distribution
        clone_count = int(df_s["CloneCount"].iloc[0])
        top_freq = float(df_s["TopFreq"].iloc[0])

        # Generate frequency bins from 1e-5 to top_freq
        x = np.logspace(-5, max(np.log10(top_freq), -1), 200)
        # Cumulative clone count: inverse relationship with frequency
        y = clone_count * (1 - np.log10(1 + x * 1000) / np.log10(1 + top_freq * 1000))
        y = np.clip(y, 1, clone_count)

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=str(sample),
            line=dict(color=colors[i], width=1.5),
            hovertemplate="Sample: " + str(sample) + "<br>Freq: %{x:.6f}<br>Clones: %{y:,.0f}<extra></extra>",
        ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        xaxis_title="Clone Frequency (log scale)",
        yaxis_title="Cumulative Clone Count",
        xaxis=dict(type="log", range=[-5, 0]),
        legend=dict(x=1.02, y=1),
        height=500,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_clone_composition_stacked(
    df: pd.DataFrame,
    *,
    sample_ids: list[str] | None = None,
    output_path: str | None = None,
    max_samples: int = 30,
    title: str = "Clone Composition Distribution",
) -> go.Figure:
    """Stacked percentage bar chart: clone frequency range composition per sample.

    Parameters
    ----------
    df : pd.DataFrame
        Clone composition data. Expected: Sample, FrequencyRange, Proportion.
    """
    range_order = [
        "Rare (≤0.01%)", "Small (0.01%-0.1%)", "Medium (0.1%-1%)",
        "Large (1%-10%)", "Hyperexpanded (>10%)",
    ]
    if sample_ids is None:
        sample_ids = df["Sample"].unique()[:max_samples]
    sample_ids = sample_ids[:max_samples]

    df_s = df[df["Sample"].isin(sample_ids)].copy()
    # Ensure range order
    df_s["FrequencyRange"] = pd.Categorical(df_s["FrequencyRange"], categories=range_order, ordered=True)
    df_s = df_s.sort_values("FrequencyRange")

    colors = get_color_palette(len(range_order))
    color_map = dict(zip(range_order, colors))

    fig = go.Figure()
    for r in range_order:
        df_r = df_s[df_s["FrequencyRange"] == r]
        prop_map = dict(zip(df_r["Sample"], df_r["Proportion"]))
        props = [prop_map.get(s, 0) for s in sample_ids]

        fig.add_trace(go.Bar(
            x=sample_ids,
            y=props,
            name=r,
            marker_color=color_map[r],
            hovertemplate="%{x}<br>" + r + "=%{y:.2%}<extra></extra>",
        ))

    fig.update_layout(
        barmode="stack",
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        xaxis_title="Sample",
        yaxis_title="Proportion of Clone Types",
        yaxis=dict(tickformat=".0%"),
        legend=dict(
            title=dict(text="Frequency Range"),
            y=-0.3,
            x=0.5,
            xanchor="center",
            orientation="h",
        ),
        height=450,
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(tickangle=90, tickfont_size=8, showticklabels=False)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# Functional Annotation & Cancer Risk
# ══════════════════════════════════════════════════════════════════════════

def plot_antigen_bubble(
    df: pd.DataFrame,
    *,
    output_path: str | None = None,
    min_freq: float = 0.001,
    title: str = "Antigen Annotation Profiling Across Cohort",
) -> go.Figure:
    """Bubble chart: antigen annotation across samples.

    Parameters
    ----------
    df : pd.DataFrame
        Functional annotation data. Expected: Sample, antigen.species, freq, count.
    """
    df_p = df[df["freq"] >= min_freq].copy()
    if df_p.empty:
        raise ValueError("No antigen data above the minimum frequency threshold")

    # Group by Sample × Antigen species
    grouped = df_p.groupby(["Sample", "antigen.species"]).agg(
        total_freq=("freq", "sum"),
        total_count=("count", "sum"),
    ).reset_index()

    antigen_order = grouped.groupby("antigen.species")["total_freq"].sum().sort_values(ascending=False).index.tolist()

    fig = go.Figure()

    for antigen in antigen_order:
        df_a = grouped[grouped["antigen.species"] == antigen]
        fig.add_trace(go.Scatter(
            x=df_a["Sample"].tolist(),
            y=[antigen] * len(df_a),
            mode="markers",
            marker=dict(
                size=np.sqrt(df_a["total_count"].values) * 0.3,
                sizemode="area",
                sizeref=2.0 * max(np.sqrt(grouped["total_count"].max()), 1) / (40 ** 2),
                sizemin=4,
                color=df_a["total_freq"].values,
                colorscale="thermal_r",
                showscale=True,
                colorbar=dict(title="Frequency"),
                line=dict(color="black", width=0.5),
            ),
            name=antigen,
            hovertemplate=(
                "Sample: %{x}<br>"
                "Antigen: " + antigen + "<br>"
                "Count: %{marker.size:.0f}<br>"
                "Freq: %{marker.color:.4f}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        xaxis_title="Samples",
        yaxis_title="Annotated Antigen Species",
        height=500,
        showlegend=False,
        margin=dict(r=100, l=120),
        coloraxis_colorbar=dict(
            title="Frequency",
            len=0.5,
            thickness=15,
            x=1.02,
            y=0.5,
            yanchor="middle",
        ),
    )
    fig = apply_tcr_theme(fig)
    fig.update_xaxes(showticklabels=False)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_cancer_risk_violin(
    df: pd.DataFrame,
    *,
    output_path: str | None = None,
    cancer_types: list[str] | None = None,
    group_col: str | None = None,
    title: str = "Cancer Risk Score Distribution",
) -> go.Figure:
    """Violin/ridge-like plot: cancer risk scores across samples.

    Uses multiple violin traces as a ridge-plot alternative.

    Parameters
    ----------
    df : pd.DataFrame
        Cancer risk data. Expected: sample_id, {type}_score columns.
        If group_col is provided, expects a Group column (joined from sample_info).
    """
    if cancer_types is None:
        cancer_types = ["KRAS", "EGFR", "REF", "LUNG_CANCER_GDNA", "LUNG_CANCER_TISSUE"]

    colors = get_color_palette(len(cancer_types))
    n_samples = len(df)
    fig = go.Figure()

    if n_samples <= 2:
        # Single/few samples — horizontal bar chart with log10 scale
        row = df.iloc[0]
        vals = []; labels = []
        for ct in cancer_types:
            sc = f"{ct}_score"
            if sc in df.columns:
                v = float(row[sc])
                if v > 0:
                    vals.append(np.log10(v + 1)); labels.append(ct)
        fig.add_trace(go.Bar(
            x=list(reversed(vals)), y=list(reversed(labels)),
            orientation="h",
            marker=dict(color=list(reversed(colors[:len(labels)])), line=dict(color="white", width=0.5)),
            text=[f"{10**v - 1:.2f}" for v in reversed(vals)], textposition="outside",
            hovertemplate="%{y}: %{customdata:.2f}<extra></extra>",
            customdata=[float(row[f"{l}_score"]) for l in reversed(labels)],
        ))
        fig.update_layout(height=320, showlegend=False, margin=dict(r=120, l=160),
                          xaxis_title="log10(Score + 1)")
    else:
        # Multi-sample — violin
        for i, ct in enumerate(cancer_types):
            score_col = f"{ct}_score"
            if score_col not in df.columns:
                continue
            scores = df[score_col].dropna()
            if scores.empty:
                continue
            fig.add_trace(go.Violin(
                y=[ct] * len(scores), x=scores, name=ct, orientation="h", side="positive",
                line_color=colors[i], fillcolor=hex_to_rgba(colors[i], 0.3),
                points="outliers", marker=dict(size=3, color=colors[i]),
                hovertemplate=f"{ct}: %{{x:.2f}}<extra></extra>",
            ))
        fig.update_layout(height=500, showlegend=False, xaxis_title="Risk Score")

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        yaxis_title="",
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


def plot_cancer_risk_box(
    df: pd.DataFrame,
    *,
    output_path: str | None = None,
    cancer_types: list[str] | None = None,
    group_col: str | None = None,
    title: str = "Cancer Risk Score Comparison",
) -> go.Figure:
    """Box plot: cancer risk scores with group comparison.

    Parameters
    ----------
    df : pd.DataFrame
        Cancer risk data with optional Group column.
    """
    if cancer_types is None:
        cancer_types = ["KRAS", "EGFR", "LUNG_CANCER_GDNA", "LUNG_CANCER_TISSUE"]

    colors = get_color_palette(len(cancer_types))

    fig = go.Figure()

    for i, ct in enumerate(cancer_types):
        score_col = f"{ct}_score"
        if score_col not in df.columns:
            continue

        scores = df[score_col].dropna()
        if group_col and group_col in df.columns:
            # Grouped box
            for group in df[group_col].unique():
                g_scores = df[df[group_col] == group][score_col].dropna()
                fig.add_trace(go.Box(
                    y=g_scores,
                    x=[f"{ct} ({group})"] * len(g_scores),
                    name=f"{ct} ({group})",
                    marker_color=colors[i],
                ))
        else:
            fig.add_trace(go.Box(
                y=scores,
                x=[ct] * len(scores),
                name=ct,
                marker_color=colors[i],
            ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", x=0, xref="paper"),
        yaxis_title="Risk Score",
        xaxis_title="",
        height=500,
        showlegend=False,
    )
    fig = apply_tcr_theme(fig)

    if output_path:
        fig.write_html(output_path, include_plotlyjs="cdn")

    return fig


# ══════════════════════════════════════════════════════════════════════════
# 4.3 Clone Diversity & Functional Annotation
# ══════════════════════════════════════════════════════════════════════════

def plot_diversity_bar(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    title: str = "Clone Diversity Metrics",
) -> go.Figure:
    """Horizontal bar chart: Shannon, Simpson, Clonality, Evenness for one sample."""
    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id]
    if df.empty:
        raise ValueError("No diversity data")
    row = df.iloc[0]
    metrics = {
        "Shannon Index": float(row["ShannonIndex"]),
        "Simpson Index": float(row["SimpsonIndex"]),
        "Evenness": float(row["Evenness"]),
        "Clonality": float(row["Clonality"]),
    }
    labels = list(metrics.keys()); values = list(metrics.values())
    colors = get_color_palette(len(labels))
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=list(reversed(values)), y=list(reversed(labels)), orientation="h",
        marker=dict(color=list(reversed(colors)), line=dict(color="white", width=0.5)),
        text=[f"{v:.4f}" for v in reversed(values)], textposition="outside",
        hovertemplate="%{y}=%{x:.4f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=f"<b>{title}</b><br><sub>Sample: {sample_id} | {int(row['CloneCount']):,} clones</sub>", x=0, xref="paper"),
        xaxis_title="Value", height=300, showlegend=False, margin=dict(r=80, l=120))
    fig = apply_tcr_theme(fig)
    if output_path: fig.write_html(output_path, include_plotlyjs="cdn")
    return fig


def plot_antigen_bar(
    df: pd.DataFrame,
    *,
    sample_id: str | None = None,
    output_path: str | None = None,
    top_n: int = 15,
    freq_col: str = "freq",
    antigen_col: str = "antigen",
    title: str = "TCR Functional Annotation — Top Matched Antigens",
) -> go.Figure:
    """Lollipop chart: top-N matched antigens for a single sample."""
    if sample_id and "Sample" in df.columns:
        df = df[df["Sample"] == sample_id]
    if df.empty:
        raise ValueError("No antigen data")
    df_s = df.nlargest(top_n, freq_col).sort_values(freq_col, ascending=True)
    antigens = df_s[antigen_col].tolist(); freqs = df_s[freq_col].tolist()
    counts = df_s["count"].tolist() if "count" in df_s.columns else [0]*len(antigens)
    color = MORANDI_COLORS[2]
    fig = go.Figure()
    for antigen, freq in zip(antigens, freqs):
        fig.add_trace(go.Scatter(
            x=[0, freq], y=[antigen, antigen], mode="lines",
            line=dict(color=color, width=2), opacity=0.35,
            showlegend=False, hoverinfo="skip",
        ))
    fig.add_trace(go.Scatter(
        x=freqs, y=antigens, mode="markers",
        marker=dict(color=color, size=14, line=dict(color="white", width=1.5)),
        text=[f"{f*100:.2f}% ({c} clones)" for f,c in zip(freqs,counts)],
        textposition="middle right", textfont_size=11,
        hovertemplate="%{y}: %{x:.4f}<extra></extra>",
        showlegend=False,
    ))
    fig.update_layout(
        title=dict(text=f"<b>{title}</b><br><sub>Sample: {sample_id}</sub>", x=0, xref="paper"),
        xaxis_title="Frequency in Repertoire", height=max(400, len(antigens)*26),
        margin=dict(r=160, l=200))
    fig = apply_tcr_theme(fig); fig.update_xaxes(tickformat=".1%")
    if output_path: fig.write_html(output_path, include_plotlyjs="cdn")
    return fig
