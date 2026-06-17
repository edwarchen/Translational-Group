"""Common utilities for TCR report Plotly charts.

Provides Morandi color palettes, Plotly theme/layout presets,
and helper functions shared across all chart modules.
"""

from __future__ import annotations

import plotly.graph_objects as go
import plotly.io as pio
import numpy as np


# ── Morandi Color Palette (20 distinct colors) ──────────────────────────
# Derived from seaborn muted + pastel palettes, suitable for
# distinguishing amino acids, V genes, or sample groups.

MORANDI_COLORS = [
    "#7FB3D5",  # muted blue
    "#F1948A",  # muted red
    "#85C1E9",  # pastel blue
    "#F5B7B1",  # pastel red
    "#82E0AA",  # muted green
    "#D7BDE2",  # muted purple
    "#73C6B6",  # pastel green
    "#E8DAEF",  # pastel purple
    "#F8C471",  # muted orange
    "#AED6F1",  # pastel cyan
    "#F0B27A",  # pastel orange
    "#BB8FCE",  # muted violet
    "#A3E4D7",  # pastel teal
    "#DC7633",  # muted brown
    "#F9E79F",  # pastel yellow
    "#AAB7B8",  # muted gray
    "#D5DBDF",  # pastel silver
    "#B03A2E",  # muted dark red
    "#1A5276",  # muted dark blue
    "#117A65",  # muted dark green
]

# Grey for "Other" / background
OTHER_GREY = "#E0E0E0"


# ── Plotly Theme ────────────────────────────────────────────────────────

PLOTLY_THEME: dict = dict(
    font_family="Arial, sans-serif",
    font_size=12,
    font_color="#333333",
    title_font_size=16,
    legend_font_size=11,
    plot_bgcolor="white",
    paper_bgcolor="white",
    margin=dict(l=60, r=40, t=60, b=60),
    hoverlabel=dict(
        bgcolor="white",
        font_size=12,
        font_family="Arial, sans-serif",
    ),
    colorway=MORANDI_COLORS,
)


def apply_tcr_theme(fig: go.Figure) -> go.Figure:
    """Apply the default TCR report theme to a Plotly figure."""
    fig.update_layout(**PLOTLY_THEME)
    fig.update_xaxes(
        showgrid=True,
        gridcolor="#E9ECEF",
        gridwidth=0.8,
        zeroline=False,
        linecolor="#CCCCCC",
        linewidth=1,
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="#E9ECEF",
        gridwidth=0.8,
        zeroline=False,
        linecolor="#CCCCCC",
        linewidth=1,
    )
    return fig


def get_color_palette(n: int) -> list[str]:
    """Return ``n`` colors from the Morandi palette, cycling if needed."""
    if n <= len(MORANDI_COLORS):
        return MORANDI_COLORS[:n]
    # Cycle through palette
    repeats = (n + len(MORANDI_COLORS) - 1) // len(MORANDI_COLORS)
    palette = (MORANDI_COLORS * repeats)[:n]
    return palette


def hex_to_rgba(hex_color: str, alpha: float = 0.4) -> str:
    """Convert hex color to rgba string with given alpha."""
    hex_color = hex_color.lstrip("#")
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    return f"rgba({r}, {g}, {b}, {alpha})"


def save_figure(fig: go.Figure, output_path: str) -> None:
    """Save a Plotly figure as a standalone HTML file.

    Uses ``include_plotlyjs='cdn'`` for report embedding — the
    report template loads plotly.js once from CDN.
    """
    fig.write_html(
        output_path,
        full_html=True,
        include_plotlyjs="cdn",
    )


def figure_to_html_div(fig: go.Figure, div_id: str | None = None) -> str:
    """Convert a Plotly figure to an HTML div+script snippet for embedding.

    Parameters
    ----------
    fig : go.Figure
        The Plotly figure to convert.
    div_id : str or None
        HTML div id. Auto-generated if None.

    Returns
    -------
    str
        HTML snippet containing <div> and <script> tags.
    """
    if div_id is None:
        import uuid
        div_id = f"plotly-{uuid.uuid4().hex[:8]}"

    return fig.to_html(
        full_html=False,
        include_plotlyjs=False,
        div_id=div_id,
    )
