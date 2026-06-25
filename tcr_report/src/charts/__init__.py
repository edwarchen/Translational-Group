"""TCR report charting module.

All chart functions follow a unified signature:

    def plot_xxx(df: pd.DataFrame, ..., output_path: str | None = None) -> go.Figure:
        ...

They accept a DataFrame with expected columns, return a Plotly Figure,
and optionally save an HTML snippet to disk.
"""

from src.charts.utils import (
    MORANDI_COLORS,
    PLOTLY_THEME,
    get_color_palette,
    hex_to_rgba,
    apply_tcr_theme,
)

__all__ = [
    "MORANDI_COLORS",
    "PLOTLY_THEME",
    "get_color_palette",
    "hex_to_rgba",
    "apply_tcr_theme",
]
