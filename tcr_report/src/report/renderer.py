"""Jinja2 HTML report renderer.

Renders the base template with chart HTML snippets, data tables,
and report metadata into a standalone HTML file.
"""

from __future__ import annotations

import base64
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape


class ReportRenderer:
    """Jinja2-based HTML report renderer.

    Usage::

        renderer = ReportRenderer(template_dir="templates", logo_path="海普洛斯-蓝橙.png")
        html = renderer.render(
            report_title="TCR 免疫组库测序报告",
            company="My Lab",
            report_date="2025-09-18",
            sample_info_table="<table>...</table>",
            qc_table="<table>...</table>",
            single_sample_charts=[{"title": "...", "html": "...", "caption": "..."}],
            multi_sample_charts=[...],
        )
        renderer.write(html, "output/report.html")
    """

    def __init__(
        self,
        template_dir: str | Path = "templates",
        logo_path: str | None = None,
    ):
        self.template_dir = Path(template_dir)
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=select_autoescape(["html", "xml"]),
        )

        # Pre-load logo as base64
        self.logo_base64 = ""
        if logo_path:
            logo = Path(logo_path)
            if logo.exists():
                b64 = base64.b64encode(logo.read_bytes()).decode("utf-8")
                ext = logo.suffix.lstrip(".").lower()
                self.logo_base64 = f"data:image/{ext};base64,{b64}"

    def render(
        self,
        report_title: str = "TCR 免疫组库测序报告",
        company: str = "实验室",
        report_date: str | None = None,
        report_year: str | None = None,
        sample_info_table: str = "",
        qc_table: str = "",
        single_sample_charts: dict[str, list[dict]] | None = None,
        multi_sample_charts: list[dict] | None = None,
        version: str = "1.0.0",
        **extra_context,
    ) -> str:
        """Render the full HTML report."""
        import datetime

        if report_date is None or report_date == "auto":
            report_date = datetime.date.today().strftime("%Y-%m-%d")

        if report_year is None:
            report_year = str(datetime.date.today().year)

        if single_sample_charts is None:
            single_sample_charts = {"sec41": [], "sec42": [], "sec43": []}
        if multi_sample_charts is None:
            multi_sample_charts = []

        template = self.env.get_template("base.html")

        return template.render(
            report_title=report_title,
            company=company,
            report_date=report_date,
            report_year=report_year,
            logo_base64=self.logo_base64,
            sample_info_table=sample_info_table,
            qc_table=qc_table,
            single_sample_charts=single_sample_charts,
            multi_sample_charts=multi_sample_charts,
            version=version,
            **extra_context,
        )

    def write(self, html: str, output_path: str | Path) -> None:
        """Write rendered HTML to a file."""
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(html, encoding="utf-8")
