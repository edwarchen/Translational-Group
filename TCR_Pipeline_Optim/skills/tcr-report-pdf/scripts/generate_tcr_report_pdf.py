#!/usr/bin/env python3
"""Convert Haplox TCR immune repertoire HTML reports to standardized PDFs."""

from __future__ import annotations

import argparse
import time
from pathlib import Path


PRINT_CSS = """
@media print {
    .report-header {
        position: static !important;
        max-width: 100% !important;
        padding: 12px 28px !important;
        margin: 0 auto 14px !important;
    }

    body {
        font-size: 15px !important;
        line-height: 1.55 !important;
    }

    .cover-page {
        page-break-after: always !important;
        break-after: page !important;
    }

    .page {
        max-width: 100% !important;
        padding: 0 16px !important;
    }

    .section {
        break-inside: auto !important;
        page-break-inside: auto !important;
    }

    .figure-block, .table-block {
        break-inside: avoid !important;
        page-break-inside: avoid !important;
    }

    .workflow-heading {
        break-before: page !important;
        page-break-before: always !important;
        break-after: avoid !important;
        page-break-after: avoid !important;
    }

    #sec-background .sop-workflow {
        break-before: auto !important;
        page-break-before: auto !important;
        break-inside: avoid !important;
        page-break-inside: avoid !important;
        margin: 6px 0 8px !important;
        padding: 10px 12px !important;
    }

    .sop-workflow-title {
        font-size: 16px !important;
        margin-bottom: 6px !important;
    }

    .sop-input-block {
        padding: 8px 10px 9px !important;
    }

    .sop-block-title {
        font-size: 14px !important;
        margin-bottom: 6px !important;
        padding-bottom: 4px !important;
    }

    .sop-flow-row {
        gap: 5px !important;
    }

    .sop-step {
        min-height: 58px !important;
        padding: 6px 6px !important;
    }

    .sop-step .sop-index {
        width: 20px !important;
        height: 20px !important;
        margin-bottom: 3px !important;
        font-size: 11px !important;
    }

    .sop-step .sop-name {
        font-size: 13px !important;
        line-height: 1.2 !important;
    }

    .sop-step .sop-desc {
        font-size: 10px !important;
        line-height: 1.2 !important;
        margin-top: 3px !important;
    }

    .sop-arrow {
        flex-basis: 14px !important;
        font-size: 14px !important;
    }

    .sop-analysis-link,
    .sop-dashed-link {
        height: 18px !important;
    }

    .sop-analysis-link::before,
    .sop-dashed-link::before {
        bottom: 5px !important;
    }

    .sop-analysis-link::after,
    .sop-dashed-link::after {
        bottom: 1px !important;
        left: calc(50% - 4px) !important;
        border-left-width: 4px !important;
        border-right-width: 4px !important;
        border-top-width: 5px !important;
    }

    .sop-analysis-section {
        padding: 8px 10px 9px !important;
    }

    .sop-analysis-section h4 {
        font-size: 14px !important;
        margin-bottom: 6px !important;
        padding-bottom: 4px !important;
    }

    .sop-module-row {
        gap: 7px !important;
    }

    .sop-module-card {
        min-height: 54px !important;
        padding: 6px 7px !important;
    }

    .sop-module-title {
        font-size: 12px !important;
        line-height: 1.2 !important;
        margin-bottom: 3px !important;
    }

    .sop-module-card ul {
        font-size: 10px !important;
        line-height: 1.22 !important;
        padding-left: 12px !important;
    }

    .sop-module-card.compact {
        min-height: 32px !important;
        padding: 5px 6px !important;
    }

    .sop-workflow-footnote {
        margin-top: 6px !important;
        font-size: 10.5px !important;
        line-height: 1.3 !important;
    }

    .table-wrapper {
        overflow-x: visible !important;
        overflow-y: visible !important;
        max-height: none !important;
    }

    .section h2 {
        font-size: 22px !important;
        line-height: 1.35 !important;
        margin: 14px 0 12px !important;
        padding-bottom: 6px !important;
    }

    .section h3 {
        font-size: 18px !important;
        line-height: 1.35 !important;
        margin: 12px 0 9px !important;
    }

    .section p {
        font-size: 15px !important;
        line-height: 1.55 !important;
        margin-bottom: 10px !important;
    }

    .figure-block .figure-title,
    .table-block .table-title {
        font-size: 16px !important;
        line-height: 1.35 !important;
        margin-bottom: 9px !important;
    }

    .table-note,
    .table-note p {
        font-size: 13px !important;
        line-height: 1.45 !important;
    }

    .table-note {
        margin-bottom: 10px !important;
    }

    .figure-note {
        font-size: 13px !important;
        line-height: 1.45 !important;
        padding: 6px 10px !important;
    }

    table.data-table {
        font-size: 12px !important;
    }

    table.data-table tbody td, table.data-table thead th {
        white-space: normal !important;
        word-break: break-word !important;
        padding: 4px 6px !important;
    }
}
"""


def default_pdf_path(html_path: Path, out_dir: Path | None) -> Path:
    if out_dir is None:
        return html_path.with_suffix(".pdf")
    return out_dir / f"{html_path.stem}.pdf"


def html_to_pdf(html_path: Path, pdf_path: Path, wait_seconds: float) -> None:
    from playwright.sync_api import sync_playwright

    html_path = html_path.resolve()
    pdf_path = pdf_path.resolve()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with sync_playwright() as p:
        browser = p.chromium.launch(
            args=[
                "--enable-unsafe-swiftshader",
                "--use-gl=swiftshader",
                "--ignore-gpu-blocklist",
            ]
        )
        page = browser.new_page()
        page.goto(f"file://{html_path}", wait_until="networkidle")
        time.sleep(wait_seconds)
        page.add_style_tag(content=PRINT_CSS)
        page.evaluate(
            """
            () => {
                const workflow = document.querySelector('#sec-background .sop-workflow');
                if (workflow && workflow.previousElementSibling) {
                    workflow.previousElementSibling.classList.add('workflow-heading');
                }
            }
            """
        )
        page.pdf(
            path=str(pdf_path),
            format="A4",
            print_background=True,
            margin={"top": "15mm", "bottom": "15mm", "left": "10mm", "right": "10mm"},
        )
        browser.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert one or more Haplox TCR report HTML files to standardized PDFs."
    )
    parser.add_argument("html", nargs="+", type=Path, help="Input report HTML file(s).")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output PDF path. Only valid when one HTML file is provided.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory for batch conversion. Defaults to each HTML file directory.",
    )
    parser.add_argument(
        "--wait-seconds",
        type=float,
        default=2.0,
        help="Seconds to wait after network idle so charts finish rendering. Default: 2.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output and len(args.html) != 1:
        raise SystemExit("--output can only be used with a single HTML file.")

    for html_path in args.html:
        if not html_path.exists():
            raise SystemExit(f"HTML not found: {html_path}")
        pdf_path = args.output if args.output else default_pdf_path(html_path, args.out_dir)
        html_to_pdf(html_path, pdf_path, args.wait_seconds)
        size_mb = pdf_path.stat().st_size / 1024 / 1024
        print(f"PDF saved: {pdf_path.resolve()} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
