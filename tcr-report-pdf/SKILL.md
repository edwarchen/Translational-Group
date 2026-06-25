---
name: tcr-report-pdf
description: Convert standard-version Haplox TCR immune repertoire HTML reports into print-ready PDF reports. Use when the user already has one or more standardized Haplox/TCR report HTML files and asks to generate, regenerate, export, convert, or verify PDF reports without header overlap, chart/table truncation, or inconsistent font hierarchy.
---

# TCR Report PDF

## Overview

Use this skill when the user already has standard-version Haplox TCR immune repertoire HTML reports and needs matching standardized A4 PDFs.

This skill does not create the report content or repair arbitrary HTML templates. Its job is the final HTML-to-PDF production step: apply the established print CSS, generate PDF files, render previews, and verify that the output is visually deliverable.

The current standard intentionally does not use Playwright PDF header templates. Earlier repeated-header attempts caused report headers to overlap body content. Keep the HTML report header as normal page content and inject print CSS instead.

## Input Assumptions

Proceed when the input is a standard Haplox TCR report HTML containing the expected report structure, including selectors such as `.report-header`, `.cover-page`, `.section`, `.figure-block`, `.table-block`, and `#sec-background .sop-workflow`.

If the user provides a non-standard HTML report, inspect the template first and adjust the print CSS only after confirming the structural differences. Do not promise the standard PDF output for unrelated HTML templates.

## Standard Workflow

1. Confirm the input `.html` files are standard-version Haplox TCR reports.
2. Identify the target output directory. By default, write each `.pdf` next to its source `.html`.
3. Generate PDFs with `scripts/generate_tcr_report_pdf.py`.
4. Render the generated PDFs to page images or contact sheets.
5. Inspect the rendered output visually before saying the PDFs are complete.
6. Report the PDF paths and the preview paths to the user.

## Generation

Run the bundled converter from the skill directory:

```bash
python3 scripts/generate_tcr_report_pdf.py /path/to/report.html
```

For several reports:

```bash
python3 scripts/generate_tcr_report_pdf.py \
  /path/to/S014_report.html \
  /path/to/S015_report.html \
  /path/to/S016_report.html
```

To write PDFs into a specific directory:

```bash
python3 scripts/generate_tcr_report_pdf.py --out-dir /path/to/pdf_dir /path/to/*_report.html
```

The converter uses Playwright Chromium, A4 paper, printed backgrounds, and margins of `15mm` top/bottom and `10mm` left/right.

## Print Standard

Maintain these typography and layout rules unless the user explicitly asks for a change:

- Body text: `15px`, line-height `1.55`.
- Section title `.section h2`: `22px`.
- Subsection title `.section h3`: `18px`.
- Figure and table titles: `16px`.
- Notes: `13px`.
- Table cells: `12px`.
- Keep the 1.3 workflow diagram on one page and prevent it from being clipped.
- Allow normal sections to break across pages to avoid blank heading-only pages.
- Keep figure and table blocks together where possible.
- Avoid fixed/repeated page headers that overlap content.

## Verification

Never finish with only the generated PDF path. Render and inspect the PDF first.

Recommended checks:

- Page count is plausible and no page is blank except intentional spacing.
- Report header does not overlap titles or body text.
- `1.3 信息分析流程图` is complete and not cut off.
- Major section headings are larger than table/figure titles.
- Body text is larger than notes.
- Tables fit horizontally; cell text wraps instead of being clipped.
- Plots, Sankey diagrams, pie charts, and legends are visible.
- Final page footer/reference/software sections remain readable.

If rendering tools are already present in the workspace, reuse them. In this project, `work/render_pdf_pages.swift` can render each PDF page to PNG, and an existing contact-sheet helper may be used to inspect all pages together.

## Delivery

In the final response, include:

- Each generated PDF path as a clickable file link.
- The preview/contact-sheet path if one was created.
- A concise statement of what was checked, especially header overlap, workflow truncation, and font hierarchy.

If any issue remains, state it directly and keep the PDF as a draft rather than calling it final.
