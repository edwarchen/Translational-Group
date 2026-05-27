# TCR/BCR Primer Evaluation Pipeline 🧬

![License](https://img.shields.io/badge/License-MIT-blue.svg)
![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![openPrimeR](https://img.shields.io/badge/openPrimeR-Supported-green)

A comprehensive bioinformatics toolkit for evaluating and visualizing the coverage of T-cell receptor (TCR) and B-cell receptor (BCR) multiplex PCR primers.

Provides both a **command-line interface** (`tcreval`) and individual R/Python scripts for flexible integration into immune repertoire sequencing (Rep-Seq) workflows.

---

## Key Features

- **One-Command Pipeline**: `tcreval run -i primers.xlsx -t TRB` runs the entire workflow end-to-end.
- **In Silico PCR Evaluation**: Uses `openPrimeR` to simulate primer binding against reference templates with up to 3 mismatches (configurable via `config/defaults.yaml`).
- **High-Resolution Visualization**:
  - **Primer Binding Heatmaps**: Nucleotide-level alignment of each primer against its reference templates.
  - **Uncovered Sequence Analysis**: Plots reference sequences missed by the primer set, grouped by gene family, to guide redesign.
- **Multi-Target Support**: Built-in reference data for `TRB`, `IGH`, `IGK`, and `IGL` loci.
- **Modular Design**: Run individual steps (`split`, `eval`, `heatmap`, `uncovered`) or the full pipeline.

---

## Quick Start

### 1. Install R Dependencies

```bash
Rscript R/install_deps.R
```

Or manually:
```R
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("openPrimeR", "Biostrings"))
install.packages(c("ggplot2", "tidyr"))
```

### 2. Install Python Package

```bash
pip install -e .          # editable install from local clone
# or
pip install git+https://github.com/edwarchen/TCR_Evaluation.git
```

Python dependencies: `click>=8.0`, `pyyaml>=6.0`, `pandas>=1.3`, `openpyxl>=3.0`.

### 3. Run

```bash
tcreval run -i data/primer_set/tcr引物序列.xlsx -t TRB
```

This executes all four steps and writes results to `results/`:

| Step | What it does |
|------|-------------|
| Split | Parses primer Excel → V-forward + J-reverse FASTA files |
| Evaluate | Runs `openPrimeR` coverage analysis → binding sites CSV, coverage stats, uncovered templates |
| Heatmap | Generates per-primer nucleotide alignment PNGs |
| Uncovered | Generates heatmaps for missed templates |

---

## CLI Reference

```
tcreval run       -i <input.xlsx> -t <TARGET>     Full pipeline
tcreval split     -i <input.xlsx> -t <TARGET>     Split primers only
tcreval eval      -t <TARGET>                      Coverage evaluation only
tcreval heatmap   -b <bindings.csv> -r <ref.fasta> -o <out/>   Binding heatmaps
tcreval uncovered -d <uncovered_dir/>              Uncovered template plots
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i`, `--input` | Primer Excel/CSV file (first two columns: ID, Sequence) | Required |
| `-t`, `--target` | Target locus: `TRB`, `IGH`, `IGK`, or `IGL` | Required |
| `-o`, `--output-dir` | Results directory | `results/` |
| `--max-mismatch` | Maximum allowed mismatches | `3` |
| `--binding-ratio` | Allowed off-target binding ratio | `1.0` |
| `--v-pattern` | Substring to identify V primers in ID | `V` |
| `--j-pattern` | Substring to identify J primers in ID | `J` |

### Input File Format

Excel (`.xlsx`) or CSV with at least two columns:
1. **Primer ID** — e.g. `TRBV1`, `IGHJ4`
2. **Sequence** — e.g. `GCTACTTCGGAGCCTCGG`

Primers are classified as V (forward) or J (reverse) by the presence of `V` or `J` in their ID.

---

## Project Structure

```
├── tcreval/                       # Python CLI package
│   ├── cli.py                     # Click-based command-line interface
│   ├── splitter.py                # Primer Excel → FASTA splitting
│   ├── orchestrator.py            # R subprocess scheduling + full pipeline
│   ├── config.py                  # YAML config loader
│   └── utils.py                   # Logging, path resolution, R preflight
├── R/
│   ├── core.R                     # Shared R functions (coverage, binding, plotting)
│   └── install_deps.R             # One-shot R dependency installer
├── config/
│   └── defaults.yaml              # Max mismatches, target paths, trim settings
├── scripts/                       # Standalone R/Python scripts (also callable directly)
│   ├── split_primers.py           # Parameterized primer splitter
│   ├── evaluation.R               # Generic evaluator (Rscript evaluation.R TRB)
│   ├── IGH_eval.R                 # IGH-specific thin wrapper
│   ├── IGK_eval.R                 # IGK-specific thin wrapper
│   ├── IGL_eval.R                 # IGL-specific thin wrapper
│   ├── TRB_eval.R                 # TRB-specific thin wrapper
│   ├── Primer_coverage_heatmap.R  # Per-primer alignment heatmaps
│   ├── Uncovered_plot.R           # Uncovered template visualization
│   └── official_pipeline.R        # openPrimeR vignette reference
├── data/                          # Reference FASTAs + primer sets
│   ├── {IGH,IGK,IGL,TRB}{V,J}.fasta
│   └── primer_set/
├── pyproject.toml                 # pip installation metadata
├── requirements.txt               # Python dependencies
└── results/                       # Generated output (gitignored)
    ├── coverage_tables/           # Binding sites + coverage stats (CSV)
    ├── coverage_plots/            # Summary coverage heatmaps (PNG)
    ├── primer_reference_heatmaps/ # Per-primer nucleotide alignments (PNG)
    └── uncovered_templates/       # Uncovered sequences + heatmaps by gene family
```

---

## Configuration

Edit `config/defaults.yaml` to change:

```yaml
evaluation:
  max_mismatch: 3              # PCR mismatch tolerance
  allowed_binding_ratio: 1.0   # Off-target binding allowance
  trim_fw: 0                   # Bases to trim from forward primers
  trim_rev: 5                  # Bases to trim from reverse primers

targets:
  TRB: ...
  IGH: ...
  IGK: ...                     # Easily add custom targets
  IGL: ...
```

---

## Use Cases

1. **Primer Design Validation**: Evaluate a new multiplex primer panel before ordering — identify coverage gaps at nucleotide resolution.
2. **Blind Spot Detection**: Review `uncovered_templates/` heatmaps to see which conserved alleles are missed.
3. **Mismatch Optimization**: Check `primer_reference_heatmaps/` — if a mismatch is consistently at the 3' end, the primer should be redesigned.
4. **IP Documentation**: The binding site CSVs and comprehensive heatmaps serve as technical disclosure materials for patent applications.

---

## License

MIT License.

---

*Developed for advanced Immune Repertoire Sequencing (Rep-Seq) primer design and evaluation.*
