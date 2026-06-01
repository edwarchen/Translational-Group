# PISTE TEST Results

External evaluation of PISTE (Physics Inspired Sliding TransformEr) on independent VDJdb test sets.

## Evaluation Setup

- **Data source**: VDJdb filtered triples (TCR CDR3, peptide, HLA)
- **Negative sampling**: Random shuffle, 10:1 neg:pos ratio
- **Metrics**: AUC, AUPR, PPV@N (precision among top-N predictions)
- **Seeds**: Fixed random seed 42 for reproducibility

## Results

### Dataset: VDJdb score > 0 (1,480 pos / 14,800 neg)

| Model | AUC | AUPR | PPV@N |
|-------|-----|------|-------|
| PISTE-random | 0.8927 | 0.5052 | 0.5318 |
| PISTE-unipep | 0.8609 | 0.5807 | 0.5399 |
| PISTE-reftcr | 0.5614 | 0.1258 | 0.1115 |

### Dataset: VDJdb score >= 2 (607 pos / 6,070 neg)

| Model | AUC | AUPR | PPV@N |
|-------|-----|------|-------|
| PISTE-random | 0.8278 | 0.4243 | 0.4135 |
| PISTE-unipep | 0.7966 | 0.4674 | 0.4069 |
| PISTE-reftcr | 0.5406 | 0.1231 | 0.1318 |

## Key Findings

- **PISTE-random** achieves the highest AUC (0.89) on external VDJdb data
- **PISTE-unipep** has the best AUPR (0.58), indicating better precision-recall trade-off under class imbalance
- **PISTE-reftcr** performs near random on external VDJdb (AUC ~0.54-0.56), suggesting a generalization gap vs its strong internal test performance (AUC 0.95 on internal test set)

## Files

- `summary.csv` — Summary metrics table
- `piste_evaluation_curves.pdf/png` — ROC and PR curves (Nature-style 2x2 panel)
- `*_preds.npz` — Raw prediction scores (y_true, y_prob) for each dataset/model
- `run_piste_eval.py` — Evaluation pipeline script
- `evaluate_piste.py` — Alternative evaluation script with multiple negative sampling methods
- `plot_piste_eval.py` — Plotting script for Nature-style figures
