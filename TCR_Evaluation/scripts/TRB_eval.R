# TRB evaluation — thin wrapper sourcing shared core.R
source("R/core.R")

cfg <- list(
  TRBV = list(
    ref_path = "data/TRBV.fasta",
    primer_path = "data/primer_set/TRBV_primers.fasta",
    direction = "fw",
    primer_id_suffix = "_fw",
    trim = 0
  ),
  TRBJ = list(
    ref_path = "data/TRBJ.fasta",
    primer_path = "data/primer_set/TRBJ_primers.fasta",
    direction = "rev",
    primer_id_suffix = "_rev",
    trim = 5
  )
)

run_full_pipeline(cfg)
