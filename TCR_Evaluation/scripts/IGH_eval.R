# IGH evaluation — thin wrapper sourcing shared core.R
source("R/core.R")

cfg <- list(
  IGHV = list(
    ref_path = "data/IGHV.fasta",
    primer_path = "data/primer_set/IGHV_primers.fasta",
    direction = "fw",
    primer_id_suffix = "_fw",
    trim = 0
  ),
  IGHJ = list(
    ref_path = "data/IGHJ.fasta",
    primer_path = "data/primer_set/IGHJ_primers.fasta",
    direction = "rev",
    primer_id_suffix = "_rev",
    trim = 5
  )
)

run_full_pipeline(cfg)
