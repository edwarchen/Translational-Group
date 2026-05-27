# IGL evaluation — thin wrapper sourcing shared core.R
source("R/core.R")

cfg <- list(
  IGLV = list(
    ref_path = "data/IGLV.fasta",
    primer_path = "data/primer_set/IGLV_primers_New.fasta",
    direction = "fw",
    primer_id_suffix = "_fw",
    trim = 0
  ),
  IGLJ = list(
    ref_path = "data/IGLJ.fasta",
    primer_path = "data/primer_set/IGLJ_primers.fasta",
    direction = "rev",
    primer_id_suffix = "_rev",
    trim = 5
  )
)

run_full_pipeline(cfg)
