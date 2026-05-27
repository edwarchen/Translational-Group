# IGK evaluation — thin wrapper sourcing shared core.R
source("R/core.R")

cfg <- list(
  IGKV = list(
    ref_path = "data/IGKV.fasta",
    primer_path = "data/primer_set/IGKV_primers.fasta",
    direction = "fw",
    primer_id_suffix = "_fw",
    trim = 0
  ),
  IGKJ = list(
    ref_path = "data/IGKJ.fasta",
    primer_path = "data/primer_set/IGKJ_primers.fasta",
    direction = "rev",
    primer_id_suffix = "_rev",
    trim = 5
  )
)

run_full_pipeline(cfg)
