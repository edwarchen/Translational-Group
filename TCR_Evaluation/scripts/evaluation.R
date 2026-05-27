# Generic evaluation script — takes target locus from CLI argument.
# Usage: Rscript scripts/evaluation.R IGH
source("R/core.R")

target <- toupper(commandArgs(trailingOnly = TRUE)[1])
if (is.na(target) || !nzchar(target)) target <- "IGH"

target_v <- paste0(target, "V")
target_j <- paste0(target, "J")

cfg <- setNames(
  list(
    list(
      ref_path = paste0("data/", target_v, ".fasta"),
      primer_path = paste0("data/primer_set/", target_v, "_primers.fasta"),
      direction = "fw",
      primer_id_suffix = "_fw",
      trim = 0
    ),
    list(
      ref_path = paste0("data/", target_j, ".fasta"),
      primer_path = paste0("data/primer_set/", target_j, "_primers.fasta"),
      direction = "rev",
      primer_id_suffix = "_rev",
      trim = 5
    )
  ),
  c(target_v, target_j)
)

run_full_pipeline(cfg)
