# TCR/BCR Primer Evaluation — shared R functions
# Source this file from any evaluation script: source("R/core.R")
library(openPrimeR)
library(Biostrings)

# ---- Constants ----
HDR_STRUCTURE <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
DEFAULT_OUTPUT_BASE <- "results"

# ---- Settings ----
create_default_settings <- function(max_mismatch = 3, allowed_binding_ratio = 1.0) {
  xml_file <- system.file("extdata", "settings", "A_Taq_PCR_design.xml", package = "openPrimeR")
  eval_settings <- read_settings(xml_file)
  conOptions(eval_settings)$allowed_other_binding_ratio <- c(max = allowed_binding_ratio)
  conOptions(eval_settings)$allowed_mismatches <- c(max = max_mismatch)
  eval_settings
}

# ---- Load reference templates ----
load_ref <- function(ref_path, direction) {
  ref <- read_templates(ref_path, HDR_STRUCTURE, delim = "|", id.column = "GROUP")
  ref$Group <- gsub("[-*].*", "", ref$Group)

  if (direction == "fw") {
    ref$Allowed_Start_fw <- 1
    ref$Allowed_End_fw <- nchar(ref$Sequence)
    ref$Allowed_fw <- ref$Sequence
  }
  ref
}

# ---- Load and trim primers ----
load_primers_by_direction <- function(primer_path, direction, primer_id_suffix, trim_len = 0) {
  primers <- if (direction == "fw") {
    read_primers(primer_path, fw.id = primer_id_suffix)
  } else {
    read_primers(primer_path, rev.id = primer_id_suffix)
  }

  if (trim_len > 0) {
    new_seqs <- sub(paste0("^.{", trim_len, "}"), "", primers$Reverse)
    primers$Reverse <- new_seqs
  } else if (trim_len > 0 && direction == "fw") {
    new_seqs <- sub(paste0("^.{", trim_len, "}"), "", primers$Forward)
    primers$Forward <- new_seqs
  }
  primers$Length <- nchar(ifelse(direction == "fw", primers$Forward, primers$Reverse))
  return(primers)
}

# ---- Run coverage evaluation ----
run_coverage <- function(one_cfg, eval_settings) {
  ref <- load_ref(one_cfg$ref_path, one_cfg$direction)
  primers <- load_primers_by_direction(
    one_cfg$primer_path,
    one_cfg$direction,
    one_cfg$primer_id_suffix,
    trim_len = if (is.null(one_cfg$trim)) 0 else one_cfg$trim
  )
  constraint.df <- check_constraints(
    primers,
    ref,
    eval_settings,
    active.constraints = "primer_coverage"
  )
  template.df <- update_template_cvg(ref, constraint.df)
  list(
    direction = one_cfg$direction,
    primers = primers,
    constraint.df = constraint.df,
    template.df = template.df
  )
}

# ---- Export binding site coordinates ----
export_binding_sites <- function(name, one_result, output_dir = "results/coverage_tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(output_dir, paste0(name, "_binding_sites.csv"))

  primers <- one_result$primers
  template.df <- one_result$template.df
  direction <- one_result$direction
  primer_col <- if (direction == "fw") "Forward" else "Reverse"

  clean_nt <- function(x) {
    x <- toupper(as.character(x))
    gsub("[^ACGTRYSWKMBDHVN]", "N", x)
  }

  records <- list()
  rec_i <- 1

  for (i in seq_len(nrow(primers))) {
    primer_id <- as.character(primers$ID[i])
    primer_seq_input <- clean_nt(primers[[primer_col]][i])
    if (!nzchar(primer_seq_input)) next

    primer_seq_match <- if (direction == "fw") {
      primer_seq_input
    } else {
      as.character(reverseComplement(DNAString(primer_seq_input)))
    }
    pattern <- DNAString(primer_seq_match)

    for (j in seq_len(nrow(template.df))) {
      template_id <- as.character(template.df$ID[j])
      template_group <- as.character(template.df$Group[j])
      template_seq <- clean_nt(template.df$Sequence[j])
      if (!nzchar(template_seq)) next

      hits <- matchPattern(pattern, DNAString(template_seq), max.mismatch = 3, fixed = FALSE)
      if (length(hits) == 0) next

      for (k in seq_along(hits)) {
        hit_start <- start(hits)[k]
        hit_end <- end(hits)[k]
        safe_match_seq <- substr(template_seq, max(1, hit_start), min(nchar(template_seq), hit_end))

        records[[rec_i]] <- data.frame(
          target = name,
          template_group = template_group,
          template_id = template_id,
          primer_id = primer_id,
          primer_direction = direction,
          primer_sequence_input = primer_seq_input,
          primer_sequence_match = primer_seq_match,
          start = hit_start,
          end = hit_end,
          matched_sequence = safe_match_seq,
          stringsAsFactors = FALSE
        )
        rec_i <- rec_i + 1
      }
    }
  }

  out_df <- if (length(records) > 0) do.call(rbind, records) else data.frame()
  write.csv(out_df, out_file, row.names = FALSE)
  normalizePath(out_file, winslash = "/", mustWork = FALSE)
}

# ---- Plot coverage heatmap ----
summarize_and_plot <- function(name, one_result, plot_dir = "results/coverage_plots") {
  constraint.df <- one_result$constraint.df
  template.df <- one_result$template.df
  ratio <- as.numeric(get_cvg_ratio(constraint.df, template.df))
  cvg.stats <- get_cvg_stats(constraint.df, template.df, for.viewing = TRUE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_file <- file.path(plot_dir, paste0(name, "_cov.png"))

  grDevices::png(plot_file, width = 1400, height = 900, res = 120)
  print(plot_template_cvg(constraint.df, template.df))
  grDevices::dev.off()

  list(name = name, ratio = ratio, cvg.stats = cvg.stats, plot_file = plot_file)
}

# ---- Save coverage statistics CSV ----
save_cov_stats <- function(name, one_summary, output_dir = "results/coverage_tables") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stats_file <- file.path(output_dir, paste0(name, "_cov.csv"))
  write.csv(one_summary$cvg.stats, stats_file, row.names = FALSE)
  normalizePath(stats_file, winslash = "/", mustWork = FALSE)
}

# ---- Save uncovered sequences grouped by gene family ----
save_uncovered_sequences <- function(name, one_result, base_dir = "results/uncovered_templates") {
  template.df <- one_result$template.df

  uncovered_df <- template.df[template.df$primer_coverage == 0, ]

  if (nrow(uncovered_df) == 0) {
    return("No uncovered sequences found.")
  }

  target_dir <- file.path(base_dir, paste0(name, "_unc"))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  family_list <- split(uncovered_df, uncovered_df$Group)

  for (family_name in names(family_list)) {
    sub_df <- family_list[[family_name]][, c("Group", "ID", "Sequence")]
    safe_family <- gsub("[/]", "_", family_name)
    file_path <- file.path(target_dir, paste0(safe_family, ".csv"))
    write.csv(sub_df, file_path, row.names = FALSE)
  }

  return(normalizePath(target_dir, winslash = "/", mustWork = FALSE))
}

# ---- Run full pipeline for a given config ----
run_full_pipeline <- function(cfg, target_name = NULL, eval_settings = NULL, base_dir = "results") {
  if (is.null(eval_settings)) {
    eval_settings <- create_default_settings()
  }

  results <- lapply(cfg, run_coverage, eval_settings = eval_settings)
  if (is.null(names(results))) {
    names(results) <- paste0("target_", seq_along(results))
  }

  saved_binding_sites <- mapply(export_binding_sites, names(results), results,
                                MoreArgs = list(output_dir = file.path(base_dir, "coverage_tables")),
                                SIMPLIFY = FALSE)

  coverage_summary <- mapply(summarize_and_plot, names(results), results,
                             MoreArgs = list(plot_dir = file.path(base_dir, "coverage_plots")),
                             SIMPLIFY = FALSE)

  saved_stats <- mapply(save_cov_stats, names(coverage_summary), coverage_summary,
                        MoreArgs = list(output_dir = file.path(base_dir, "coverage_tables")),
                        SIMPLIFY = FALSE)

  saved_uncovered_paths <- mapply(save_uncovered_sequences, names(results), results,
                                  MoreArgs = list(base_dir = file.path(base_dir, "uncovered_templates")),
                                  SIMPLIFY = FALSE)

  for (target_name in names(results)) {
    cat("\n===", target_name, "primer coverage ===\n")
    cat(target_name, "coverage ratio:", coverage_summary[[target_name]]$ratio, "\n")
    print(coverage_summary[[target_name]]$cvg.stats)
    cat(target_name, "coverage plot saved:", coverage_summary[[target_name]]$plot_file, "\n")
    cat(target_name, "coverage table saved:", saved_stats[[target_name]], "\n")
    cat(target_name, "binding sites saved:", saved_binding_sites[[target_name]], "\n")
    cat(target_name, "uncovered templates dir:", saved_uncovered_paths[[target_name]], "\n")
  }
}
