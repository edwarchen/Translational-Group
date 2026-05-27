library(Biostrings)
library(ggplot2)

# Parse optional CLI args: --binding= --ref= --out= --max_plots=
parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(defaults)

  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) next
    key <- kv[1]
    val <- kv[2]
    if (key %in% names(defaults)) defaults[[key]] <- val
  }
  defaults
}

# Convert sequence string into long-format table for ggplot.
seq_to_long <- function(label, seq_string) {
  chars <- strsplit(seq_string, "", fixed = TRUE)[[1]]
  data.frame(
    Track = label,
    Position = seq_along(chars),
    Base = chars,
    stringsAsFactors = FALSE
  )
}

# Make safe filename component.
safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

opts <- parse_args(list(
  binding = "results/coverage_tables/IGLV_binding_sites.csv",
  ref = "data/IGLV.fasta",
  out = "results/primer_reference_heatmaps/IGLV",
  max_plots = "Inf"
))

binding_file <- opts$binding
reference_fasta <- opts$ref
output_dir <- opts$out
max_plots <- suppressWarnings(as.numeric(opts$max_plots))
if (is.na(max_plots)) max_plots <- Inf

if (!file.exists(binding_file)) stop("binding CSV not found: ", binding_file)
if (!file.exists(reference_fasta)) stop("reference FASTA not found: ", reference_fasta)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

binding_df <- read.csv(binding_file, stringsAsFactors = FALSE)
if (nrow(binding_df) == 0) stop("No rows in binding CSV: ", binding_file)

required_cols <- c("template_id", "primer_id", "primer_sequence_match", "start", "end")

missing_cols <- setdiff(required_cols, colnames(binding_df))
if (length(missing_cols) > 0) {
  stop("binding CSV missing columns: ", paste(missing_cols, collapse = ", "))
}

ref_set <- readDNAStringSet(reference_fasta)
ref_headers <- names(ref_set)
ref_map <- list()
for (idx in seq_along(ref_headers)) {
  seq_i <- as.character(ref_set[[idx]])
  tokens <- strsplit(ref_headers[idx], "\\|")[[1]]
  # 支持用 ACCESSION 或等位基因 ID（如 IGHJ4*01）来匹配模板。
  key_candidates <- unique(tokens[seq_len(min(2, length(tokens)))])
  for (k in key_candidates) {
    if (nzchar(k)) ref_map[[k]] <- seq_i
  }
}

base_colors <- c(
  "A" = "#477C54",
  "T" = "#FF7078",
  "C" = "#5465AB",
  "G" = "#F6C971",
  "N" = "#999999",
  "-" = "#F5F5F5"
)

plot_count <- 0
primer_ids <- unique(binding_df$primer_id)

for (primer_id in primer_ids) {
  if (plot_count >= max_plots) break

  primer_df <- binding_df[binding_df$primer_id == primer_id, , drop = FALSE]
  if (nrow(primer_df) == 0) next

  # 统一以该 primer 的匹配起点对齐（同一列开始）。
  valid_idx <- which(
    !is.na(primer_df$start) &
      !is.na(primer_df$end) &
      as.integer(primer_df$start) >= 1 &
      as.integer(primer_df$end) >= as.integer(primer_df$start)
  )
  if (length(valid_idx) == 0) next
  primer_df <- primer_df[valid_idx, , drop = FALSE]

  primer_seq <- toupper(as.character(primer_df$primer_sequence_match[1]))
  anchor_start <- max(as.integer(primer_df$start))

  ref_rows <- list()
  ref_labels <- character()
  row_i <- 1

  max_len <- anchor_start - 1 + nchar(primer_seq)
  for (i in seq_len(nrow(primer_df))) {
    template_id <- primer_df$template_id[i]
    start_pos <- as.integer(primer_df$start[i])
    end_pos <- as.integer(primer_df$end[i])

    if (!template_id %in% names(ref_map)) {
      warning("template_id not found in FASTA, skip: ", template_id)
      next
    }

    ref_seq <- toupper(ref_map[[template_id]])
    ref_len <- nchar(ref_seq)
    if (end_pos > ref_len) next

    left_shift <- anchor_start - start_pos
    aligned_ref <- paste0(strrep("-", left_shift), ref_seq)
    ref_rows[[row_i]] <- aligned_ref
    ref_labels[row_i] <- paste0(template_id, " (", start_pos, "-", end_pos, ")")
    row_i <- row_i + 1
    max_len <- max(max_len, nchar(aligned_ref))
  }

  if (length(ref_rows) == 0) next

  # 末行只放一条 primer，并与 anchor_start 对齐。
  primer_row <- paste0(strrep("-", anchor_start - 1), primer_seq)
  max_len <- max(max_len, nchar(primer_row))

  rows_long <- list()
  for (i in seq_along(ref_rows)) {
    ref_rows[[i]] <- paste0(ref_rows[[i]], strrep("-", max_len - nchar(ref_rows[[i]])))
    rows_long[[i]] <- seq_to_long(ref_labels[i], ref_rows[[i]])
  }
  primer_row <- paste0(primer_row, strrep("-", max_len - nchar(primer_row)))
  primer_label <- paste0("Primer ", primer_id)
  rows_long[[length(rows_long) + 1]] <- seq_to_long(primer_label, primer_row)

  long_df <- do.call(rbind, rows_long)
  # 参考序列在上方，primer 在最下方。
  long_df$Track <- factor(long_df$Track, levels = c(primer_label, rev(ref_labels)))
  long_df$Base[!long_df$Base %in% names(base_colors)] <- "N"

  p <- ggplot(long_df, aes(x = Position, y = Track, fill = Base)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = Base), size = 2.5, family = "mono") +
    scale_fill_manual(values = base_colors) +
    scale_x_continuous(expand = c(0, 0), position = "top") +
    labs(
      title = paste0(primer_id, " aligned on binding region (n=", length(ref_rows), ")"),
      x = "Position",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(face = "bold", size = 8),
      axis.text.x = element_text(size = 8),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11)
    )

  out_name <- paste0(safe_name(primer_id), "_covered_region.png")
  out_file <- file.path(output_dir, out_name)
  ggsave(
    out_file,
    plot = p,
    width = max(9, max_len * 0.12),
    height = max(3, (length(ref_rows) + 1) * 0.35),
    dpi = 300,
    limitsize = FALSE
  )

  plot_count <- plot_count + 1
}

cat("Saved", plot_count, "primer-level plots to:", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
