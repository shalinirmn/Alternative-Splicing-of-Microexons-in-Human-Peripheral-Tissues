#!/usr/bin/env Rscript

# Usage:
# Rscript scripts/04_differential_microexons.R \
#   --input vast_align/INCLUSION_LEVELS_FULL-Hsa.tab \
#   --metadata config/sample_metadata.tsv \
#   --baseline Control \
#   --case Disease \
#   --name cohort_wilcox \
#   --cutoff_baseline 3 \
#   --cutoff_case 20 \
#   --outdir results/differential_splicing

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    cutoff_baseline = 3,
    cutoff_case = NULL,
    outdir = "results/differential_splicing"
  )
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    opts[[key]] <- args[i + 1]
    i <- i + 2
  }
  if (is.null(opts$input) || is.null(opts$metadata) || is.null(opts$baseline) || is.null(opts$case) || is.null(opts$name)) {
    stop("Required arguments: --input <INCLUSION_LEVELS_FULL.tab> --metadata <meta.tsv> --baseline <grp1> --case <grp2> --name <run_label>")
  }
  opts$cutoff_baseline <- as.numeric(opts$cutoff_baseline)
  opts$cutoff_case <- if (is.null(opts$cutoff_case)) opts$cutoff_baseline else as.numeric(opts$cutoff_case)
  return(opts)
}

opts <- parse_args(args)
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

# Function 1: Identify microexons (<= 27 nt) from VAST-TOOLS diff output
microfilt <- function(total_table, diff_table, direction) {
  microexon_list <- total_table[
    (total_table$COMPLEX %in% c("S", "C1", "C2", "C3", "MIC", "ANN")) &
      (total_table$LENGTH <= 27),
    "EVENT"
  ]
  microexon_list <- microexon_list[!is.na(microexon_list)]

  diff_microexon_table <- subset(diff_table, EVENT %in% microexon_list)
  diff_microexon_table2 <- subset(diff_microexon_table, MV.dPsi._at_0.95 > 0)

  if (missing(direction)) {
    dec <- diff_microexon_table2[which(diff_microexon_table2$E.dPsi. < 0), ]
    inc <- diff_microexon_table2[which(diff_microexon_table2$E.dPsi. > 0), ]
  } else {
    dec <- diff_microexon_table2[which(diff_microexon_table2$E.dPsi. > 0), ]
    inc <- diff_microexon_table2[which(diff_microexon_table2$E.dPsi. < 0), ]
  }

  out_list <- list(
    ME_total = diff_microexon_table,
    ME_dec = dec,
    ME_inc = inc
  )
  return(out_list)
}

# Function 2: Perform Wilcoxon rank-sum test across PSI matrix and filter events
microfilt_wilcox <- function(total_table, config_rep, baseline, case, name, cutoff_baseline, cutoff_case, outdir) {
  baseline_samples <- subset(config_rep, condition == baseline)$sample
  case_samples <- subset(config_rep, condition == case)$sample

  missing_cols <- setdiff(c(baseline_samples, case_samples), colnames(total_table))
  if (length(missing_cols) > 0) {
    stop("The following samples from metadata were not found in the inclusion table: ", paste(missing_cols, collapse = ", "))
  }

  baseline_psi <- total_table[, baseline_samples, drop = FALSE]
  case_psi <- total_table[, case_samples, drop = FALSE]

  # Remove events completely missing in either group
  rm_rows1 <- which(rowSums(!is.na(case_psi)) == 0 | rowSums(!is.na(baseline_psi)) == 0)
  if (length(rm_rows1) > 0) {
    case_psi <- case_psi[-rm_rows1, , drop = FALSE]
    baseline_psi <- baseline_psi[-rm_rows1, , drop = FALSE]
    diff_table <- total_table[-rm_rows1, 1:6, drop = FALSE]
  } else {
    diff_table <- total_table[, 1:6, drop = FALSE]
  }

  # Keep events with non-zero variance and >= 2 observations in at least one group
  tmp1 <- which(apply(case_psi, 1, function(x) {
    x <- na.omit(x)
    length(x) > 1 && var(x) != 0
  }))
  tmp2 <- which(apply(baseline_psi, 1, function(x) {
    x <- na.omit(x)
    length(x) > 1 && var(x) != 0
  }))

  whitelist1 <- unique(c(tmp1, tmp2))
  case_psi <- case_psi[whitelist1, , drop = FALSE]
  baseline_psi <- baseline_psi[whitelist1, , drop = FALSE]
  diff_table <- diff_table[whitelist1, , drop = FALSE]

  # Filter by minimum non-NA sample count per group
  whitelist2 <- which(
    rowSums(!is.na(case_psi)) > cutoff_case &
      rowSums(!is.na(baseline_psi)) > cutoff_baseline
  )

  case_psi <- case_psi[whitelist2, , drop = FALSE]
  baseline_psi <- baseline_psi[whitelist2, , drop = FALSE]
  diff_table <- diff_table[whitelist2, , drop = FALSE]

  rnames <- paste(diff_table$GENE, diff_table$EVENT, sep = "_")
  rownames(baseline_psi) <- rnames
  rownames(case_psi) <- rnames
  rownames(diff_table) <- rnames

  diff_table$mean <- 0
  diff_table$p.val <- 0
  diff_table$padj.val <- 0

  for (i in seq_len(nrow(diff_table))) {
    case_i <- na.omit(unlist(case_psi[i, ]))
    baseline_i <- na.omit(unlist(baseline_psi[i, ]))

    test <- suppressWarnings(wilcox.test(baseline_i, case_i))
    delta <- mean(case_i, na.rm = TRUE) - mean(baseline_i, na.rm = TRUE)
    diff_table[i, c("mean", "p.val")] <- c(delta, test$p.value)
  }

  diff_table$padj.val <- p.adjust(diff_table$p.val, method = "BH")

  rds_path <- file.path(outdir, paste0("diff_table_", name, ".RDS"))
  saveRDS(diff_table, rds_path)
  write.table(diff_table, file.path(outdir, paste0("diff_table_", name, "_all.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  return(diff_table)
}

# Load inputs
message(">>> Reading VAST-TOOLS inclusion table: ", opts$input)
total_table <- read.delim(opts$input, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

message(">>> Reading sample metadata: ", opts$metadata)
meta <- read.delim(opts$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta$condition <- factor(meta$condition)

# Run Wilcoxon differential splicing
out_diff <- microfilt_wilcox(
  total_table = total_table,
  config_rep = meta,
  baseline = opts$baseline,
  case = opts$case,
  name = opts$name,
  cutoff_baseline = opts$cutoff_baseline,
  cutoff_case = opts$cutoff_case,
  outdir = opts$outdir
)

# Stratify significant differential splicing events (|dPSI| > 10, p < 0.05) by exon length
valid_complex <- c("S", "C1", "C2", "C3", "MIC", "ANN")

sig_1_15 <- subset(out_diff, COMPLEX %in% valid_complex & LENGTH >= 1 & LENGTH <= 15 & p.val < 0.05 & abs(mean) > 10)
sig_16_27 <- subset(out_diff, COMPLEX %in% valid_complex & LENGTH >= 16 & LENGTH <= 27 & p.val < 0.05 & abs(mean) > 10)
sig_28_51 <- subset(out_diff, COMPLEX %in% valid_complex & LENGTH >= 28 & LENGTH <= 51 & p.val < 0.05 & abs(mean) > 10)
sig_le_51 <- subset(out_diff, COMPLEX %in% valid_complex & LENGTH <= 51 & p.val < 0.05 & abs(mean) > 10)

write.table(sig_1_15, file.path(opts$outdir, paste0(opts$name, "_sig_1_15nt.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_16_27, file.path(opts$outdir, paste0(opts$name, "_sig_16_27nt.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_28_51, file.path(opts$outdir, paste0(opts$name, "_sig_28_51nt.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_le_51, file.path(opts$outdir, paste0(opts$name, "_sig_le_51nt.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("Summary (%s vs %s):", opts$case, opts$baseline))
message(sprintf("  1-15 nt microexons:  %d skipped (dPSI < -10), %d included (dPSI > 10)", sum(sig_1_15$mean < 0), sum(sig_1_15$mean > 0)))
message(sprintf("  16-27 nt microexons: %d skipped (dPSI < -10), %d included (dPSI > 10)", sum(sig_16_27$mean < 0), sum(sig_16_27$mean > 0)))
message(sprintf("  28-51 nt short exons: %d skipped (dPSI < -10), %d included (dPSI > 10)", sum(sig_28_51$mean < 0), sum(sig_28_51$mean > 0)))
message(sprintf("  <= 51 nt combined:   %d skipped (dPSI < -10), %d included (dPSI > 10)", sum(sig_le_51$mean < 0), sum(sig_le_51$mean > 0)))
