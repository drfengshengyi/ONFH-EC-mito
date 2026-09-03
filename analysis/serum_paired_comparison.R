#!/usr/bin/env Rscript

# Paired comparison of two archived, aggregated cross-fitted serum predictions.
# This test compares predictions for the same 40 participants. It does not
# include uncertainty from refitting either nested-CV pipeline and therefore
# remains an internal, exploratory comparison.

suppressPackageStartupMessages({
  library(pROC)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
root <- dirname(dirname(script_file))
input_dir <- Sys.getenv(
  "ONFH_OUTPUT_DIR",
  unset = file.path(root, "results", "figure_inputs")
)
input_file <- file.path(input_dir, "diag_oof_predictions_aggregated_v8.csv")
output_file <- file.path(input_dir, "diag_delong_model_comparison_v8.csv")

if (!file.exists(input_file)) {
  stop("Missing paired prediction file: ", input_file)
}

x <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c(
  "sample", "label", "candidate_mean_oof_probability", "ma_mean_oof_probability"
)
missing <- setdiff(required, names(x))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(x$sample)) stop("Paired comparison requires one row per participant")
if (!identical(sort(unique(x$label)), c(0L, 1L))) stop("label must be binary 0/1")

roc_candidate <- roc(
  response = x$label,
  predictor = x$candidate_mean_oof_probability,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)
roc_ma <- roc(
  response = x$label,
  predictor = x$ma_mean_oof_probability,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)
test <- roc.test(roc_candidate, roc_ma, method = "delong", paired = TRUE)

out <- data.frame(
  comparison = "mitochondrial_candidate_space_minus_Ma_four_gene_model",
  n_participants = nrow(x),
  n_sonfh = sum(x$label == 1),
  n_steroid_exposed_non_sonfh = sum(x$label == 0),
  candidate_auc = as.numeric(auc(roc_candidate)),
  ma_auc = as.numeric(auc(roc_ma)),
  delta_auc = as.numeric(auc(roc_candidate) - auc(roc_ma)),
  delong_z = as.numeric(test$statistic),
  delong_p_two_sided = as.numeric(test$p.value),
  paired = TRUE,
  scope = paste(
    "paired DeLong test on fixed probabilities averaged across five repeated",
    "outer-CV runs; model-refitting variability excluded"
  ),
  stringsAsFactors = FALSE
)
write.csv(out, output_file, row.names = FALSE)
cat(normalizePath(output_file, winslash = "/"), "\n")
