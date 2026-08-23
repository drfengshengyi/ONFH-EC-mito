#!/usr/bin/env Rscript

# Export only the official scTenifoldKnk network and manifold objects needed
# for the publication figure. No new model is fitted here.

suppressPackageStartupMessages(library(Matrix))

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "virtual_knockout/export_figure_data.R"
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
input_dir <- file.path(results_dir, "official_r_vko_manuscript")
output_dir <- file.path(results_dir, "official_r_vko_figure_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

target <- "SQSTM1"
donors <- c("hoa2", "hoa3")

for (donor in donors) {
  model_path <- file.path(input_dir, paste0("vko_sqstm1_", donor, "_official_r.rds"))
  model <- readRDS(model_path)

  wt <- as.matrix(model$tensorNetworks$WT)
  ko <- as.matrix(model$tensorNetworks$KO)
  stopifnot(target %in% rownames(wt), target %in% rownames(ko))
  stopifnot(all(abs(ko[target, ]) < .Machine$double.eps^0.5))

  edge_table <- data.frame(
    donor = donor,
    gene = colnames(wt),
    wt_outgoing_weight = as.numeric(wt[target, ]),
    ko_outgoing_weight = as.numeric(ko[target, ]),
    absolute_wt_weight = abs(as.numeric(wt[target, ])),
    stringsAsFactors = FALSE
  )
  edge_table <- edge_table[order(-edge_table$absolute_wt_weight, edge_table$gene), ]
  utils::write.csv(
    edge_table,
    file.path(output_dir, paste0("sqstm1_outgoing_edges_", donor, ".csv")),
    row.names = FALSE
  )

  manifold <- as.matrix(model$manifoldAlignment)
  stopifnot(nrow(manifold) == 2L * nrow(wt))
  embedding <- stats::prcomp(manifold, center = TRUE, scale. = FALSE)$x[, 1:2, drop = FALSE]
  state <- ifelse(startsWith(rownames(manifold), "X_"), "WT", "KO")
  gene <- sub("^[XY]_", "", rownames(manifold))
  manifold_table <- data.frame(
    donor = donor,
    state = state,
    gene = gene,
    PC1 = embedding[, 1],
    PC2 = embedding[, 2],
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifold_table,
    file.path(output_dir, paste0("wt_ko_manifold_", donor, ".csv")),
    row.names = FALSE
  )
}

message("Figure data exported to: ", output_dir)
