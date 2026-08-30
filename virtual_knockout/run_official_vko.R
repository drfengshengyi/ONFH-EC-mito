#!/usr/bin/env Rscript

# Donor-separated SQSTM1 virtual knockout using the authors' official
# scTenifoldKnk R package.
#
# This is an independent implementation audit of the existing Python-port
# analysis. It intentionally freezes the same 300-gene feature set and the
# same endothelial-cell barcode lists so differences primarily reflect the
# official R workflow. Important official-package behavior is retained:
# CPM normalization, cell subsampling with replacement, and internal seed 1.

options(stringsAsFactors = FALSE, warn = 1)

parse_cli <- function(args) {
  config <- list(
    profile = "manuscript",
    cores = 1L,
    donors = c("hoa2", "hoa3"),
    output_dir = NULL,
    target = "SQSTM1",
    exclude_mt_encoded = FALSE
  )
  for (arg in args) {
    if (grepl("^--profile=", arg)) {
      config$profile <- sub("^--profile=", "", arg)
    } else if (grepl("^--cores=", arg)) {
      config$cores <- as.integer(sub("^--cores=", "", arg))
    } else if (grepl("^--donors=", arg)) {
      value <- sub("^--donors=", "", arg)
      config$donors <- strsplit(value, ",", fixed = TRUE)[[1L]]
    } else if (grepl("^--output-dir=", arg)) {
      config$output_dir <- sub("^--output-dir=", "", arg)
    } else if (grepl("^--target=", arg)) {
      config$target <- sub("^--target=", "", arg)
    } else if (identical(arg, "--exclude-mt-encoded")) {
      config$exclude_mt_encoded <- TRUE
    } else if (arg %in% c("--help", "-h")) {
      cat(
        "Usage:\n",
        "  Rscript virtual_knockout/run_official_vko.R --profile=manuscript --cores=1\n\n",
        "Profiles:\n",
        "  manuscript       20 networks, q=0.95, tensor rank 3, MA d=30\n",
        "  official-default 10 networks, q=0.90, tensor rank 3, MA d=2\n\n",
        "Optional:\n",
        "  --donors=hoa2,hoa3\n",
        "  --output-dir=results/custom-vko-output\n",
        "  --target=SQSTM1\n",
        "  --exclude-mt-encoded  Refit after removing genes beginning MT-\n",
        sep = ""
      )
      quit(save = "no", status = 0L)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  if (is.na(config$cores) || config$cores < 1L) {
    stop("--cores must be a positive integer", call. = FALSE)
  }
  config
}

profile_parameters <- function(profile) {
  if (identical(profile, "manuscript")) {
    return(list(
      n_net = 20L,
      n_cells = 500L,
      n_comp = 3L,
      q = 0.95,
      td_k = 3L,
      td_max_iter = 500L,
      td_max_error = 1e-5,
      td_n_decimal = 3L,
      ma_n_dim = 30L
    ))
  }
  if (identical(profile, "official-default")) {
    return(list(
      n_net = 10L,
      n_cells = 500L,
      n_comp = 3L,
      q = 0.90,
      td_k = 3L,
      td_max_iter = 1000L,
      td_max_error = 1e-5,
      td_n_decimal = 3L,
      ma_n_dim = 2L
    ))
  }
  stop(
    "Unknown profile '", profile,
    "'. Use manuscript or official-default.",
    call. = FALSE
  )
}

script_location <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(
      dirname(sub("^--file=", "", file_arg[[1L]])),
      winslash = "/",
      mustWork = TRUE
    ))
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    return(normalizePath(
      dirname(frame_files[[length(frame_files)]]),
      winslash = "/",
      mustWork = TRUE
    ))
  }

  project_from_env <- Sys.getenv("ONFH_PROJECT_ROOT", unset = "")
  candidates <- c(
    if (nzchar(project_from_env)) file.path(project_from_env, "virtual_knockout") else NULL,
    file.path(getwd(), "virtual_knockout"),
    getwd()
  )
  expected_file <- "vko_selected_features_v5.csv"
  valid <- candidates[file.exists(file.path(candidates, expected_file))]
  if (length(valid) > 0L) {
    return(normalizePath(valid[[1L]], winslash = "/", mustWork = TRUE))
  }
  stop(
    "Cannot locate the virtual_knockout directory. Set ONFH_PROJECT_ROOT or run the ",
    "script from the repository root.",
    call. = FALSE
  )
}

invocation_args <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    commandArgs(trailingOnly = TRUE)
  } else {
    character(0)
  }
}

assert_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Run environment/install_official_sctenifoldknk.R first.",
      call. = FALSE
    )
  }
}

standardize_dr <- function(frame, donor) {
  required <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
  absent <- setdiff(required, colnames(frame))
  if (length(absent) > 0L) {
    stop(
      "Unexpected official-package output; missing columns: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  frame <- frame[order(frame$p.value, frame$gene), , drop = FALSE]
  data.frame(
    donor = donor,
    rank = seq_len(nrow(frame)),
    Gene = as.character(frame$gene),
    Distance = as.numeric(frame$distance),
    Z = as.numeric(frame$Z),
    FC = as.numeric(frame$FC),
    p_value = as.numeric(frame$p.value),
    adjusted_p_value = as.numeric(frame$p.adj),
    check.names = FALSE
  )
}

make_consensus <- function(results, donors) {
  rank_tables <- lapply(donors, function(donor) {
    frame <- results[[donor]][, c("Gene", "rank", "p_value", "adjusted_p_value")]
    names(frame)[-1L] <- paste0(names(frame)[-1L], "_", donor)
    frame
  })
  merged <- Reduce(
    function(x, y) merge(x, y, by = "Gene", all = FALSE, sort = TRUE),
    rank_tables
  )
  rank_columns <- paste0("rank_", donors)
  merged$mean_rank <- rowMeans(merged[, rank_columns, drop = FALSE])
  merged$max_rank <- apply(merged[, rank_columns, drop = FALSE], 1L, max)
  merged[order(merged$mean_rank, merged$max_rank, merged$Gene), , drop = FALSE]
}

compare_python_ranks <- function(r_consensus, python_path, donors) {
  if (!file.exists(python_path)) {
    return(NULL)
  }
  python <- utils::read.csv(python_path, check.names = FALSE)
  rows <- lapply(donors, function(donor) {
    column <- paste0("rank_", donor)
    if (!column %in% names(python) || !column %in% names(r_consensus)) {
      return(NULL)
    }
    joined <- merge(
      r_consensus[, c("Gene", column)],
      python[, c("Gene", column)],
      by = "Gene",
      suffixes = c("_R", "_Python")
    )
    test <- suppressWarnings(stats::cor.test(
      joined[[paste0(column, "_R")]],
      joined[[paste0(column, "_Python")]],
      method = "spearman",
      exact = FALSE
    ))
    data.frame(
      donor = donor,
      n_genes = nrow(joined),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

cli <- parse_cli(invocation_args())
params <- profile_parameters(cli$profile)
analysis_dir <- script_location()
project_root <- normalizePath(
  file.path(analysis_dir, ".."), winslash = "/", mustWork = TRUE
)
data_dir <- file.path(project_root, "data", "liao2022")
features_path <- file.path(analysis_dir, "vko_selected_features_v5.csv")
barcodes_path <- file.path(analysis_dir, "vko_ec_barcodes_v5.csv")
python_consensus_path <- file.path(
  analysis_dir, "vko_sqstm1_consensus_v5.csv"
)

output_dir <- cli$output_dir
if (is.null(output_dir)) {
  suffix <- if (isTRUE(cli$exclude_mt_encoded)) "_no_mt_encoded" else ""
  output_dir <- file.path(
    project_root,
    "results",
    paste0("official_r_vko_", gsub("-", "_", cli$profile), suffix)
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

display_path <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(project_root, "/")
  if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else normalized
}

assert_packages(c("Matrix", "scTenifoldKnk", "jsonlite"))

for (path in c(features_path, barcodes_path)) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }
}

target <- cli$target
features <- utils::read.csv(features_path, check.names = FALSE)
if (!all(c("gene", "selected_order") %in% names(features))) {
  stop("Malformed selected-feature file: ", features_path, call. = FALSE)
}
selected_genes <- as.character(features$gene)
if (length(selected_genes) != 300L || anyDuplicated(selected_genes)) {
  stop("Expected exactly 300 unique genes before sensitivity filtering", call. = FALSE)
}
mt_encoded_genes <- selected_genes[grepl("^MT-", selected_genes)]
if (isTRUE(cli$exclude_mt_encoded)) {
  selected_genes <- setdiff(selected_genes, mt_encoded_genes)
  if (length(mt_encoded_genes) == 0L) {
    stop("No MT-encoded genes were available for exclusion", call. = FALSE)
  }
}
if (!target %in% selected_genes) {
  stop(target, " is absent from the frozen feature set", call. = FALSE)
}

barcode_manifest <- utils::read.csv(barcodes_path, check.names = FALSE)
required_barcode_columns <- c("donor", "ec_obs_name", "source_barcode")
if (!all(required_barcode_columns %in% names(barcode_manifest))) {
  stop("Malformed barcode manifest: ", barcodes_path, call. = FALSE)
}

results <- list()
run_audit <- list()

for (donor in cli$donors) {
  message("\n===== ", toupper(donor), " =====")
  source_path <- file.path(data_dir, paste0(donor, ".rds"))
  if (!file.exists(source_path)) {
    stop("Missing source count matrix: ", source_path, call. = FALSE)
  }
  count_matrix <- readRDS(source_path)
  if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
    stop("Source matrix lacks gene/cell names: ", source_path, call. = FALSE)
  }
  missing_genes <- setdiff(selected_genes, rownames(count_matrix))
  if (length(missing_genes) > 0L) {
    stop(
      donor, " is missing selected genes: ",
      paste(utils::head(missing_genes, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  donor_barcodes <- barcode_manifest$source_barcode[
    barcode_manifest$donor == donor
  ]
  donor_barcodes <- as.character(donor_barcodes)
  if (length(donor_barcodes) == 0L || anyDuplicated(donor_barcodes)) {
    stop("Missing or duplicated barcodes for ", donor, call. = FALSE)
  }
  missing_barcodes <- setdiff(donor_barcodes, colnames(count_matrix))
  if (length(missing_barcodes) > 0L) {
    stop(
      donor, " is missing source barcodes: ",
      paste(utils::head(missing_barcodes, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  expression <- count_matrix[selected_genes, donor_barcodes, drop = FALSE]
  if (ncol(expression) < params$n_cells) {
    stop(
      donor, " has only ", ncol(expression),
      " ECs, fewer than requested nc_nCells=", params$n_cells,
      call. = FALSE
    )
  }
  target_prevalence <- as.numeric(
    Matrix::rowSums(expression[target, , drop = FALSE] > 0) / ncol(expression)
  )

  started <- proc.time()[["elapsed"]]
  model <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = expression,
    gKO = target,
    qc = FALSE,
    nc_lambda = 0,
    nc_nNet = params$n_net,
    nc_nCells = params$n_cells,
    nc_nComp = params$n_comp,
    nc_scaleScores = TRUE,
    nc_symmetric = FALSE,
    nc_q = params$q,
    nc_priorNetwork = NULL,
    td_K = params$td_k,
    td_maxIter = params$td_max_iter,
    td_maxError = params$td_max_error,
    td_nDecimal = params$td_n_decimal,
    ma_nDim = params$ma_n_dim,
    nCores = cli$cores
  )
  elapsed <- proc.time()[["elapsed"]] - started

  donor_result <- standardize_dr(model$diffRegulation, donor)
  results[[donor]] <- donor_result
  utils::write.csv(
    donor_result,
    file.path(
      output_dir,
      paste0("vko_", tolower(target), "_", donor, "_official_r.csv")
    ),
    row.names = FALSE
  )
  saveRDS(
    model,
    file.path(
      output_dir,
      paste0("vko_", tolower(target), "_", donor, "_official_r.rds")
    ),
    compress = "xz"
  )

  target_rank <- donor_result$rank[donor_result$Gene == target]
  downstream_fdr_count <- sum(
    donor_result$adjusted_p_value < 0.05 & donor_result$Gene != target,
    na.rm = TRUE
  )
  run_audit[[donor]] <- list(
    donor = donor,
    n_ec = ncol(expression),
    n_genes = nrow(expression),
    target_prevalence = target_prevalence,
    target_rank = unname(target_rank),
    downstream_genes_fdr_lt_0_05 = downstream_fdr_count,
    elapsed_seconds = unname(elapsed)
  )
}

consensus <- make_consensus(results, cli$donors)
utils::write.csv(
  consensus,
  file.path(
    output_dir,
    paste0("vko_", tolower(target), "_consensus_official_r.csv")
  ),
  row.names = FALSE
)

donor_rank_test <- NULL
if (length(cli$donors) == 2L) {
  x_column <- paste0("rank_", cli$donors[[1L]])
  y_column <- paste0("rank_", cli$donors[[2L]])
  donor_rank_test <- suppressWarnings(stats::cor.test(
    consensus[[x_column]], consensus[[y_column]],
    method = "spearman", exact = FALSE
  ))
}

python_comparison <- if (
  identical(target, "SQSTM1") && !isTRUE(cli$exclude_mt_encoded)
) {
  compare_python_ranks(consensus, python_consensus_path, cli$donors)
} else {
  NULL
}
if (!is.null(python_comparison) && nrow(python_comparison) > 0L) {
  utils::write.csv(
    python_comparison,
    file.path(output_dir, "official_r_vs_python_rank_agreement.csv"),
    row.names = FALSE
  )
}

package_versions <- list(
  R = paste(R.version$major, R.version$minor, sep = "."),
  scTenifoldKnk = as.character(utils::packageVersion("scTenifoldKnk")),
  scTenifoldNet = as.character(utils::packageVersion("scTenifoldNet")),
  Matrix = as.character(utils::packageVersion("Matrix")),
  jsonlite = as.character(utils::packageVersion("jsonlite"))
)

provenance <- list(
  analysis = paste(
    "donor-separated", target,
    "virtual knockout with official R package"
  ),
  target = target,
  profile = cli$profile,
  donors = cli$donors,
  frozen_feature_file = display_path(features_path),
  frozen_barcode_file = display_path(barcodes_path),
  feature_sensitivity = list(
    exclude_mt_encoded = isTRUE(cli$exclude_mt_encoded),
    excluded_genes = if (isTRUE(cli$exclude_mt_encoded)) {
      mt_encoded_genes
    } else {
      character(0)
    },
    retained_gene_count = length(selected_genes),
    design = if (isTRUE(cli$exclude_mt_encoded)) {
      paste(
        "Complete network refit after removing all frozen-feature",
        "symbols beginning MT-."
      )
    } else {
      "Original frozen 300-gene feature set."
    }
  ),
  parameters = params,
  n_cores = cli$cores,
  official_package_behavior = list(
    normalization = "CPM normalization is always applied by scTenifoldKnk",
    sampling = "scTenifoldNet::makeNetworks samples cells with replacement",
    seed = "the official scTenifoldKnk wrapper internally calls set.seed(1)"
  ),
  runs = unname(run_audit),
  donor_rank_spearman = if (is.null(donor_rank_test)) NULL else list(
    rho = unname(donor_rank_test$estimate),
    p_value = donor_rank_test$p.value
  ),
  software = package_versions,
  interpretation_guardrail = paste(
    "Network perturbation ranks are exploratory predictions, not",
    "experimental knockout validation or participant-level causal tests."
  )
)

jsonlite::write_json(
  provenance,
  file.path(output_dir, "vko_provenance_official_r.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

summary_lines <- c(
  paste0("Profile: ", cli$profile),
  paste0("Official scTenifoldKnk version: ", package_versions$scTenifoldKnk),
  vapply(run_audit, function(x) {
    paste0(
      toupper(x$donor), ": n=", x$n_ec,
      ", ", target, " prevalence=", sprintf("%.1f%%", 100 * x$target_prevalence),
      ", ", target, " rank=", x$target_rank,
      ", downstream FDR<0.05 genes=", x$downstream_genes_fdr_lt_0_05
    )
  }, character(1)),
  if (is.null(donor_rank_test)) NULL else paste0(
    "Cross-donor Spearman rho=", sprintf("%.4f", donor_rank_test$estimate),
    ", p=", format(donor_rank_test$p.value, scientific = TRUE, digits = 4)
  ),
  paste0("Output directory: ", display_path(output_dir))
)
writeLines(summary_lines, file.path(output_dir, "RUN_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
