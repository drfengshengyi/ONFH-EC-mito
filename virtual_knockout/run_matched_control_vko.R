#!/usr/bin/env Rscript

# Matched-comparator virtual knockouts for the official-R SQSTM1 analysis.
#
# The donor-specific WT networks are frozen before comparator selection. The
# official wrapper constructs WT independently of gKO and then zeroes the
# selected gene's outgoing row. Reusing each frozen WT therefore isolates the
# perturbation target and avoids introducing a new network-resampling draw for
# every comparator. These genes are matched comparators, not experimentally
# validated negative controls.

options(stringsAsFactors = FALSE, warn = 1)

parse_cli <- function(args) {
  config <- list(
    cores = 1L,
    output_dir = NULL
  )
  for (arg in args) {
    if (grepl("^--cores=", arg)) {
      config$cores <- as.integer(sub("^--cores=", "", arg))
    } else if (grepl("^--output-dir=", arg)) {
      config$output_dir <- sub("^--output-dir=", "", arg)
    } else if (arg %in% c("--help", "-h")) {
      cat(
        "Usage:\n",
        "  Rscript virtual_knockout/run_matched_control_vko.R --cores=1\n\n",
        "Optional:\n",
        "  --output-dir=results/official_r_vko_matched_controls\n",
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

script_location <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(
      dirname(sub("^--file=", "", file_arg[[1L]])),
      winslash = "/",
      mustWork = TRUE
    ))
  }
  normalizePath(file.path(getwd(), "virtual_knockout"), winslash = "/", mustWork = TRUE)
}

invocation_args <- function() {
  if (any(grepl("^--file=", commandArgs(trailingOnly = FALSE)))) {
    commandArgs(trailingOnly = TRUE)
  } else {
    character(0)
  }
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

standardize_dr <- function(frame, donor, profile, perturbation_gene, role) {
  required <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
  if (!all(required %in% names(frame))) {
    stop("Unexpected dRegulation output", call. = FALSE)
  }
  frame <- frame[order(frame$p.value, frame$gene), , drop = FALSE]
  data.frame(
    profile = profile,
    donor = donor,
    perturbation_gene = perturbation_gene,
    perturbation_role = role,
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

run_perturbation <- function(wt, gene, d, cores, donor, profile, role) {
  if (!gene %in% rownames(wt)) {
    stop(gene, " is absent from ", profile, " ", donor, " WT", call. = FALSE)
  }
  ko <- wt
  ko[gene, ] <- 0
  stopifnot(all(ko[gene, ] == 0))
  set.seed(1)
  started <- proc.time()[["elapsed"]]
  ma <- scTenifoldNet::manifoldAlignment(wt, ko, d = d, nCores = cores)
  dr <- scTenifoldKnk::dRegulation(ma)
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    result = standardize_dr(dr, donor, profile, gene, role),
    elapsed_seconds = unname(elapsed)
  )
}

rank_summary <- function(frame, perturbation_gene) {
  donors <- c("hoa2", "hoa3")
  split_result <- split(frame, frame$donor)
  if (!all(donors %in% names(split_result))) {
    stop("Both donors are required for matched-control calibration", call. = FALSE)
  }
  ranks <- lapply(split_result[donors], function(x) {
    x <- x[x$Gene != perturbation_gene, c("Gene", "rank", "adjusted_p_value")]
    rownames(x) <- x$Gene
    x
  })
  common <- intersect(ranks[[1L]]$Gene, ranks[[2L]]$Gene)
  r2 <- ranks[[1L]][common, "rank"]
  r3 <- ranks[[2L]][common, "rank"]
  rho_test <- suppressWarnings(stats::cor.test(
    r2, r3, method = "spearman", exact = FALSE
  ))
  top_genes <- lapply(ranks, function(x) {
    x$Gene[order(x$rank, x$Gene)][seq_len(min(20L, nrow(x)))]
  })
  significant <- lapply(ranks, function(x) {
    x$Gene[x$adjusted_p_value < 0.05]
  })
  data.frame(
    cross_donor_spearman = unname(rho_test$estimate),
    cross_donor_spearman_p = rho_test$p.value,
    top20_overlap_count = length(intersect(top_genes[[1L]], top_genes[[2L]])),
    top20_overlap_fraction = length(intersect(top_genes[[1L]], top_genes[[2L]])) / 20,
    downstream_fdr_hoa2 = length(significant[[1L]]),
    downstream_fdr_hoa3 = length(significant[[2L]]),
    replicated_downstream_fdr_count = length(intersect(significant[[1L]], significant[[2L]])),
    stringsAsFactors = FALSE
  )
}

cli <- parse_cli(invocation_args())
assert_packages(c("Matrix", "jsonlite", "scTenifoldKnk", "scTenifoldNet"))

analysis_dir <- script_location()
project_root <- normalizePath(file.path(analysis_dir, ".."), winslash = "/", mustWork = TRUE)
features_path <- file.path(analysis_dir, "vko_selected_features_v5.csv")
output_dir <- cli$output_dir
if (is.null(output_dir)) {
  output_dir <- file.path(project_root, "results", "official_r_vko_matched_controls")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

profiles <- list(
  manuscript = list(directory = "official_r_vko_manuscript", d = 30L),
  official_default = list(directory = "official_r_vko_official_default", d = 2L)
)
donors <- c("hoa2", "hoa3")
target <- "SQSTM1"

models <- lapply(profiles, function(spec) {
  setNames(lapply(donors, function(donor) {
    path <- file.path(
      project_root,
      "results",
      spec$directory,
      paste0("vko_sqstm1_", donor, "_official_r.rds")
    )
    if (!file.exists(path)) stop("Missing frozen model: ", path, call. = FALSE)
    readRDS(path)
  }), donors)
})

features <- utils::read.csv(features_path, check.names = FALSE)
if (nrow(features) != 300L || anyDuplicated(features$gene)) {
  stop("Expected 300 unique frozen features", call. = FALSE)
}

# Match only on pre-perturbation control data and the manuscript-profile WT
# topology. No disease effect or virtual-knockout result enters selection.
for (donor in donors) {
  wt <- as.matrix(models$manuscript[[donor]]$tensorNetworks$WT)
  wt <- wt[features$gene, features$gene, drop = FALSE]
  features[[paste0("wt_out_degree_", donor)]] <- rowSums(abs(wt) > 1e-12)
  features[[paste0("wt_out_strength_", donor)]] <- rowSums(abs(wt))
}

target_row <- features[features$gene == target, , drop = FALSE]
candidate <- features[
  features$gene != target &
    !grepl("^MT-", features$gene) &
    features$prevalence_hoa2 > 0 &
    features$prevalence_hoa3 > 0,
  ,
  drop = FALSE
]

caliper_pass <- rep(TRUE, nrow(candidate))
for (donor in donors) {
  caliper_pass <- caliper_pass &
    abs(candidate[[paste0("prevalence_", donor)]] -
          target_row[[paste0("prevalence_", donor)]]) <= 0.20 &
    abs(candidate[[paste0("mean_log1p_cpm_", donor)]] -
          target_row[[paste0("mean_log1p_cpm_", donor)]]) <= 0.60 &
    abs(candidate[[paste0("wt_out_degree_", donor)]] -
          target_row[[paste0("wt_out_degree_", donor)]]) <= 40 &
    candidate[[paste0("wt_out_strength_", donor)]] /
      target_row[[paste0("wt_out_strength_", donor)]] >= 1 / 3 &
    candidate[[paste0("wt_out_strength_", donor)]] /
      target_row[[paste0("wt_out_strength_", donor)]] <= 3
}
matched <- candidate[caliper_pass, , drop = FALSE]

expected_controls <- c(
  "AKAP12", "ATP5F1A", "C11orf96", "CALCOCO2", "EGR1", "FIS1",
  "GPX4", "HSPE1", "LGALS1", "MAP1LC3B", "MT2A", "NDUFA1",
  "NFE2L2", "PRDX1", "PRDX2", "PRDX5", "RHEB", "SOD1", "UQCRB",
  "VDAC2"
)
if (!setequal(matched$gene, expected_controls)) {
  stop(
    "Frozen matched-comparator set changed. Observed: ",
    paste(sort(matched$gene), collapse = ", "),
    call. = FALSE
  )
}

distance_columns <- c(
  "prevalence_hoa2", "prevalence_hoa3",
  "mean_log1p_cpm_hoa2", "mean_log1p_cpm_hoa3",
  "dispersion_pct_hoa2", "dispersion_pct_hoa3",
  "wt_out_degree_hoa2", "wt_out_degree_hoa3",
  "wt_out_strength_hoa2", "wt_out_strength_hoa3"
)
distance_pool <- rbind(target_row, candidate)
distance_values <- distance_pool[, distance_columns, drop = FALSE]
distance_values$wt_out_strength_hoa2 <- log1p(distance_values$wt_out_strength_hoa2)
distance_values$wt_out_strength_hoa3 <- log1p(distance_values$wt_out_strength_hoa3)
scaled <- scale(distance_values)
target_scaled <- scaled[1L, ]
distance_pool$standardized_distance <- sqrt(rowSums(
  (scaled - matrix(
    target_scaled,
    nrow = nrow(scaled),
    ncol = ncol(scaled),
    byrow = TRUE
  ))^2
))
matched <- merge(
  matched,
  distance_pool[, c("gene", "standardized_distance")],
  by = "gene",
  sort = FALSE
)
matched <- matched[order(matched$standardized_distance, matched$gene), , drop = FALSE]
matched$closest_five <- seq_len(nrow(matched)) <= 5L

matching_output <- rbind(
  transform(target_row, standardized_distance = 0, closest_five = FALSE),
  matched[, names(transform(target_row, standardized_distance = 0, closest_five = FALSE)), drop = FALSE]
)
matching_output$role <- ifelse(matching_output$gene == target, "target", "matched comparator")
matching_output <- matching_output[, c(
  "role", "gene", "standardized_distance", "closest_five",
  "prevalence_hoa2", "prevalence_hoa3",
  "mean_log1p_cpm_hoa2", "mean_log1p_cpm_hoa3",
  "dispersion_pct_hoa2", "dispersion_pct_hoa3",
  "wt_out_degree_hoa2", "wt_out_degree_hoa3",
  "wt_out_strength_hoa2", "wt_out_strength_hoa3",
  "in_prespecified_panel"
)]
utils::write.csv(
  matching_output,
  file.path(output_dir, "matched_control_selection.csv"),
  row.names = FALSE
)

design <- data.frame(
  item = c(
    "candidate universe", "excluded symbols", "prevalence caliper",
    "mean log1p CPM caliper", "WT out-degree caliper",
    "WT absolute out-strength ratio", "number of comparators",
    "selection blinding", "interpretation"
  ),
  specification = c(
    "Frozen 300-gene feature set; detectable in both donors",
    "SQSTM1 and symbols beginning MT-",
    "Absolute donor-specific difference <= 0.20",
    "Absolute donor-specific difference <= 0.60",
    "Absolute donor-specific difference <= 40 edges",
    "Between 1/3 and 3 in each donor",
    "20 (minimum attainable empirical p = 1/21)",
    "Control covariates only; no disease effect or KO outcome used",
    "Matched comparators, not validated biological negative controls"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(design, file.path(output_dir, "matched_control_design.csv"), row.names = FALSE)

all_genes <- c(target, matched$gene)
all_results <- list()
runtime_rows <- list()
validation_rows <- list()
result_index <- 0L
runtime_index <- 0L
validation_index <- 0L

for (profile in names(profiles)) {
  spec <- profiles[[profile]]
  for (donor in donors) {
    wt <- as.matrix(models[[profile]][[donor]]$tensorNetworks$WT)
    for (gene in all_genes) {
      role <- if (gene == target) "target" else "matched comparator"
      rerun <- run_perturbation(
        wt = wt,
        gene = gene,
        d = spec$d,
        cores = cli$cores,
        donor = donor,
        profile = profile,
        role = role
      )
      result_index <- result_index + 1L
      all_results[[result_index]] <- rerun$result
      runtime_index <- runtime_index + 1L
      runtime_rows[[runtime_index]] <- data.frame(
        profile = profile,
        donor = donor,
        perturbation_gene = gene,
        elapsed_seconds = rerun$elapsed_seconds
      )

      if (gene == target) {
        stored <- standardize_dr(
          models[[profile]][[donor]]$diffRegulation,
          donor,
          profile,
          target,
          "target"
        )
        joined <- merge(
          rerun$result[, c("Gene", "rank", "p_value")],
          stored[, c("Gene", "rank", "p_value")],
          by = "Gene",
          suffixes = c("_rerun", "_stored")
        )
        validation_index <- validation_index + 1L
        validation_rows[[validation_index]] <- data.frame(
          profile = profile,
          donor = donor,
          n_genes = nrow(joined),
          ranks_identical = identical(joined$rank_rerun, joined$rank_stored),
          max_abs_p_value_difference = max(abs(
            joined$p_value_rerun - joined$p_value_stored
          ))
        )
      }
    }
  }
}

gene_ranks <- do.call(rbind, all_results)
utils::write.csv(
  gene_ranks,
  file.path(output_dir, "matched_control_vko_gene_ranks.csv"),
  row.names = FALSE
)
utils::write.csv(
  do.call(rbind, runtime_rows),
  file.path(output_dir, "matched_control_vko_runtime.csv"),
  row.names = FALSE
)
validation <- do.call(rbind, validation_rows)
utils::write.csv(
  validation,
  file.path(output_dir, "frozen_wt_reuse_validation.csv"),
  row.names = FALSE
)
if (!all(validation$ranks_identical) || any(validation$max_abs_p_value_difference > 1e-12)) {
  stop("Frozen-WT rerun did not reproduce the stored SQSTM1 result", call. = FALSE)
}

summary_rows <- list()
summary_index <- 0L
for (profile in names(profiles)) {
  profile_frame <- gene_ranks[gene_ranks$profile == profile, , drop = FALSE]
  for (gene in all_genes) {
    summary_index <- summary_index + 1L
    one <- profile_frame[profile_frame$perturbation_gene == gene, , drop = FALSE]
    metrics <- rank_summary(one, gene)
    summary_rows[[summary_index]] <- cbind(
      data.frame(
        profile = profile,
        perturbation_gene = gene,
        role = if (gene == target) "target" else "matched comparator",
        stringsAsFactors = FALSE
      ),
      metrics
    )
  }
}
summary_table <- do.call(rbind, summary_rows)
utils::write.csv(
  summary_table,
  file.path(output_dir, "matched_control_vko_summary.csv"),
  row.names = FALSE
)

metric_directions <- c(
  cross_donor_spearman = "higher",
  top20_overlap_fraction = "higher",
  replicated_downstream_fdr_count = "higher",
  downstream_fdr_hoa2 = "higher",
  downstream_fdr_hoa3 = "higher"
)
calibration_rows <- list()
calibration_index <- 0L
for (profile in names(profiles)) {
  one_profile <- summary_table[summary_table$profile == profile, , drop = FALSE]
  target_summary <- one_profile[one_profile$role == "target", , drop = FALSE]
  controls <- one_profile[one_profile$role == "matched comparator", , drop = FALSE]
  for (metric in names(metric_directions)) {
    target_value <- target_summary[[metric]]
    comparator_values <- controls[[metric]]
    calibration_index <- calibration_index + 1L
    calibration_rows[[calibration_index]] <- data.frame(
      profile = profile,
      metric = metric,
      direction = metric_directions[[metric]],
      sqstm1_value = target_value,
      comparator_median = stats::median(comparator_values),
      comparator_min = min(comparator_values),
      comparator_max = max(comparator_values),
      sqstm1_percentile = 100 * mean(comparator_values <= target_value),
      exact_empirical_p = (1 + sum(comparator_values >= target_value)) /
        (1 + length(comparator_values)),
      stringsAsFactors = FALSE
    )
  }
}
calibration <- do.call(rbind, calibration_rows)
utils::write.csv(
  calibration,
  file.path(output_dir, "sqstm1_matched_control_calibration.csv"),
  row.names = FALSE
)

selection_rationale <- data.frame(
  step = 1:6,
  criterion = c(
    "Prespecified biology",
    "Directional disease-oriented contrast",
    "Independent-donor contrast",
    "Control-EC detectability",
    "Network eligibility",
    "Interpretive boundary"
  ),
  observed_evidence = c(
    "SQSTM1 belonged to the prespecified selective-clearance/mitophagy panel.",
    "Estimated SONFH-versus-HOA EC log2 fold change = -1.265 (descriptive cross-cohort contrast).",
    "Estimated Liao ARCO 3A-versus-HOA EC log2 fold change = -0.450; direction matched but FDR was not significant.",
    "Detected in 70.8% of HOA2 ECs and 72.2% of HOA3 ECs.",
    "HOA2 (1,043 ECs) and HOA3 (760 ECs) were modeled separately; HOA1 (26 ECs) was excluded.",
    "Selection supports a loss-hypothesis perturbation only; it does not establish SQSTM1 causality."
  ),
  decision = c(
    "retain", "retain as descriptive support", "retain as directional support",
    "eligible", "two-donor analysis", "exploratory target"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  selection_rationale,
  file.path(output_dir, "sqstm1_selection_rationale.csv"),
  row.names = FALSE
)

provenance <- list(
  analysis = "matched-comparator virtual knockouts using frozen official-R WT networks",
  target = target,
  comparator_genes = matched$gene,
  closest_five_by_standardized_distance = matched$gene[matched$closest_five],
  profiles = profiles,
  donors = donors,
  selection = as.list(stats::setNames(design$specification, design$item)),
  seed = 1,
  n_cores = cli$cores,
  validation = validation,
  software = list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    scTenifoldKnk = as.character(utils::packageVersion("scTenifoldKnk")),
    scTenifoldNet = as.character(utils::packageVersion("scTenifoldNet")),
    Matrix = as.character(utils::packageVersion("Matrix")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))
  ),
  guardrail = paste(
    "Comparator genes are matched computational perturbations, not",
    "experimentally established negative controls."
  )
)
jsonlite::write_json(
  provenance,
  file.path(output_dir, "matched_control_vko_provenance.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

summary_lines <- c(
  "Matched-comparator virtual-knockout audit",
  paste0("Comparators (n=", nrow(matched), "): ", paste(matched$gene, collapse = ", ")),
  paste0("Closest five: ", paste(matched$gene[matched$closest_five], collapse = ", ")),
  apply(calibration, 1L, function(x) {
    paste0(
      x[["profile"]], " | ", x[["metric"]],
      ": SQSTM1=", signif(as.numeric(x[["sqstm1_value"]]), 4),
      ", empirical p=", signif(as.numeric(x[["exact_empirical_p"]]), 4)
    )
  })
)
writeLines(summary_lines, file.path(output_dir, "RUN_SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
