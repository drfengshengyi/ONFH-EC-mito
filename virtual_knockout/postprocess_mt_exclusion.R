#!/usr/bin/env Rscript

# Compare original and complete-refit SQSTM1 virtual knockouts after excluding
# all frozen-feature symbols beginning MT-. Comparisons use the same 294
# shared, non-target, nuclear-encoded genes and rank percentiles so the
# 300-versus-295 model dimensions are not compared on their raw rank scales.

options(stringsAsFactors = FALSE, warn = 1)

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

rank_on_universe <- function(frame, universe) {
  frame <- frame[frame$Gene %in% universe, , drop = FALSE]
  frame <- frame[order(frame$p_value, frame$Gene), , drop = FALSE]
  frame$common_rank <- seq_len(nrow(frame))
  frame$rank_percentile <- if (nrow(frame) == 1L) {
    0
  } else {
    (frame$common_rank - 1) / (nrow(frame) - 1)
  }
  frame$common_adjusted_p_value <- stats::p.adjust(frame$p_value, method = "BH")
  frame
}

unique_sample_matrix <- function(pool, k, n_permutations) {
  if (k < 1L || k > length(pool)) stop("Invalid stratum sample size", call. = FALSE)
  if (k == length(pool)) {
    return(matrix(rep(pool, each = n_permutations), nrow = n_permutations))
  }
  sampled <- matrix(
    sample(pool, n_permutations * k, replace = TRUE),
    nrow = n_permutations,
    ncol = k
  )
  if (k > 1L) {
    duplicate_rows <- function(x) {
      pairs <- utils::combn(seq_len(ncol(x)), 2L)
      duplicate <- rep(FALSE, nrow(x))
      for (pair_index in seq_len(ncol(pairs))) {
        duplicate <- duplicate |
          x[, pairs[1L, pair_index]] == x[, pairs[2L, pair_index]]
      }
      duplicate
    }
    duplicate <- duplicate_rows(sampled)
    while (any(duplicate)) {
      n_duplicate <- sum(duplicate)
      sampled[duplicate, ] <- matrix(
        sample(pool, n_duplicate * k, replace = TRUE),
        nrow = n_duplicate,
        ncol = k
      )
      duplicate <- duplicate_rows(sampled)
    }
  }
  sampled
}

make_null_indices <- function(universe, strata, genes, n_permutations, seed) {
  set.seed(seed)
  gene_indices <- match(genes, universe)
  counts <- table(strata[gene_indices])
  matrices <- list()
  matrix_index <- 0L
  for (stratum_name in names(counts)) {
    k <- unname(counts[[stratum_name]])
    pool <- which(strata == as.integer(stratum_name))
    matrix_index <- matrix_index + 1L
    matrices[[matrix_index]] <- unique_sample_matrix(pool, k, n_permutations)
  }
  do.call(cbind, matrices)
}

analysis_dir <- script_location()
project_root <- normalizePath(file.path(analysis_dir, ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
target <- "SQSTM1"
donors <- c("hoa2", "hoa3")
n_permutations <- 20000L
seed <- 20260829L

profiles <- list(
  manuscript = list(
    original = "official_r_vko_manuscript",
    mt_excluded = "official_r_vko_manuscript_no_mt_encoded"
  ),
  official_default = list(
    original = "official_r_vko_official_default",
    mt_excluded = "official_r_vko_official_default_no_mt_encoded"
  )
)

features <- utils::read.csv(
  file.path(analysis_dir, "vko_selected_features_v5.csv"),
  check.names = FALSE
)
mt_encoded <- sort(features$gene[grepl("^MT-", features$gene)])
expected_mt <- sort(c("MT-ATP6", "MT-CO1", "MT-CO2", "MT-ND1", "MT-ND4"))
if (!identical(mt_encoded, expected_mt)) {
  stop("Unexpected frozen MT-encoded feature set", call. = FALSE)
}
universe <- sort(setdiff(features$gene, c(target, mt_encoded)))
if (length(universe) != 294L) stop("Expected 294 shared downstream genes", call. = FALSE)

read_result <- function(directory, donor) {
  path <- file.path(
    results_dir,
    directory,
    paste0("vko_sqstm1_", donor, "_official_r.csv")
  )
  if (!file.exists(path)) stop("Missing analysis output: ", path, call. = FALSE)
  frame <- utils::read.csv(path, check.names = FALSE)
  required <- c("Gene", "rank", "p_value", "adjusted_p_value")
  if (!all(required %in% names(frame))) stop("Malformed result: ", path, call. = FALSE)
  rank_on_universe(frame, universe)
}

all_results <- list()
for (profile in names(profiles)) {
  all_results[[profile]] <- list()
  for (variant in names(profiles[[profile]])) {
    all_results[[profile]][[variant]] <- setNames(
      lapply(donors, function(donor) {
        read_result(profiles[[profile]][[variant]], donor)
      }),
      donors
    )
  }
}

gene_rows <- list()
gene_index <- 0L
summary_rows <- list()
summary_index <- 0L
for (profile in names(all_results)) {
  for (donor in donors) {
    original <- all_results[[profile]]$original[[donor]]
    excluded <- all_results[[profile]]$mt_excluded[[donor]]
    joined <- merge(
      original[, c(
        "Gene", "common_rank", "rank_percentile", "p_value",
        "common_adjusted_p_value"
      )],
      excluded[, c(
        "Gene", "common_rank", "rank_percentile", "p_value",
        "common_adjusted_p_value"
      )],
      by = "Gene",
      suffixes = c("_original", "_mt_excluded")
    )
    joined$profile <- profile
    joined$donor <- donor
    joined$rank_percentile_shift <-
      joined$rank_percentile_mt_excluded - joined$rank_percentile_original
    joined <- joined[, c("profile", "donor", setdiff(names(joined), c("profile", "donor")))]
    gene_index <- gene_index + 1L
    gene_rows[[gene_index]] <- joined

    rho <- suppressWarnings(stats::cor.test(
      joined$common_rank_original,
      joined$common_rank_mt_excluded,
      method = "spearman",
      exact = FALSE
    ))
    top20_original <- joined$Gene[order(joined$common_rank_original)][1:20]
    top20_excluded <- joined$Gene[order(joined$common_rank_mt_excluded)][1:20]
    top50_original <- joined$Gene[order(joined$common_rank_original)][1:50]
    top50_excluded <- joined$Gene[order(joined$common_rank_mt_excluded)][1:50]
    significant_original <- joined$Gene[joined$common_adjusted_p_value_original < 0.05]
    significant_excluded <- joined$Gene[joined$common_adjusted_p_value_mt_excluded < 0.05]
    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- data.frame(
      profile = profile,
      donor = donor,
      n_shared_downstream_genes = nrow(joined),
      spearman_rho = unname(rho$estimate),
      spearman_p_value = rho$p.value,
      top20_overlap_count = length(intersect(top20_original, top20_excluded)),
      top20_jaccard = length(intersect(top20_original, top20_excluded)) /
        length(union(top20_original, top20_excluded)),
      top50_overlap_count = length(intersect(top50_original, top50_excluded)),
      top50_jaccard = length(intersect(top50_original, top50_excluded)) /
        length(union(top50_original, top50_excluded)),
      common_bh_fdr_original_count = length(significant_original),
      common_bh_fdr_mt_excluded_count = length(significant_excluded),
      common_bh_fdr_overlap_count = length(intersect(
        significant_original,
        significant_excluded
      )),
      stringsAsFactors = FALSE
    )
  }
}

gene_comparison <- do.call(rbind, gene_rows)
summary_table <- do.call(rbind, summary_rows)
utils::write.csv(
  gene_comparison,
  file.path(results_dir, "official_r_vko_no_mt_gene_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_table,
  file.path(results_dir, "official_r_vko_no_mt_summary.csv"),
  row.names = FALSE
)

# Freeze the cross-donor common-nuclear audit used in the Results text and
# supplementary Table S10g. This comparison is distinct from summary_table:
# it compares HOA2 with HOA3 after the mtDNA-feature-excluded refit rather than
# comparing original and refitted models within each donor.
cross_donor_rows <- lapply(names(all_results), function(profile) {
  hoa2 <- all_results[[profile]]$mt_excluded$hoa2
  hoa3 <- all_results[[profile]]$mt_excluded$hoa3
  joined <- merge(
    hoa2[, c("Gene", "common_rank", "common_adjusted_p_value")],
    hoa3[, c("Gene", "common_rank", "common_adjusted_p_value")],
    by = "Gene",
    suffixes = c("_hoa2", "_hoa3")
  )
  joined <- joined[order(joined$Gene), , drop = FALSE]
  rho <- suppressWarnings(stats::cor.test(
    joined$common_rank_hoa2,
    joined$common_rank_hoa3,
    method = "spearman",
    exact = FALSE
  ))
  top20_hoa2 <- joined$Gene[order(joined$common_rank_hoa2)][1:20]
  top20_hoa3 <- joined$Gene[order(joined$common_rank_hoa3)][1:20]
  significant_hoa2 <- sort(joined$Gene[joined$common_adjusted_p_value_hoa2 < 0.05])
  significant_hoa3 <- sort(joined$Gene[joined$common_adjusted_p_value_hoa3 < 0.05])
  replicated <- intersect(significant_hoa2, significant_hoa3)
  data.frame(
    profile = profile,
    model_variant = "mt_excluded",
    donor_pair = "hoa2_vs_hoa3",
    n_common_nuclear_genes = nrow(joined),
    rank_metric = "common_rank",
    spearman_rho = unname(rho$estimate),
    spearman_p_value = rho$p.value,
    top20_overlap_count = length(intersect(top20_hoa2, top20_hoa3)),
    top20_jaccard = length(intersect(top20_hoa2, top20_hoa3)) /
      length(union(top20_hoa2, top20_hoa3)),
    common_bh_fdr_hoa2_count = length(significant_hoa2),
    common_bh_fdr_hoa3_count = length(significant_hoa3),
    common_bh_fdr_replicated_count = length(replicated),
    common_bh_fdr_hoa2_genes = paste(significant_hoa2, collapse = ";"),
    common_bh_fdr_hoa3_genes = paste(significant_hoa3, collapse = ";"),
    common_bh_fdr_replicated_genes = paste(replicated, collapse = ";"),
    stringsAsFactors = FALSE
  )
})
cross_donor_audit <- do.call(rbind, cross_donor_rows)
utils::write.csv(
  cross_donor_audit,
  file.path(results_dir, "official_r_vko_no_mt_cross_donor_audit.csv"),
  row.names = FALSE
)

# Record the exact downstream HOA3 FDR-positive identities under the primary
# mtDNA-feature-excluded refit. Both BH families are retained: the official
# 295-row output (including SQSTM1) and the 294-gene non-target common universe.
hoa3_primary <- all_results$manuscript$mt_excluded$hoa3
hoa3_fdr_audit <- hoa3_primary[
  hoa3_primary$common_adjusted_p_value < 0.05,
  c(
    "Gene", "rank", "common_rank", "rank_percentile", "Distance", "Z", "FC",
    "p_value", "adjusted_p_value", "common_adjusted_p_value"
  ),
  drop = FALSE
]
hoa3_fdr_audit <- hoa3_fdr_audit[order(hoa3_fdr_audit$common_rank), , drop = FALSE]
names(hoa3_fdr_audit) <- c(
  "gene", "official_rank", "common_rank", "common_rank_percentile", "distance",
  "z_score", "fold_change", "raw_p_value", "official_bh_fdr_295_family",
  "common_nuclear_bh_fdr_294_family"
)
hoa3_fdr_audit <- data.frame(
  profile = "manuscript",
  model_variant = "mt_excluded",
  donor = "hoa3",
  hoa3_fdr_audit,
  common_nuclear_fdr_0_05 = hoa3_fdr_audit$common_nuclear_bh_fdr_294_family < 0.05,
  stringsAsFactors = FALSE
)
utils::write.csv(
  hoa3_fdr_audit,
  file.path(results_dir, "official_r_vko_no_mt_hoa3_fdr_audit.csv"),
  row.names = FALSE
)

matching_columns <- c(
  "prevalence_hoa2", "prevalence_hoa3",
  "mean_log1p_cpm_hoa2", "mean_log1p_cpm_hoa3"
)
feature_meta <- features[match(universe, features$gene), c("gene", matching_columns)]
percentile_columns <- lapply(matching_columns, function(column) {
  rank(feature_meta[[column]], ties.method = "average") / nrow(feature_meta)
})
matching_score <- rowMeans(do.call(cbind, percentile_columns))
ordered_score <- rank(matching_score, ties.method = "first")
strata <- pmin(5L, ceiling(5 * ordered_score / length(ordered_score)))

gene_sets <- jsonlite::read_json(
  file.path(project_root, "analysis", "genesets_final.json"),
  simplifyVector = TRUE
)
pathway_names <- names(gene_sets)
null_indices <- list()
pathway_members <- list()
for (pathway_index in seq_along(pathway_names)) {
  pathway <- pathway_names[[pathway_index]]
  genes <- sort(intersect(as.character(gene_sets[[pathway]]), universe))
  if (length(genes) < 5L) next
  pathway_members[[pathway]] <- genes
  null_indices[[pathway]] <- make_null_indices(
    universe = universe,
    strata = strata,
    genes = genes,
    n_permutations = n_permutations,
    seed = seed + pathway_index
  )
}

pathway_rows <- list()
pathway_row_index <- 0L
for (profile in names(all_results)) {
  for (variant in names(all_results[[profile]])) {
    for (donor in donors) {
      frame <- all_results[[profile]][[variant]][[donor]]
      rank_values <- frame$common_rank[match(universe, frame$Gene)]
      for (pathway in names(pathway_members)) {
        genes <- pathway_members[[pathway]]
        observed <- mean(rank_values[match(genes, universe)])
        sampled <- null_indices[[pathway]]
        permuted <- rowMeans(matrix(
          rank_values[sampled],
          nrow = nrow(sampled),
          ncol = ncol(sampled)
        ))
        empirical_p <- (1 + sum(permuted <= observed)) / (n_permutations + 1)
        pathway_row_index <- pathway_row_index + 1L
        pathway_rows[[pathway_row_index]] <- data.frame(
          profile = profile,
          model_variant = variant,
          donor = donor,
          pathway = pathway,
          n_nuclear_genes = length(genes),
          nuclear_genes = paste(genes, collapse = ";"),
          observed_mean_common_rank = observed,
          empirical_p = empirical_p,
          null_model = paste(
            "20,000 sets matched by quintiles of control detection",
            "prevalence and mean log1p CPM on the shared nuclear universe"
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
pathways <- do.call(rbind, pathway_rows)
pathways$empirical_fdr <- NA_real_
for (profile in unique(pathways$profile)) {
  for (variant in unique(pathways$model_variant)) {
    for (donor in donors) {
      take <- pathways$profile == profile &
        pathways$model_variant == variant &
        pathways$donor == donor
      pathways$empirical_fdr[take] <- stats::p.adjust(
        pathways$empirical_p[take],
        method = "BH"
      )
    }
  }
}

pathway_summary <- list()
pathway_summary_index <- 0L
for (profile in unique(pathways$profile)) {
  for (variant in unique(pathways$model_variant)) {
    subset_frame <- pathways[
      pathways$profile == profile & pathways$model_variant == variant,
      ,
      drop = FALSE
    ]
    for (pathway in unique(subset_frame$pathway)) {
      one <- subset_frame[subset_frame$pathway == pathway, , drop = FALSE]
      p_values <- one$empirical_p[match(donors, one$donor)]
      pathway_summary_index <- pathway_summary_index + 1L
      pathway_summary[[pathway_summary_index]] <- data.frame(
        profile = profile,
        model_variant = variant,
        pathway = pathway,
        n_nuclear_genes = one$n_nuclear_genes[[1L]],
        replicated_nominal_0_05 = all(p_values < 0.05),
        fisher_p = stats::pchisq(-2 * sum(log(p_values)), df = 4, lower.tail = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }
}
pathway_summary <- do.call(rbind, pathway_summary)
pathway_summary$fisher_fdr <- NA_real_
for (profile in unique(pathway_summary$profile)) {
  for (variant in unique(pathway_summary$model_variant)) {
    take <- pathway_summary$profile == profile & pathway_summary$model_variant == variant
    pathway_summary$fisher_fdr[take] <- stats::p.adjust(
      pathway_summary$fisher_p[take],
      method = "BH"
    )
  }
}

utils::write.csv(
  pathways,
  file.path(results_dir, "official_r_vko_no_mt_pathway_by_donor.csv"),
  row.names = FALSE
)
utils::write.csv(
  pathway_summary,
  file.path(results_dir, "official_r_vko_no_mt_pathway_summary.csv"),
  row.names = FALSE
)

manifest <- data.frame(
  item = c(
    "sensitivity design", "original feature count", "excluded feature count",
    "excluded symbols", "refit feature count", "comparison universe",
    "profiles", "donors", "rank scale", "pathway null"
  ),
  value = c(
    "Complete donor-specific network refit; no replacement genes",
    "300", "5", paste(mt_encoded, collapse = ";"), "295",
    "294 shared, non-target, nuclear-encoded genes",
    "manuscript;official_default", "hoa2;hoa3",
    "Within-universe rank percentile",
    "20,000 expression/prevalence-matched sets"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  file.path(results_dir, "vko_mt_encoded_exclusion_manifest.csv"),
  row.names = FALSE
)

writeLines(
  sub("[[:space:]]+$", "", capture.output(utils::sessionInfo())),
  file.path(results_dir, "official_r_vko_no_mt_sessionInfo.txt")
)
cat("Wrote mtDNA-feature exclusion comparisons for ", nrow(summary_table), " donor-profile pairs.\n", sep = "")
