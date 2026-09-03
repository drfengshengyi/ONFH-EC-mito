#!/usr/bin/env Rscript

# Hallmark leading-edge and leave-one-participant-out sensitivity analysis.
# The fgsea P values remain conditional on the competitive gene-set null. The
# six omitted-participant fits quantify rank stability and must not be read as
# six independent cohorts or as participant-resampled confirmatory evidence.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
root <- dirname(dirname(script_file))
out <- Sys.getenv(
  "ONFH_FGSEA_STABILITY_DIR",
  unset = file.path(root, "results", "participant_fgsea_stability")
)
ranks_dir <- file.path(out, "de_ranks")
if (!dir.exists(ranks_dir)) stop("Missing rank directory: ", ranks_dir)

msig <- msigdbr(species = "Homo sapiens")
hallmark <- msig[msig$gs_collection == "H", ]
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

run_one <- function(path) {
  de <- fread(path)
  required <- c("gene", "stat")
  if (length(setdiff(required, names(de)))) {
    stop("Rank file missing gene/stat columns: ", path)
  }
  de <- de[!is.na(stat) & !is.na(gene) & gene != ""]
  de <- de[!duplicated(gene)]
  ranks <- de$stat
  names(ranks) <- de$gene
  ranks <- sort(ranks, decreasing = TRUE)
  # Reuse the archived primary analysis seed so the full-fit values remain
  # byte-for-byte comparable to the manuscript's frozen Hallmark table.
  set.seed(42)
  fg <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )
  fg <- as.data.frame(fg[order(fg$padj), ])
  fg$leading_edge_count <- lengths(fg$leadingEdge)
  fg$leading_edge_genes <- vapply(
    fg$leadingEdge,
    function(x) paste(x, collapse = ";"),
    character(1)
  )
  fg$leadingEdge <- NULL
  fg
}

rank_files <- list.files(ranks_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
if (!length(rank_files)) stop("No rank files in ", ranks_dir)
all_results <- lapply(rank_files, function(path) {
  fit_id <- sub("\\.csv\\.gz$", "", basename(path))
  x <- run_one(path)
  x$fit_id <- fit_id
  x
})
all_results <- bind_rows(all_results) %>%
  select(fit_id, pathway, NES, pval, padj, size, leading_edge_count, leading_edge_genes)
write.csv(all_results, file.path(out, "fgsea_hallmark_all_fits.csv"), row.names = FALSE)

full <- all_results %>% filter(fit_id == "full") %>% arrange(padj)
selected <- c(
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
selected_results <- all_results %>% filter(pathway %in% selected)
write.csv(
  selected_results,
  file.path(out, "fgsea_selected_pathways_by_omission.csv"),
  row.names = FALSE
)

lopo_summary <- selected_results %>%
  filter(fit_id != "full") %>%
  group_by(pathway) %>%
  summarise(
    n_omissions = n(),
    n_positive_nes = sum(NES > 0, na.rm = TRUE),
    n_fdr_below_0_05 = sum(padj < 0.05, na.rm = TRUE),
    nes_min = min(NES, na.rm = TRUE),
    nes_median = median(NES, na.rm = TRUE),
    nes_max = max(NES, na.rm = TRUE),
    padj_min = min(padj, na.rm = TRUE),
    padj_median = median(padj, na.rm = TRUE),
    padj_max = max(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    full %>%
      filter(pathway %in% selected) %>%
      select(pathway, full_nes = NES, full_padj = padj, full_leading_edge_count = leading_edge_count),
    by = "pathway"
  )
write.csv(lopo_summary, file.path(out, "fgsea_lopo_summary.csv"), row.names = FALSE)

full_ranks <- fread(file.path(ranks_dir, "full.csv.gz")) %>%
  select(gene, stat)
leading_edge_top20 <- full %>%
  filter(pathway %in% selected) %>%
  select(pathway, leading_edge_genes) %>%
  tidyr::separate_rows(leading_edge_genes, sep = ";") %>%
  rename(gene = leading_edge_genes) %>%
  left_join(full_ranks, by = "gene") %>%
  group_by(pathway) %>%
  arrange(desc(stat), .by_group = TRUE) %>%
  mutate(driver_rank = row_number()) %>%
  filter(driver_rank <= 20) %>%
  ungroup() %>%
  select(pathway, driver_rank, gene, stat)
write.csv(
  leading_edge_top20,
  file.path(out, "fgsea_leading_edge_top20.csv"),
  row.names = FALSE
)

pairwise_overlap <- function(sets, type) {
  pairs <- combn(names(sets), 2, simplify = FALSE)
  bind_rows(lapply(pairs, function(pair) {
    a <- unique(sets[[pair[[1]]]])
    b <- unique(sets[[pair[[2]]]])
    data.frame(
      overlap_type = type,
      pathway_a = pair[[1]],
      pathway_b = pair[[2]],
      n_a = length(a),
      n_b = length(b),
      n_intersection = length(intersect(a, b)),
      jaccard = length(intersect(a, b)) / length(union(a, b)),
      stringsAsFactors = FALSE
    )
  }))
}
gene_sets_selected <- pathways[selected]
leading_sets_selected <- setNames(
  lapply(selected, function(p) {
    value <- full$leading_edge_genes[match(p, full$pathway)]
    if (is.na(value)) character(0) else strsplit(value, ";", fixed = TRUE)[[1]]
  }),
  selected
)
overlap <- bind_rows(
  pairwise_overlap(gene_sets_selected, "Hallmark gene-set membership"),
  pairwise_overlap(leading_sets_selected, "Full-fit leading edge")
)
write.csv(overlap, file.path(out, "fgsea_selected_pathway_overlap.csv"), row.names = FALSE)

method <- data.frame(
  ranking_statistic = "participant-pseudobulk DESeq2 Wald statistic",
  full_design = "three ONFH ARCO 3A versus three HOA participants",
  sensitivity_design = "six leave-one-participant-out refits",
  enrichment_null = "competitive gene-set null; participants are not permuted by fgsea",
  interpretive_scope = paste(
    "leave-one-participant-out sensitivity to a single influential participant;",
    "not independent cohort validation"
  ),
  fgsea_version = as.character(packageVersion("fgsea")),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  stringsAsFactors = FALSE
)
write.csv(method, file.path(out, "fgsea_lopo_method_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat(normalizePath(out, winslash = "/"), "\n")
