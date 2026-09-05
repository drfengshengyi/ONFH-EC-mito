#!/usr/bin/env Rscript

# Genes-revision redraw of Figures 2, 3, 4 and 7.
# Core conclusion: participant-aware analyses define reproducibility limits
# and testable gene/pathway hypotheses rather than confirmatory mechanisms.
# Backend: R only. All labels and statistics are read from versioned tables.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(patchwork)
  library(pROC)
  library(ragg)
  library(tidyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
ROOT <- dirname(dirname(script_file))
INPUT <- file.path(ROOT, "results", "figure_inputs")
ROBUST <- file.path(ROOT, "results", "participant_fgsea_stability")
OUT <- file.path(ROOT, "figures", "final")
QA <- file.path(ROOT, "qa")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(QA, recursive = TRUE, showWarnings = FALSE)

FONT <- "sans"
GROUP_ORDER <- c("Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH")
GROUP_LABELS <- c("Healthy", "HOA", "FNF", "ARCO 3A", "ARCO 4", "SONFH\n(libraries)")
GROUP_COLORS <- c(
  Healthy = "#2C7FB8", HOA = "#56B4E9", FNF = "#1B9E77",
  ONFH_3A = "#E6A51A", ONFH_4 = "#D95F0E", SONFH = "#CC79A7"
)
SUBTYPE_ORDER <- c(
  "Lymphatic", "Type H / EMCN-KDR", "Type R / bone-remodeling", "Venous / ACKR1"
)
SUBTYPE_COLORS <- c(
  "Lymphatic" = "#CC79A7", "Type H / EMCN-KDR" = "#0072B2",
  "Type R / bone-remodeling" = "#009E73", "Venous / ACKR1" = "#E69F00"
)
SUBTYPE_COLUMNS <- c(
  "Lymphatic" = "lymphatic", "Type H / EMCN-KDR" = "typeH_EMCN_KDR",
  "Type R / bone-remodeling" = "typeR_bone_remodel", "Venous / ACKR1" = "venous_ACKR1"
)
MODULE_LABELS <- c(
  Mito_fission = "Mito fission", Mito_fusion = "Mito fusion",
  Mitophagy_core = "Mitophagy core", cGAS_STING = "cGAS-STING",
  EC_inflammation = "EC inflammation", YAP_mTOR = "YAP-mTOR"
)

theme_pub <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.7, colour = "#4B5563", margin = margin(b = 3)),
      plot.tag = element_text(size = 8.8, face = "bold"),
      axis.title = element_text(size = 6.4),
      axis.text = element_text(size = 5.8, colour = "#222222"),
      axis.line = element_line(linewidth = 0.3, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.1, face = "bold"),
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.4),
      legend.key.height = unit(2.6, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "#E5E7EB"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_pub <- function(plot, name, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::cairo_pdf(
    file.path(OUT, paste0(name, ".pdf")), width = width_in, height = height_in,
    family = FONT, onefile = TRUE, bg = "white"
  )
  print(plot)
  dev.off()
  ragg::agg_png(
    file.path(OUT, paste0(name, ".png")), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  dev.off()
}

read_modules <- function() {
  x <- fread(file.path(INPUT, "module_scores_by_library_v4.csv"))
  x$group <- factor(x$group, levels = GROUP_ORDER)
  x
}

sampling_panel <- function(data, columns, ncol, y_label, free_y = TRUE) {
  labels <- unname(MODULE_LABELS[columns])
  long <- data %>%
    select(sample, group, all_of(columns)) %>%
    pivot_longer(all_of(columns), names_to = "measure", values_to = "value") %>%
    mutate(
      measure = factor(measure, levels = columns, labels = labels),
      group = factor(group, levels = GROUP_ORDER)
    )
  summary <- long %>%
    group_by(measure, group) %>%
    summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), .groups = "drop")
  ggplot(long, aes(group, value, colour = group)) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0, seed = 20260903),
      size = 1.15, alpha = 0.88
    ) +
    geom_errorbar(
      data = summary,
      inherit.aes = FALSE,
      mapping = aes(x = group, ymin = mean - sd, ymax = mean + sd),
      width = 0.15, linewidth = 0.30, colour = "#333333"
    ) +
    geom_point(
      data = summary,
      aes(x = group, y = mean),
      inherit.aes = FALSE,
      shape = 21, size = 1.55, fill = "white", colour = "#111111", stroke = 0.35
    ) +
    facet_wrap(~measure, ncol = ncol, scales = if (free_y) "free_y" else "fixed") +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    scale_x_discrete(labels = setNames(GROUP_LABELS, GROUP_ORDER)) +
    labs(x = NULL, y = y_label) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))
}

make_figure2 <- function() {
  umap <- fread(file.path(INPUT, "figure2_umap.csv.gz")) %>%
    mutate(subtype = factor(subtype, levels = SUBTYPE_ORDER))
  dot <- fread(file.path(INPUT, "figure2_marker_dotplot.csv"))
  dot$gene <- factor(dot$gene, levels = unique(dot$gene))
  dot$subtype <- factor(dot$subtype, levels = rev(SUBTYPE_ORDER))
  composition <- fread(file.path(INPUT, "ec_subtype_composition_v4.csv"))
  composition$group <- factor(composition$group, levels = GROUP_ORDER)
  confidence <- fread(file.path(ROOT, "results", "ec_annotation_confidence_v4.csv")) %>%
    mutate(
      cluster = factor(cluster, levels = rev(cluster[order(top_two_margin)])),
      subtype_label = recode(
        assigned_subtype,
        lymphatic = "Lymphatic", typeH_EMCN_KDR = "Type H / EMCN-KDR",
        typeR_bone_remodel = "Type R / bone-remodeling", venous_ACKR1 = "Venous / ACKR1"
      )
    )

  p_a <- ggplot(umap, aes(UMAP1, UMAP2, colour = subtype)) +
    geom_point(size = 0.16, alpha = 0.68) +
    scale_colour_manual(values = SUBTYPE_COLORS, drop = FALSE) +
    coord_equal() +
    labs(title = "Marker-panel-defined endothelial states", x = NULL, y = NULL, colour = NULL) +
    theme_pub() +
    theme(
      axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )

  p_b <- ggplot(dot, aes(gene, subtype)) +
    geom_point(aes(size = fraction, colour = scaled_mean), alpha = 0.92) +
    scale_size_continuous(range = c(0.35, 3.2), labels = scales::percent_format(accuracy = 1)) +
    scale_colour_gradient(low = "#E8F1F7", high = "#08519C", limits = c(0, 1)) +
    labs(
      title = "Marker support for endothelial-state labels",
      x = NULL, y = NULL, size = "Cells expressing", colour = "Scaled mean"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.x = element_line(linewidth = 0.18, colour = "#F0F0F0"),
      legend.position = "right"
    )

  comp_long <- composition %>%
    select(sample, group, all_of(unname(SUBTYPE_COLUMNS))) %>%
    pivot_longer(all_of(unname(SUBTYPE_COLUMNS)), names_to = "measure", values_to = "value") %>%
    mutate(
      measure = factor(measure, levels = unname(SUBTYPE_COLUMNS), labels = names(SUBTYPE_COLUMNS)),
      group = factor(group, levels = GROUP_ORDER)
    )
  comp_summary <- comp_long %>%
    group_by(measure, group) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  p_c <- ggplot(comp_long, aes(group, value, colour = group)) +
    geom_point(
      position = position_jitter(width = 0.10, seed = 20260903), size = 1.1, alpha = 0.88
    ) +
    geom_errorbar(
      data = comp_summary,
      mapping = aes(x = group, ymin = pmax(0, mean - sd), ymax = mean + sd),
      inherit.aes = FALSE, width = 0.14, linewidth = 0.3, colour = "#333333"
    ) +
    geom_point(
      data = comp_summary, aes(x = group, y = mean), inherit.aes = FALSE,
      shape = 21, size = 1.5, fill = "white", colour = "#111111", stroke = 0.35
    ) +
    facet_wrap(~measure, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    scale_x_discrete(labels = setNames(GROUP_LABELS, GROUP_ORDER)) +
    labs(title = "Subtype composition by sampling unit", x = NULL, y = "Endothelial cells (%)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 28, hjust = 1))

  p_d <- ggplot(confidence, aes(top_two_margin, cluster, fill = subtype_label)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.20) +
    geom_vline(xintercept = 0.15, linetype = "dashed", colour = "#B23A33", linewidth = 0.45) +
    annotate(
      "text", x = 0.15, y = Inf, label = "Prespecified margin = 0.15",
      hjust = -0.06, vjust = 1.4, size = 1.75, colour = "#B23A33", family = FONT
    ) +
    scale_fill_manual(values = SUBTYPE_COLORS, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = "Annotation sensitivity",
      subtitle = "Seven of nine clusters had a top-two panel-score margin <0.15",
      x = "Top-panel minus second-panel score", y = "Leiden cluster", fill = "Assigned state"
    ) +
    theme_pub()

  figure <- ((p_a | p_b) + plot_layout(widths = c(0.90, 1.35))) /
    ((p_c | p_d) + plot_layout(widths = c(1.45, 0.80))) +
    plot_annotation(tag_levels = "A") &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_pub(figure, "Figure2", 183, 138)
}

make_figure3 <- function() {
  de <- fread(file.path(INPUT, "de_ec_SONFH_vs_HOA_descriptive_v4.csv")) %>%
    rename(gene = V1) %>%
    mutate(abundance = log10(pmax(baseMean, 0) + 1))
  effects <- fread(file.path(INPUT, "fig3_key_gene_effects_v4.csv"))
  names(effects) <- c("gene", "Liao ARCO 3A vs HOA", "SONFH vs HOA")
  effects_long <- effects %>%
    pivot_longer(-gene, names_to = "contrast", values_to = "log2_fold_change") %>%
    mutate(contrast = factor(contrast, levels = c("Liao ARCO 3A vs HOA", "SONFH vs HOA")))
  targets <- effects$gene

  p_a <- ggplot(de, aes(abundance, log2FoldChange)) +
    geom_point(size = 0.36, colour = "#9FB3C1", alpha = 0.22) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#333333") +
    geom_point(data = filter(de, gene %in% targets), size = 1.45, colour = "#D9472F") +
    geom_text_repel(
      data = filter(de, gene %in% targets), aes(label = gene),
      seed = 20260903, size = 1.85, family = FONT, box.padding = 0.28,
      point.padding = 0.18, min.segment.length = 0, max.overlaps = Inf,
      segment.size = 0.20, segment.colour = "#737373"
    ) +
    labs(
      title = "Descriptive SONFH library-level endothelial effects",
      subtitle = "Cross-cohort effect sizes only; participant mapping was unavailable",
      x = "Mean endothelial pseudobulk abundance, log10(baseMean + 1)",
      y = "Descriptive log2 fold change"
    ) +
    theme_pub()

  p_b <- ggplot(effects_long, aes(contrast, factor(gene, levels = rev(effects$gene)))) +
    geom_tile(aes(fill = log2_fold_change), colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", log2_fold_change)), size = 1.65, family = FONT) +
    scale_fill_gradient2(low = "#3B7EA1", mid = "white", high = "#D95F59", midpoint = 0) +
    labs(
      title = "Prespecified gene-level effects",
      x = NULL, y = NULL, fill = "log2 fold change"
    ) +
    scale_x_discrete(labels = c(
      "Liao ARCO 3A vs HOA" = "Liao ARCO 3A vs HOA\n(independent participants)",
      "SONFH vs HOA" = "SONFH vs HOA\n(descriptive libraries)"
    )) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

  p_c <- sampling_panel(
    read_modules(), c("Mito_fission", "Mito_fusion", "Mitophagy_core"),
    ncol = 3, y_label = "Mean module score"
  ) +
    labs(title = "Mitochondrial module scores by sampling unit")

  figure <- ((p_a | p_b) + plot_layout(widths = c(1.42, 0.78))) / p_c +
    plot_layout(heights = c(1.12, 0.88)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_pub(figure, "Figure3", 183, 127)
}

format_fdr <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x))
}

make_figure4 <- function() {
  gsea <- fread(file.path(INPUT, "gsea_ONFH3A_vs_HOA_H.csv")) %>%
    filter(!is.na(NES), !is.na(padj), padj < 0.05) %>%
    mutate(label = gsub("_", " ", sub("HALLMARK_", "", pathway)))
  gsea_plot <- bind_rows(slice_min(gsea, NES, n = 6), slice_max(gsea, NES, n = 6)) %>%
    distinct(pathway, .keep_all = TRUE) %>% arrange(NES) %>%
    mutate(label = factor(label, levels = label))
  p_a <- ggplot(gsea_plot, aes(NES, label, fill = NES > 0)) +
    geom_col(width = 0.70) +
    geom_vline(xintercept = 0, linewidth = 0.32, colour = "#333333") +
    geom_text(
      aes(label = paste0("FDR ", format_fdr(padj))),
      hjust = ifelse(gsea_plot$NES > 0, 1.04, -0.04), size = 1.55,
      colour = "white", fontface = "bold", family = FONT
    ) +
    scale_fill_manual(values = c(`TRUE` = "#DF7621", `FALSE` = "#2B8CBE"), guide = "none") +
    labs(
      title = "Competitive Hallmark enrichment",
      subtitle = "Gene-set-membership null; participants are not resampled",
      x = "Normalized enrichment score", y = NULL
    ) +
    theme_pub()

  lopo <- fread(file.path(ROBUST, "fgsea_selected_pathways_by_omission.csv")) %>%
    mutate(
      label = recode(
        pathway,
        HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection",
        HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon-gamma response",
        HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon-alpha response",
        HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
      )
    )
  lopo_ranges <- lopo %>% filter(fit_id != "full") %>%
    group_by(pathway, label) %>%
    summarise(
      lo = min(NES), hi = max(NES), median = median(NES),
      n_fdr = sum(padj < 0.05), .groups = "drop"
    )
  lopo_full <- gsea %>%
    filter(pathway %in% unique(lopo$pathway)) %>%
    select(pathway, full_nes = NES)
  lopo_ranges <- left_join(lopo_ranges, lopo_full, by = "pathway") %>%
    mutate(label = factor(label, levels = rev(c(
      "Allograft rejection", "Interferon-gamma response",
      "Interferon-alpha response", "Inflammatory response"
    ))))
  p_b <- ggplot(lopo_ranges, aes(lo, label)) +
    geom_segment(aes(xend = hi, yend = label), linewidth = 1.1, colour = "#BFD3E6", lineend = "round") +
    geom_point(aes(x = median), shape = 21, size = 2.1, fill = "white", colour = "#4C78A8") +
    geom_point(aes(x = full_nes), shape = 18, size = 2.2, colour = "#D95F59") +
    geom_text(aes(x = median, label = paste0(n_fdr, "/6")), vjust = -1.05, size = 1.55, colour = "#3F4C59") +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(
      title = "Leave-one-participant-out sensitivity",
      subtitle = paste0(
        "Six omissions; <span style='color:#D95F59;'>&#9670;</span> ",
        "full fit<br>Text = FDR&lt;0.05 count"
      ),
      x = "Normalized enrichment score", y = NULL
    ) +
    theme_pub() +
    theme(
      plot.subtitle = ggtext::element_markdown(
        size = 5.7, colour = "#4B5563", lineheight = 1.05,
        margin = margin(b = 3)
      )
    )

  p_c <- sampling_panel(
    read_modules(), c("EC_inflammation", "cGAS_STING", "YAP_mTOR"),
    ncol = 3, y_label = "Mean module score"
  ) +
    labs(title = "Inflammatory and YAP-axis scores by sampling unit") +
    theme(
      # Use the full-width row rather than inheriting the large left-side
      # alignment allowance from panel A.  The larger facet gap keeps the
      # three sampling-unit summaries visually distinct and evenly spaced.
      panel.spacing.x = unit(7.5, "mm"),
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, margin = margin(b = 2)),
      plot.margin = margin(0, 5, 2, 2)
    )

  focus <- c("RELA", "NFKB1", "ATF4", "FOXO3", "STAT3", "HIF1A")
  tf <- fread(file.path(INPUT, "tf_stats_ulm_v4.csv")) %>%
    filter(TF %in% focus) %>%
    mutate(TF = factor(TF, levels = rev(focus)))
  p_d <- ggplot(tf, aes(ONFH_3A_vs_HOA_hodges_lehmann, TF)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#9CA3AF") +
    geom_segment(aes(x = 0, xend = ONFH_3A_vs_HOA_hodges_lehmann, yend = TF), linewidth = 0.55, colour = "#9EC3D5") +
    geom_point(size = 2.0, colour = "#2B8CBE") +
    labs(
      title = "Signed TF-activity summary",
      subtitle = "ARCO 3A minus HOA; no tested regulon had FDR<0.05",
      x = "Hodges-Lehmann shift", y = NULL
    ) +
    theme_pub() +
    theme(
      plot.title.position = "plot",
      plot.margin = margin(0, 4, 4, 1)
    )

  top <- fread(file.path(INPUT, "comm_top_ONFH_4_vs_ONFH_3A_v4.csv"))[1:12]
  keys <- unique(top[, .(sender, receiver, pair)])
  comm <- fread(file.path(INPUT, "comm_scores_v4_long.csv.gz")) %>%
    filter(independent_for_inference, group %in% c("ONFH_3A", "ONFH_4")) %>%
    inner_join(keys, by = c("sender", "receiver", "pair")) %>%
    mutate(
      label = paste0(sender, " to ", receiver, ": ", pair),
      label = factor(label, levels = rev(paste0(top$sender, " to ", top$receiver, ": ", top$pair))),
      group = factor(group, levels = c("ONFH_3A", "ONFH_4"))
    )
  p_e <- ggplot(comm, aes(score, label, colour = group)) +
    geom_point(
      position = position_jitter(height = 0.09, width = 0, seed = 20260903),
      size = 1.2, alpha = 0.90
    ) +
    scale_colour_manual(
      values = GROUP_COLORS[c("ONFH_3A", "ONFH_4")],
      labels = c("ARCO 3A", "ARCO 4"), name = NULL
    ) +
    labs(
      title = "EC-centered communication scores",
      subtitle = "Independent Liao participants; no tested interaction had FDR<0.05",
      x = "Ligand-receptor score", y = NULL
    ) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.margin = margin(0, 5, 4, 10)
    )

  # `free()` prevents the long labels in the upper panels from forcing panels
  # C and D into the same oversized left gutter.  This pulls their plotting
  # areas toward their panel tags without changing any quantitative content.
  middle_row <- free(p_c, side = "lr")
  bottom_row <- (free(p_d, side = "l") | p_e) +
    plot_layout(widths = c(0.68, 1.32))

  figure <- ((p_a | p_b) + plot_layout(widths = c(1.18, 0.82))) /
    middle_row /
    bottom_row +
    plot_layout(heights = c(1.05, 0.76, 1.03)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_pub(figure, "Figure4", 183, 176)
}

make_figure7 <- function() {
  permutation <- fread(file.path(INPUT, "diag_permutation_v4.csv"))
  stability <- fread(file.path(INPUT, "diag_feature_stability_v4.csv"))
  names(stability)[1:2] <- c("gene", "selection_frequency")
  stability <- stability %>% arrange(desc(selection_frequency)) %>% slice_head(n = 15) %>%
    mutate(gene = factor(gene, levels = rev(gene)))
  nested <- fread(file.path(INPUT, "diag_nested_cv_performance_v7.csv")) %>%
    mutate(model = "Mitochondrial candidate space")
  ma <- fread(file.path(INPUT, "diag_ma_comparator_repeat_performance_v7.csv")) %>%
    mutate(model = "Ma four-gene comparator")
  perf <- bind_rows(nested, ma) %>%
    pivot_longer(c(AUC, average_precision), names_to = "metric", values_to = "value") %>%
    mutate(
      metric = recode(metric, AUC = "AUC", average_precision = "Average precision"),
      model = factor(
        model,
        levels = c("Mitochondrial candidate space", "Ma four-gene comparator"),
        labels = c("Candidate space", "Ma comparator")
      )
    )
  summary <- jsonlite::read_json(file.path(INPUT, "diag_summary_v4.json"), simplifyVector = TRUE)
  delong <- fread(file.path(INPUT, "diag_delong_model_comparison_v8.csv"))[1]
  paired <- fread(file.path(INPUT, "diag_paired_model_comparison_v8.csv"))[1]
  predictions <- fread(file.path(INPUT, "diag_oof_predictions_aggregated_v8.csv"))

  p_a <- ggplot(permutation, aes(perm_auc)) +
    geom_histogram(bins = 26, fill = "#9ECAE1", colour = "white", linewidth = 0.25) +
    geom_vline(xintercept = summary$aggregate_oof_auc, colour = "#CB181D", linewidth = 0.85) +
    annotate(
      "text", x = summary$aggregate_oof_auc, y = Inf,
      label = sprintf("Observed AUC %.3f\nempirical p=%.4f", summary$aggregate_oof_auc, summary$permutation_empirical_p),
      hjust = 1.05, vjust = 1.2, size = 1.75, family = FONT, colour = "#8B1A1A"
    ) +
    labs(
      title = "Matched full-pipeline label-permutation null",
      x = "Aggregate repeated-OOF AUC", y = "Permutation count"
    ) +
    theme_pub()

  p_b <- ggplot(stability, aes(selection_frequency, gene)) +
    geom_col(width = 0.70, fill = "#2B8CBE") +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(title = "Feature-selection stability", x = "Selection across 25 outer fits", y = NULL) +
    theme_pub()

  p_c <- ggplot(perf, aes(model, value, group = .data[["repeat"]], colour = model)) +
    geom_line(colour = "#BFC5CC", linewidth = 0.35) +
    geom_point(size = 1.55) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = c(
      "Candidate space" = "#0072B2",
      "Ma comparator" = "#CC79A7"
    ), guide = "none") +
    labs(
      title = "Paired repeat-level performance",
      subtitle = "Identical outer folds were used for both models",
      x = NULL, y = "Metric value"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  roc_candidate <- roc(
    predictions$label, predictions$candidate_mean_oof_probability,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  roc_ma <- roc(
    predictions$label, predictions$ma_mean_oof_probability,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  roc_frame <- bind_rows(
    data.frame(
      fpr = 1 - roc_candidate$specificities, tpr = roc_candidate$sensitivities,
      model = "Mitochondrial candidate space"
    ),
    data.frame(
      fpr = 1 - roc_ma$specificities, tpr = roc_ma$sensitivities,
      model = "Ma four-gene comparator"
    )
  )
  roc_frame$model <- factor(
    roc_frame$model,
    levels = c("Mitochondrial candidate space", "Ma four-gene comparator"),
    labels = c("Candidate space", "Ma comparator")
  )
  p_d <- ggplot(roc_frame, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9CA3AF", linewidth = 0.35) +
    geom_step(linewidth = 0.75) +
    coord_equal() +
    scale_colour_manual(values = c(
      "Candidate space" = "#0072B2",
      "Ma comparator" = "#CC79A7"
    ), labels = c(
      "Candidate space" = "Candidate space · AUC 0.870",
      "Ma comparator" = "Ma comparator · AUC 0.813"
    ), name = NULL) +
    guides(colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(linewidth = 1.0)
    )) +
    annotate(
      "label", x = 0.98, y = 0.055,
      label = sprintf(
        "Delta AUC %.3f\n95%% CI %.3f to %.3f\npaired DeLong p=%.3f",
        delong$delta_auc, paired$ci_low, paired$ci_high, delong$delong_p_two_sided
      ),
      hjust = 1, vjust = 0, size = 1.65, family = FONT,
      fill = "white", colour = "#333333", linewidth = 0.18
    ) +
    labs(
      title = "Aggregated cross-fitted ROC",
      subtitle = "Paired fixed-prediction comparison",
      x = "False-positive rate", y = "True-positive rate"
    ) +
    theme_pub() +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.text = element_text(size = 5.0),
      legend.key.width = unit(4.0, "mm"),
      legend.key.height = unit(2.0, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      legend.spacing.y = unit(0, "mm"),
      legend.margin = margin(0, 0, 1, 0)
    )

  figure <- ((p_a | p_b) / (p_c | p_d)) +
    plot_layout(widths = c(1, 1), heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  save_pub(figure, "Figure7", 183, 140)
}

figure_only <- Sys.getenv("FIGURE_ONLY", unset = "")
if (nzchar(figure_only)) {
  figure_builders <- list(
    `2` = make_figure2,
    `3` = make_figure3,
    `4` = make_figure4,
    `7` = make_figure7
  )
  if (!figure_only %in% names(figure_builders)) {
    stop("FIGURE_ONLY must be one of: 2, 3, 4, 7")
  }
  figure_builders[[figure_only]]()
} else {
  make_figure2()
  make_figure3()
  make_figure4()
  make_figure7()
}

writeLines(
  c(
    "Genes figure QA",
    "Backend: R only for Figures 2, 3, 4 and 7.",
    "Figure 2: duplicate module scores removed.",
    "Figure 3: mitochondrial modules retained here only.",
    "Figure 4: the full fit uses a red diamond; redundant FDR=1.00 labels are hidden.",
    "Figure 7: paired candidate-versus-Ma comparison includes fixed-prediction uncertainty.",
    "All plotted statistics are read from versioned CSV/JSON inputs."
  ),
  file.path(QA, "Genes_revision_figures_QA_notes.txt")
)

cat(normalizePath(OUT, winslash = "/"), "\n")
