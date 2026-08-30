#!/usr/bin/env Rscript

# Submission-width redraw of Figure 5 using the official-R scTenifoldKnk
# outputs. Core conclusion: donor-separated virtual SQSTM1 perturbation yields
# a parameter-audited candidate ranking, not experimental validation.
# Evidence chain: A, design; B, removed WT neighborhood (hero); C, gene-manifold
# displacement; D, matched-null pathway recovery. Archetype: asymmetric mixed-
# modality figure. R is the exclusive drawing/export backend. No analytical
# values, selected nodes, edge signs, ranks, p values, or FDR values are changed.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(ragg)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
ROOT <- dirname(dirname(script_file))
RESULTS <- file.path(ROOT, "results")
FIGURE_DATA <- file.path(RESULTS, "official_r_vko_figure_data")
PRIMARY <- file.path(RESULTS, "official_r_vko_manuscript")
OUT <- file.path(ROOT, "figures", "final")
QA_OUT <- file.path(ROOT, "qa")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(QA_OUT, recursive = TRUE, showWarnings = FALSE)

TARGET <- "SQSTM1"
# Generic sans maps to Arial/Helvetica-compatible metrics across the Windows
# ragg, Cairo PDF and base SVG devices without PostScript font warnings.
FONT <- "sans"

NODE_COLORS <- c(
  "EC inflammation" = "#E45756",
  "OXPHOS" = "#4C78A8",
  "Angiogenesis" = "#59A14F",
  "Mitophagy" = "#F2CF5B",
  "Mito stress" = "#B279A2",
  "Other" = "#D9D9D9"
)
EDGE_COLORS <- c(
  "positive" = "#D95F59",
  "negative" = "#4C78A8",
  "discordant" = "#9CA3AF"
)

theme_pub <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      plot.title = element_text(size = 7.2, face = "plain", margin = margin(b = 3, l = 14)),
      plot.tag = element_text(size = 8.5, face = "bold"),
      plot.tag.position = c(0.004, 0.995),
      plot.caption = element_text(size = 5.2, colour = "#4B5563", hjust = 0, margin = margin(t = 2)),
      axis.title = element_text(size = 6.2),
      axis.text = element_text(size = 5.5),
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      strip.text = element_text(size = 6.2, face = "bold"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 5.1),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(4.2, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

pathway_table <- read.csv(
  file.path(PRIMARY, "vko_sqstm1_pathway_enrichment_official_r.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pathway_sets <- setNames(
  lapply(pathway_table$genes, function(x) strsplit(x, ";", fixed = TRUE)[[1]]),
  pathway_table$pathway
)

gene_category <- function(gene) {
  if (gene %in% pathway_sets[["EC_inflammation"]]) return("EC inflammation")
  if (gene %in% pathway_sets[["OXPHOS"]]) return("OXPHOS")
  if (gene %in% pathway_sets[["Angiogenesis"]]) return("Angiogenesis")
  if (gene %in% pathway_sets[["Mitophagy_core"]]) return("Mitophagy")
  stress <- unique(c(
    pathway_sets[["ROS_defense"]], pathway_sets[["mtDNA_release"]],
    pathway_sets[["cGAS_STING"]], pathway_sets[["Mito_fission"]],
    pathway_sets[["Mito_proteostasis"]]
  ))
  if (gene %in% stress) return("Mito stress")
  "Other"
}

read_primary <- function(donor) {
  read.csv(
    file.path(PRIMARY, sprintf("vko_sqstm1_%s_official_r.csv", donor)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

edge_table <- function() {
  h2 <- read.csv(file.path(FIGURE_DATA, "sqstm1_outgoing_edges_hoa2.csv"), stringsAsFactors = FALSE) %>%
    transmute(gene, weight_hoa2 = wt_outgoing_weight, abs_hoa2 = absolute_wt_weight)
  h3 <- read.csv(file.path(FIGURE_DATA, "sqstm1_outgoing_edges_hoa3.csv"), stringsAsFactors = FALSE) %>%
    transmute(gene, weight_hoa3 = wt_outgoing_weight, abs_hoa3 = absolute_wt_weight)
  inner_join(h2, h3, by = "gene") %>%
    mutate(
      scaled_hoa2 = abs_hoa2 / max(abs_hoa2, na.rm = TRUE),
      scaled_hoa3 = abs_hoa3 / max(abs_hoa3, na.rm = TRUE),
      consensus_edge_score = (scaled_hoa2 + scaled_hoa3) / 2,
      edge_sign = case_when(
        weight_hoa2 > 0 & weight_hoa3 > 0 ~ "positive",
        weight_hoa2 < 0 & weight_hoa3 < 0 ~ "negative",
        TRUE ~ "discordant"
      )
    ) %>%
    arrange(desc(consensus_edge_score))
}

spread_labels <- function(frame, side, lower = -1.28, upper = 1.28, gap = 0.17) {
  if (nrow(frame) == 0) return(frame)
  frame <- frame %>% arrange(y)
  y_new <- pmax(lower, pmin(upper, frame$y))
  if (length(y_new) > 1) {
    for (i in 2:length(y_new)) y_new[[i]] <- max(y_new[[i]], y_new[[i - 1]] + gap)
    if (max(y_new) > upper) y_new <- y_new - (max(y_new) - upper)
    for (i in (length(y_new) - 1):1) {
      if (i >= 1) y_new[[i]] <- min(y_new[[i]], y_new[[i + 1]] - gap)
    }
    if (min(y_new) < lower) y_new <- y_new + (lower - min(y_new))
  }
  frame$label_y <- y_new
  frame$label_x <- 1.55 * side
  frame$hjust <- ifelse(side > 0, 0, 1)
  frame
}

network_components <- function() {
  edges <- edge_table()
  stable <- read.csv(
    file.path(RESULTS, "official_r_vko_cross_profile_gene_summary.csv"),
    stringsAsFactors = FALSE
  )$Gene[1:12]
  significant <- unique(c(
    read_primary("hoa2") %>% filter(Gene != TARGET, adjusted_p_value < 0.05) %>% pull(Gene),
    read_primary("hoa3") %>% filter(Gene != TARGET, adjusted_p_value < 0.05) %>% pull(Gene)
  ))

  selected <- character(0)
  add_nodes <- function(x) {
    selected <<- unique(c(selected, setdiff(x, TARGET)))
  }
  add_nodes(stable)
  add_nodes(head(edges$gene, 30))
  categorized <- edges %>% mutate(category = vapply(gene, gene_category, character(1)))
  for (cat in names(NODE_COLORS)) add_nodes(head(categorized$gene[categorized$category == cat], 4))
  selected <- head(selected, 42)

  nodes <- edges %>%
    filter(gene %in% selected) %>%
    mutate(
      category = vapply(gene, gene_category, character(1)),
      category_order = match(category, names(NODE_COLORS)),
      stable = gene %in% stable,
      significant = gene %in% significant
    ) %>%
    arrange(category_order, desc(consensus_edge_score)) %>%
    mutate(
      index = row_number() - 1,
      angle = pi / 2 + 2 * pi * index / n(),
      radius = ifelse(index %% 2 == 0, 1.02, 1.32),
      x = radius * cos(angle),
      y = radius * sin(angle),
      ux = x / sqrt(x^2 + y^2),
      uy = y / sqrt(x^2 + y^2),
      x_start = 0.36 * ux,
      y_start = 0.36 * uy,
      x_end = x - 0.070 * ux,
      y_end = y - 0.070 * uy,
      node_size = 2.2 + 2.6 * consensus_edge_score
    )

  label_genes <- unique(c(
    "NFKBIA", "ICAM1", "TNFAIP3", "SELE", "MT-ND1", "MT-ATP6",
    "MT-ND4", "MT-CO2", "MT-CO1", "VWF", "ENG", "SERPINE1",
    "CXCL2", "MCL1", "KDM6B", "TIPARP", significant
  ))
  top_genes <- c("ICAM1", "TNFAIP3", "SELE", "NFKBIA")
  label_nodes <- nodes %>%
    filter(gene %in% label_genes) %>%
    mutate(side = ifelse(x >= 0, 1, -1))
  top_labels <- label_nodes %>%
    filter(gene %in% top_genes) %>%
    mutate(
      label_x = c("ICAM1" = -0.82, "TNFAIP3" = -0.30, "SELE" = 0.18, "NFKBIA" = 0.75)[gene],
      label_y = c("ICAM1" = 1.46, "TNFAIP3" = 1.67, "SELE" = 1.46, "NFKBIA" = 1.67)[gene],
      side = 0,
      hjust = 0.5
    )
  side_labels <- filter(label_nodes, !gene %in% top_genes)
  labels <- bind_rows(
    spread_labels(filter(side_labels, side < 0), -1),
    spread_labels(filter(side_labels, side > 0), 1),
    top_labels
  ) %>%
    mutate(fontface = ifelse(stable, "bold.italic", "italic"))

  list(nodes = nodes, labels = labels)
}

network_plot <- function() {
  comp <- network_components()
  nodes <- comp$nodes
  labels <- comp$labels

  ggplot() +
    geom_segment(
      data = nodes,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
          colour = edge_sign, linewidth = consensus_edge_score),
      alpha = 0.66, lineend = "round",
      arrow = arrow(type = "closed", length = unit(1.15, "mm")),
      show.legend = c(colour = TRUE, linewidth = FALSE)
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, size = node_size, fill = category),
      shape = 21, colour = "white", stroke = 0.35,
      show.legend = c(fill = TRUE, size = FALSE)
    ) +
    geom_point(
      data = filter(nodes, significant),
      aes(x = x, y = y, size = node_size),
      shape = 21, fill = NA, colour = "#1F2937", stroke = 0.8,
      show.legend = FALSE
    ) +
    geom_segment(
      data = labels,
      aes(x = x * 1.04, y = y * 1.04,
          xend = label_x - 0.035 * side, yend = label_y),
      colour = "#7A8491", linewidth = 0.22, alpha = 0.85,
      show.legend = FALSE
    ) +
    geom_text(
      data = labels,
      aes(x = label_x, y = label_y, label = gene, hjust = hjust, fontface = fontface),
      family = FONT, size = 2.05, colour = "#111827", show.legend = FALSE
    ) +
    annotate(
      "point", x = 0, y = 0, shape = 21, size = 16.4,
      fill = "#F6D6D3", colour = "#B23A33", stroke = 1.0
    ) +
    annotate(
      "text", x = 0, y = 0, label = "SQSTM1",
      family = FONT, fontface = "bold", size = 2.62, colour = "#111827"
    ) +
    scale_fill_manual(name = "Node class", values = NODE_COLORS, breaks = names(NODE_COLORS)) +
    scale_colour_manual(
      name = "Edge sign",
      values = EDGE_COLORS,
      breaks = names(EDGE_COLORS),
      labels = c("Positive edge", "Negative edge", "Discordant edge")
    ) +
    scale_size_identity() +
    scale_linewidth_continuous(range = c(0.22, 0.78), guide = "none") +
    guides(
      fill = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = list(size = 3.0)),
      colour = guide_legend(order = 2, nrow = 1, override.aes = list(linewidth = 0.8))
    ) +
    coord_equal(xlim = c(-1.86, 1.86), ylim = c(-1.62, 1.82), clip = "off") +
    labs(title = "SQSTM1 WT neighborhood removed by virtual KO", tag = "B") +
    theme_pub() +
    theme(
      axis.line = element_blank(), axis.text = element_blank(),
      axis.title = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 5.2, face = "bold"),
      legend.spacing.y = unit(-1.2, "mm"),
      plot.margin = margin(4, 2, 1, 2)
    )
}

# Minimal revision requested after review: retain the original radial visual
# grammar and change only the two actual collision sources (central edge field
# and a handful of labels). This candidate is previewed separately and is not
# inserted into the manuscript until approved.
network_plot_minimal <- function() {
  comp <- network_components()
  nodes <- comp$nodes

  label_genes <- c(
    "NFKBIA", "ICAM1", "TNFAIP3", "SELE", "MT-ND1", "MT-ATP6",
    "MT-ND4", "MT-CO2", "MT-CO1", "VWF", "ENG", "SERPINE1",
    "CXCL2", "MCL1", "KDM6B", "TIPARP"
  )
  label_positions <- tibble::tribble(
    ~gene,      ~label_x, ~label_y, ~hjust,
    "ICAM1",      -0.18,     1.23,     0.5,
    "TNFAIP3",    -0.61,     1.10,     1.0,
    "SELE",       -0.29,     0.93,     0.5,
    "NFKBIA",      0.16,     1.08,     0.0,
    "MT-ATP6",    -1.02,     0.78,     1.0,
    "MT-ND1",     -1.02,     0.62,     1.0,
    "MT-ND4",     -1.04,     0.45,     1.0,
    "MT-CO2",     -1.14,     0.29,     1.0,
    "MT-CO1",     -1.08,     0.12,     1.0,
    "VWF",        -1.16,    -0.05,     1.0,
    "ENG",        -1.00,    -0.23,     1.0,
    "SERPINE1",    1.11,     0.00,     0.0,
    "CXCL2",       0.84,    -0.43,     0.0,
    "TIPARP",      0.69,    -0.66,     0.0,
    "KDM6B",       0.65,    -0.94,     0.0,
    "MCL1",       -0.28,    -1.16,     1.0
  )
  labels <- nodes %>%
    filter(gene %in% label_genes) %>%
    left_join(label_positions, by = "gene") %>%
    mutate(fontface = ifelse(stable, "bold.italic", "italic"))

  ggplot() +
    geom_segment(
      data = nodes,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
          colour = edge_sign, linewidth = consensus_edge_score),
      alpha = 0.62, lineend = "round",
      arrow = arrow(type = "closed", length = unit(1.05, "mm")),
      show.legend = c(colour = TRUE, linewidth = FALSE)
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, size = node_size, fill = category),
      shape = 21, colour = "white", stroke = 0.3,
      show.legend = c(fill = TRUE, size = FALSE)
    ) +
    geom_point(
      data = filter(nodes, significant),
      aes(x = x, y = y, size = node_size),
      shape = 21, fill = NA, colour = "#1F2937", stroke = 0.75,
      show.legend = FALSE
    ) +
    annotate(
      "point", x = 0, y = 0, shape = 21, size = 12.2,
      fill = "#F6D6D3", colour = "#B23A33", stroke = 0.9
    ) +
    annotate(
      "text", x = 0, y = 0, label = "SQSTM1",
      family = FONT, fontface = "bold", size = 2.75, colour = "#111827"
    ) +
    geom_segment(
      data = labels,
      aes(x = x, y = y, xend = label_x, yend = label_y),
      colour = "#A3AAB3", linewidth = 0.18, alpha = 0.8,
      lineend = "round", show.legend = FALSE
    ) +
    geom_text(
      data = labels,
      aes(x = label_x, y = label_y, label = gene, hjust = hjust,
          fontface = fontface),
      family = FONT, size = 1.95, colour = "#111827",
      vjust = 0.5, show.legend = FALSE
    ) +
    scale_fill_manual(name = "Node class", values = NODE_COLORS, breaks = names(NODE_COLORS)) +
    scale_colour_manual(
      name = "Edge sign", values = EDGE_COLORS, breaks = names(EDGE_COLORS),
      labels = c("Positive edge", "Negative edge", "Discordant edge")
    ) +
    scale_size_identity() +
    scale_linewidth_continuous(range = c(0.20, 0.72), guide = "none") +
    guides(
      fill = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = list(size = 2.8)),
      colour = guide_legend(order = 2, nrow = 1, override.aes = list(linewidth = 0.75))
    ) +
    coord_equal(xlim = c(-1.46, 1.46), ylim = c(-1.32, 1.32), clip = "off") +
    labs(title = "SQSTM1 WT neighborhood removed by virtual KO", tag = "B") +
    theme_pub() +
    theme(
      axis.line = element_blank(), axis.text = element_blank(),
      axis.title = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 5.0, face = "bold"),
      legend.text = element_text(size = 4.9),
      legend.spacing.y = unit(-1.0, "mm"),
      plot.margin = margin(4, 3, 1, 3)
    )
}

workflow_plot <- function() {
  boxes <- tibble::tribble(
    ~xmin, ~xmax, ~ymin, ~ymax, ~label, ~fill, ~text_size,
    0.4, 3.4, 5.65, 7.55, "HOA2 control ECs\nn = 1,043\nSQSTM1+ = 70.8%", "#DDEBF7", 1.52,
    0.4, 3.4, 3.25, 5.15, "HOA3 control ECs\nn = 760\nSQSTM1+ = 72.2%", "#E2F0D9", 1.52,
    4.6, 7.7, 4.50, 7.50, "Official R packages\nscTenifoldKnk v1.1\nscTenifoldNet v1.4\n300 shared genes", "#FFF2CC", 1.42,
    8.8, 11.8, 4.15, 7.85, "Virtual SQSTM1 KO\nSQSTM1 network\nrow = 0\nWT-KO manifold\nalignment\nDifferential regulation", "#FCE4D6", 1.32
  )
  arrows <- tibble::tribble(
    ~x, ~y, ~xend, ~yend,
    3.4, 6.60, 4.6, 6.20,
    3.4, 4.20, 4.6, 5.80,
    7.7, 6.00, 8.8, 6.00
  )

  ggplot() +
    geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      colour = "#4B5563", linewidth = 0.35, show.legend = FALSE
    ) +
    geom_text(
      data = boxes,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label, size = text_size),
      family = FONT, lineheight = 1.00
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#4B5563", linewidth = 0.55,
      arrow = arrow(type = "closed", length = unit(1.4, "mm"))
    ) +
    annotate(
      "rect", xmin = 0.4, xmax = 11.8, ymin = 8.35, ymax = 9.45,
      fill = "#F3F4F6", colour = "#9CA3AF", linewidth = 0.25
    ) +
    annotate(
      "text", x = 0.65, y = 9.08, label = "Eligibility screening",
      family = FONT, size = 1.75, fontface = "bold", colour = "#374151", hjust = 0
    ) +
    annotate(
      "text", x = 0.65, y = 8.66,
      label = "HOA1 excluded from network modeling (26 ECs after filtering)",
      family = FONT, size = 1.38, colour = "#4B5563", hjust = 0
    ) +
    annotate("text", x = 0.4, y = 1.92, label = "Parameter audit", family = FONT,
             size = 2.45, fontface = "bold", hjust = 0) +
    annotate(
      "text", x = 0.4, y = 0.60,
      label = paste0(
        "Primary: 20 networks, q=0.95, 30 manifold dimensions\n",
        "Package-default sensitivity: 10 networks, q=0.90, 2 dimensions\n",
        "Calibration: 295-gene mtDNA-feature-excluded refits + 20 matched-comparator KOs (Table S10)\n",
        "No downstream gene or pathway was FDR-significant in both donors under both profiles."
      ),
      family = FONT, size = 1.34, lineheight = 1.04, hjust = 0
    ) +
    coord_cartesian(xlim = c(0, 12.2), ylim = c(0, 10), clip = "off") +
    scale_fill_identity() +
    scale_size_identity() +
    labs(title = "Official-package virtual-knockout design", tag = "A") +
    theme_pub() +
    theme(
      axis.line = element_blank(), axis.text = element_blank(),
      axis.title = element_blank(), axis.ticks = element_blank(),
      plot.margin = margin(4, 4, 5, 2)
    )
}

manifold_for_donor <- function(donor) {
  manifold <- read.csv(
    file.path(FIGURE_DATA, sprintf("wt_ko_manifold_%s.csv", donor)),
    stringsAsFactors = FALSE
  ) %>%
    pivot_wider(names_from = state, values_from = c(PC1, PC2), names_sep = "_")
  results <- read_primary(donor)
  frame <- inner_join(manifold, results, by = c("gene" = "Gene")) %>%
    filter(gene != TARGET) %>%
    mutate(dx = PC1_KO - PC1_WT, dy = PC2_KO - PC2_WT)

  xq <- quantile(frame$PC1_WT, c(0.01, 0.99), na.rm = TRUE)
  yq <- quantile(frame$PC2_WT, c(0.01, 0.99), na.rm = TRUE)
  xpad <- max(diff(xq) * 0.10, 1e-4)
  ypad <- max(diff(yq) * 0.10, 1e-4)
  robust_span <- max(diff(xq), diff(yq))
  reference <- max(quantile(sqrt(frame$dx^2 + frame$dy^2), 0.95, na.rm = TRUE), 1e-12)
  factor <- 0.045 * robust_span / reference

  frame %>%
    mutate(
      donor = toupper(donor), display_factor = factor,
      xend = PC1_WT + dx * factor,
      yend = PC2_WT + dy * factor,
      visible = between(PC1_WT, xq[[1]] - xpad, xq[[2]] + xpad) &
                between(PC2_WT, yq[[1]] - ypad, yq[[2]] + ypad),
      top10 = rank <= 10,
      significant = adjusted_p_value < 0.05,
      xmin = xq[[1]] - xpad, xmax = xq[[2]] + xpad,
      ymin = yq[[1]] - ypad, ymax = yq[[2]] + ypad
    )
}

manifold_plot <- function() {
  frame <- bind_rows(manifold_for_donor("hoa2"), manifold_for_donor("hoa3"))
  visible <- filter(frame, visible)
  top <- filter(visible, top10)
  sig <- filter(visible, significant)
  sig_label_positions <- tibble::tribble(
    ~donor, ~gene,     ~label_x, ~label_y, ~hjust,
    "HOA2", "MT-ND1",  -0.012,   0.0140,     0.5,
    "HOA2", "MT-ND4",  -0.015,  -0.0100,     1.0,
    "HOA2", "MT-CO2",   0.006,  -0.0010,     0.0,
    "HOA2", "MT-ATP6", -0.008,  -0.0250,     0.5,
    "HOA3", "MT-ND1",  -0.040,   0.0180,     0.5,
    "HOA3", "VWF",      0.035,  -0.0160,     0.5
  )
  sig_labels <- sig %>%
    left_join(sig_label_positions, by = c("donor", "gene"))
  bounds <- frame %>%
    group_by(donor) %>%
    summarise(xmin = first(xmin), xmax = first(xmax), ymin = first(ymin), ymax = first(ymax), .groups = "drop")
  factors <- frame %>%
    group_by(donor) %>%
    summarise(x = first(xmin), y = first(ymin), factor = first(display_factor), .groups = "drop") %>%
    mutate(label = sprintf("Arrows enlarged %sx", format(round(factor), big.mark = ",")))

  ggplot(visible, aes(PC1_WT, PC2_WT)) +
    geom_blank(data = bounds, aes(x = xmin, y = ymin)) +
    geom_blank(data = bounds, aes(x = xmax, y = ymax)) +
    geom_point(size = 0.42, colour = "#D1D5DB", alpha = 0.62) +
    geom_segment(
      aes(xend = xend, yend = yend), colour = "#9CA3AF", alpha = 0.34,
      linewidth = 0.13, arrow = arrow(type = "closed", length = unit(0.65, "mm"))
    ) +
    geom_segment(
      data = top,
      aes(xend = xend, yend = yend, colour = donor),
      alpha = 0.92, linewidth = 0.48,
      arrow = arrow(type = "closed", length = unit(0.95, "mm"))
    ) +
    geom_point(
      data = sig, aes(x = xend, y = yend), shape = 21, fill = NA,
      colour = "#D62728", stroke = 0.55, size = 2.0
    ) +
    geom_segment(
      data = sig_labels,
      aes(x = xend, y = yend, xend = label_x, yend = label_y),
      inherit.aes = FALSE, colour = "#6B7280", linewidth = 0.18,
      show.legend = FALSE
    ) +
    geom_text(
      data = sig_labels,
      aes(x = label_x, y = label_y, label = gene, hjust = hjust),
      inherit.aes = FALSE, family = FONT, fontface = "italic", size = 1.68,
      vjust = 0.5, colour = "#111827", show.legend = FALSE
    ) +
    geom_text(
      data = factors, aes(x = x, y = y, label = label),
      family = FONT, size = 1.65, colour = "#4B5563", hjust = -0.02, vjust = -0.5
    ) +
    facet_wrap(~donor, nrow = 1, scales = "free") +
    scale_colour_manual(values = c("HOA2" = "#4C78A8", "HOA3" = "#59A14F"), guide = "none") +
    labs(
      title = "WT-to-virtual-KO gene-manifold displacement",
      tag = "C", x = "Gene-manifold PC1", y = "Gene-manifold PC2",
      caption = paste0(
        "Gray: downstream genes; color: top 10 ranks; red rings: FDR<0.05.\n",
        "Arrows are uniformly enlarged gene-manifold displacements."
      )
    ) +
    theme_pub() +
    theme(
      strip.background = element_blank(), strip.text = element_text(size = 6.0),
      plot.title = element_text(size = 6.7, margin = margin(b = 3, l = 14)),
      plot.caption = element_text(size = 4.65, lineheight = 0.95, hjust = 0),
      panel.spacing = unit(3.5, "mm"), plot.margin = margin(4, 5, 2, 2)
    )
}

recovery_frame <- function(pathway, donor, display, panel_colour) {
  results <- read_primary(donor) %>% arrange(rank)
  row <- pathway_table[pathway_table$pathway == pathway, , drop = FALSE]
  genes <- setdiff(strsplit(row$genes[[1]], ";", fixed = TRUE)[[1]], TARGET)
  member_ranks <- sort(results$rank[results$Gene %in% genes])
  n_total <- nrow(results)
  x <- c(0, member_ranks / n_total, 1)
  y <- c(0, seq_along(member_ranks) / length(member_ranks), 1)
  panel <- sprintf("%s\n%s primary", display, toupper(donor))
  p_value <- row[[sprintf("empirical_p_%s", donor)]][[1]]
  fdr_value <- row[[sprintf("empirical_fdr_%s", donor)]][[1]]
  tibble(
    panel = panel, pathway = pathway, donor = donor, colour_key = panel_colour,
    x = x, y = y,
    ymin = ifelse(y >= x, x, NA_real_), ymax = ifelse(y >= x, y, NA_real_),
    empirical_p = p_value,
    empirical_fdr = fdr_value
  )
}

recovery_plot <- function() {
  frame <- bind_rows(
    recovery_frame("EC_inflammation", "hoa3", "EC inflammation", "inflam"),
    recovery_frame("OXPHOS", "hoa2", "Oxidative phosphorylation", "oxphos"),
    recovery_frame("Mitophagy_core", "hoa3", "Mitophagy core", "mitophagy")
  )
  panel_order <- unique(frame$panel)
  frame$panel <- factor(frame$panel, levels = panel_order)
  rugs <- frame %>%
    group_by(panel, pathway, donor, colour_key) %>%
    summarise(x = list(x[x > 0 & x < 1]), .groups = "drop") %>%
    unnest(x)
  ann <- frame %>%
    group_by(panel, colour_key) %>%
    summarise(
      empirical_p = first(empirical_p), empirical_fdr = first(empirical_fdr),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("empirical p=%.3f\nBH FDR=%.3f", empirical_p, empirical_fdr))
  diagonal <- tibble(
    panel = factor(panel_order, levels = panel_order),
    x = 0, y = 0, xend = 1, yend = 1
  )
  ribbon <- filter(frame, !is.na(ymin), !is.na(ymax))

  ggplot(frame, aes(x, y, colour = colour_key, fill = colour_key, group = panel)) +
    geom_ribbon(data = ribbon, aes(ymin = ymin, ymax = ymax), alpha = 0.12, colour = NA, show.legend = FALSE) +
    geom_segment(
      data = diagonal, aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "#9CA3AF", linetype = "dashed", linewidth = 0.32
    ) +
    geom_step(linewidth = 0.62, direction = "hv", show.legend = FALSE) +
    geom_segment(
      data = rugs, aes(x = x, xend = x, y = -0.052, yend = 0),
      inherit.aes = FALSE, colour = "#111827", linewidth = 0.2
    ) +
    geom_text(
      data = ann, aes(x = 0.04, y = 0.95, label = label),
      inherit.aes = FALSE, family = FONT, size = 1.8, hjust = 0, vjust = 1
    ) +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = c(inflam = NODE_COLORS[["EC inflammation"]],
                                   oxphos = NODE_COLORS[["OXPHOS"]],
                                   mitophagy = NODE_COLORS[["Mitophagy"]])) +
    scale_fill_manual(values = c(inflam = NODE_COLORS[["EC inflammation"]],
                                 oxphos = NODE_COLORS[["OXPHOS"]],
                                 mitophagy = NODE_COLORS[["Mitophagy"]])) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.07, 1.02), clip = "off") +
    labs(
      title = "Matched-null pathway-rank recovery", tag = "D",
      x = "Perturbation rank percentile", y = "Cumulative fraction of pathway genes",
      caption = paste0(
        "Pathway-member recovery along official perturbation ranks.\n",
        "p/FDR use 20,000 expression- and prevalence-matched null sets."
      )
    ) +
    theme_pub() +
    theme(
      strip.background = element_blank(), strip.text = element_text(size = 5.9, face = "bold"),
      plot.title = element_text(size = 6.7, margin = margin(b = 3, l = 14)),
      plot.caption = element_text(size = 4.65, lineheight = 0.95, hjust = 0),
      panel.spacing = unit(3.0, "mm"), plot.margin = margin(4, 2, 2, 4)
    )
}

save_pub_r <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  if (identical(Sys.getenv("ONFH_EXPORT_ARCHIVAL", unset = "0"), "1")) {
    grDevices::svg(
      filename = paste0(stem, ".svg"), width = width_in, height = height_in,
      family = FONT, onefile = TRUE, bg = "white"
    )
    print(plot)
    dev.off()

    ragg::agg_tiff(
      filename = paste0(stem, ".tiff"), width = width_in, height = height_in,
      units = "in", res = dpi, compression = "lzw", background = "white"
    )
    print(plot)
    dev.off()
  }

  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = FONT, onefile = TRUE, bg = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    filename = paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  dev.off()
}

p_a <- workflow_plot()
p_b <- network_plot()
p_c <- manifold_plot()
p_d <- recovery_plot()

network_legend <- cowplot::get_legend(p_b + theme(legend.position = "bottom"))
p_b_core <- p_b + theme(legend.position = "none")
top_row <- (p_a | p_b_core) + plot_layout(widths = c(1.10, 1.55))
legend_row <- (plot_spacer() | wrap_elements(full = network_legend)) +
  plot_layout(widths = c(1.10, 1.55))
bottom_row <- (p_c | p_d) + plot_layout(widths = c(1.00, 1.55))

figure5 <- top_row / legend_row / bottom_row +
  plot_layout(heights = c(0.90, 0.15, 1.05)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_pub_r(figure5, file.path(OUT, "Figure5"), 183, 132, dpi = 600)

writeLines(
  c(
    "Figure 5 submission redraw",
    "Core conclusion: SQSTM1 virtual perturbation removes the WT outgoing neighborhood and yields parameter-audited candidate ranks, not experimental validation.",
    "Archetype: asymmetric mixed-modality figure.",
    "Backend: R only.",
    "Panel B: original radial grammar retained; enlarged central node contains the complete SQSTM1 label; inner/outer node radii are increased from 0.82/1.03 to 1.02/1.32; maximum node diameter is slightly reduced while preserving weight encoding; the legend is moved into a dedicated strip so the network receives more physical space.",
    "Panel A: all workflow boxes use explicit line breaks and expanded internal padding; arrows terminate at box boundaries; the four-line audit adds the mtDNA-feature-excluded and matched-comparator calibration note.",
    "Panel C: significant-gene labels use fixed donor-specific positions with short leader lines.",
    "Panel C/D captions are shortened to remain within their panel widths.",
    "Data and statistics: unchanged from official-R source tables."
  ),
  file.path(QA_OUT, "Figure5_QA_notes.txt")
)

cat(normalizePath(file.path(OUT, "Figure5.pdf"), winslash = "/"), "\n")
