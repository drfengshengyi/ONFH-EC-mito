#!/usr/bin/env Rscript

# Genes Figure 6: donor-separated, mtDNA-excluded SQSTM1 perturbation.
# Panels show design, gene ranks, donor agreement, matched controls, and pathways.
# CALCOCO2 is a secondary candidate; no experimental knockout is shown.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(patchwork)
  library(ragg)
  library(tidyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
ROOT <- dirname(dirname(script_file))
OUT <- file.path(ROOT, "figures", "final")
QA <- file.path(ROOT, "qa")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(QA, recursive = TRUE, showWarnings = FALSE)

FONT <- "sans"
COL_POS <- "#D95F59"
COL_NEG <- "#4C78A8"
COL_TEAL <- "#2A9D8F"
COL_GOLD <- "#E9A31B"
COL_GREY <- "#9AA3AD"

theme_pub <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 5.7, colour = "#4B5563", margin = margin(b = 3)),
      plot.tag = element_text(size = 8.8, face = "bold"),
      axis.title = element_text(size = 6.3),
      axis.text = element_text(size = 5.7, colour = "#222222"),
      axis.line = element_line(linewidth = 0.3, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.0, face = "bold"),
      legend.title = element_text(size = 5.7, face = "bold"),
      legend.text = element_text(size = 5.3),
      panel.grid.major = element_line(linewidth = 0.22, colour = "#E5E7EB"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 6, 4, 6)
    )
}

save_pub <- function(plot, name, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(
    file.path(OUT, paste0(name, ".svg")), width = width_in, height = height_in,
    bg = "white", system_fonts = list(sans = "Arial")
  )
  print(plot)
  dev.off()
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
  ragg::agg_tiff(
    file.path(OUT, paste0(name, ".tiff")), width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()
}

# A: separate HOA2 and HOA3 workflows.
lanes <- data.frame(
  donor = c("HOA2", "HOA3"),
  y = c(0.64, 0.34),
  n = c("1,043 ECs", "760 ECs"),
  fill = c("#EEF4FA", "#EDF7F4"),
  accent = c(COL_NEG, COL_TEAL)
)

lane_arrows <- expand.grid(
  y = lanes$y,
  start = c(0.255, 0.470, 0.735)
) %>%
  mutate(end = start + c(0.070, 0.070, 0.060)[match(start, c(0.255, 0.470, 0.735))])

network_nodes <- bind_rows(lapply(seq_len(nrow(lanes)), function(i) {
  data.frame(
    donor = lanes$donor[[i]],
    x = 0.39 + c(-0.030, 0.000, 0.031, -0.025, 0.026),
    y = lanes$y[[i]] + c(0.040, 0.064, 0.031, -0.040, -0.046),
    accent = lanes$accent[[i]]
  )
}))

network_edges <- bind_rows(lapply(seq_len(nrow(lanes)), function(i) {
  y0 <- lanes$y[[i]]
  data.frame(
    x = c(0.360, 0.390, 0.421, 0.365, 0.416, 0.390),
    y = y0 + c(0.040, 0.064, 0.031, -0.040, -0.046, 0.064),
    xend = c(0.390, 0.421, 0.416, 0.390, 0.390, 0.365),
    yend = y0 + c(0.064, 0.031, -0.046, 0.064, 0.064, -0.040)
  )
}))

rank_lines <- bind_rows(lapply(seq_len(nrow(lanes)), function(i) {
  offsets <- c(0.046, 0.024, 0.002, -0.020, -0.042)
  data.frame(
    x = 0.825,
    xend = c(0.955, 0.937, 0.918, 0.895, 0.875),
    y = lanes$y[[i]] + offsets,
    yend = lanes$y[[i]] + offsets,
    accent = lanes$accent[[i]]
  )
}))

p_a <- ggplot() +
  # Donor lanes.
  geom_rect(
    data = lanes,
    aes(xmin = 0.018, xmax = 0.982, ymin = y - 0.115, ymax = y + 0.115),
    fill = lanes$fill, colour = NA
  ) +
  # Input cards.
  geom_rect(
    data = lanes,
    aes(xmin = 0.035, xmax = 0.245, ymin = y - 0.078, ymax = y + 0.078),
    fill = "white", colour = "#C8D0D8", linewidth = 0.32
  ) +
  geom_rect(
    data = lanes,
    aes(xmin = 0.035, xmax = 0.047, ymin = y - 0.078, ymax = y + 0.078),
    fill = lanes$accent, colour = NA
  ) +
  geom_text(
    data = lanes, aes(x = 0.060, y = y + 0.026, label = donor),
    hjust = 0, size = 2.18, family = FONT, fontface = "bold",
    colour = lanes$accent
  ) +
  geom_text(
    data = lanes,
    aes(x = 0.060, y = y - 0.028, label = n),
    hjust = 0, size = 1.62, family = FONT, colour = "#3E4852"
  ) +
  # Stage flow.
  geom_segment(
    data = lane_arrows,
    aes(x = start, y = y, xend = end, yend = y),
    arrow = arrow(type = "closed", length = unit(1.35, "mm")),
    linewidth = 0.38, colour = "#7A8793"
  ) +
  # Donor-specific WT networks.
  geom_segment(
    data = network_edges, aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.30, colour = "#8C98A4", lineend = "round"
  ) +
  geom_point(
    data = network_nodes, aes(x = x, y = y),
    shape = 21, size = 2.25, fill = network_nodes$accent,
    colour = "white", stroke = 0.35
  ) +
  geom_point(
    data = lanes, aes(x = 0.39, y = y),
    shape = 21, size = 10.2, fill = "white", colour = lanes$accent, stroke = 0.70
  ) +
  geom_text(
    data = lanes, aes(x = 0.39, y = y + 0.050, label = "WT\nGRN"),
    size = 1.55, lineheight = 0.86, family = FONT, fontface = "bold",
    colour = "#303943"
  ) +
  geom_text(
    data = lanes, aes(x = 0.39, y = y - 0.038, label = "300 genes"),
    size = 1.47, family = FONT, colour = "#59636E"
  ) +
  # Perturbation step.
  geom_rect(
    data = lanes,
    aes(xmin = 0.545, xmax = 0.730, ymin = y - 0.074, ymax = y + 0.074),
    fill = "#FFF2EE", colour = "#D9897F", linewidth = 0.42
  ) +
  geom_text(
    data = lanes, aes(x = 0.604, y = y + 0.023, label = "SQSTM1"),
    hjust = 0, size = 1.87, family = FONT, fontface = "bold", colour = "#A8443E"
  ) +
  geom_text(
    data = lanes, aes(x = 0.604, y = y - 0.027, label = "outgoing row = 0"),
    hjust = 0, size = 1.40, family = FONT, colour = "#5E4E4B"
  ) +
  geom_point(
    data = lanes, aes(x = 0.573, y = y),
    shape = 21, size = 4.8, fill = "#FCE3DF", colour = COL_POS, stroke = 0.55
  ) +
  geom_text(
    data = lanes, aes(x = 0.573, y = y, label = "x"),
    size = 1.72, family = FONT, fontface = "bold", colour = COL_POS
  ) +
  # Ranked output.
  geom_segment(
    data = rank_lines, aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.62, colour = rank_lines$accent, lineend = "round"
  ) +
  geom_point(
    data = rank_lines, aes(x = x, y = y),
    size = 1.15, colour = rank_lines$accent
  ) +
  geom_text(
    data = lanes, aes(x = 0.890, y = y - 0.078, label = "nuclear ranks"),
    size = 1.46, family = FONT, colour = "#59636E"
  ) +
  annotate(
    "text", x = 0.50, y = 0.105,
    label = "Same nuclear feature space  |  independent donor matrices and networks  |  no pooling",
    size = 1.48, family = FONT, fontface = "italic", colour = "#5C6670"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(
    title = "Independent donor-resolved virtual knockout",
    subtitle = "Matched feature space; separate networks, perturbations and readouts"
  ) +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = 7.8, face = "bold", margin = margin(b = 2)),
        plot.subtitle = element_text(size = 5.9, colour = "#4B5563", margin = margin(b = 2)),
        plot.margin = margin(5, 6, 5, 6))

# B: predicted regulation ranks; distances have no direction of change.
vko_dir <- file.path(ROOT, "results", "official_r_vko_manuscript_no_mt_encoded")
leading <- bind_rows(
  fread(file.path(vko_dir, "vko_sqstm1_hoa2_official_r.csv")),
  fread(file.path(vko_dir, "vko_sqstm1_hoa3_official_r.csv"))
) %>%
  filter(Gene != "SQSTM1") %>%
  group_by(donor) %>%
  arrange(rank, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    donor_label = factor(toupper(donor), levels = c("HOA2", "HOA3")),
    gene_key = paste(donor_label, Gene, sep = "::"),
    fdr_hit = adjusted_p_value < 0.05,
    fdr_label = ifelse(
      fdr_hit,
      ifelse(adjusted_p_value < 0.001, "FDR<0.001",
             paste0("FDR=", sprintf("%.3f", adjusted_p_value))),
      NA_character_
    )
  )
gene_levels <- leading %>%
  arrange(donor_label, desc(rank)) %>%
  pull(gene_key) %>%
  unique()
leading$gene_key <- factor(leading$gene_key, levels = gene_levels)
panel_b_key <- data.frame(
  x = NA_real_,
  y = NA_character_,
  role = factor(
    "red rings denote donor-specific FDR<0.05",
    levels = "red rings denote donor-specific FDR<0.05"
  )
)

p_b <- ggplot(leading, aes(x = Z, y = gene_key)) +
  geom_segment(aes(x = 0, xend = Z, yend = gene_key),
               linewidth = 0.45, colour = "#CBD2D9") +
  geom_point(aes(fill = donor_label), shape = 21, size = 2.15,
             colour = "white", stroke = 0.35) +
  geom_point(data = filter(leading, fdr_hit), shape = 21, size = 3.0,
             fill = NA, colour = COL_POS, stroke = 0.70) +
  # FDR-hit key.
  geom_point(
    data = panel_b_key, aes(x = x, y = y, shape = role),
    alpha = 0, size = 0.1, inherit.aes = FALSE, show.legend = TRUE,
    na.rm = TRUE
  ) +
  facet_wrap(~donor_label, nrow = 1, scales = "free_y") +
  scale_y_discrete(
    labels = function(x) sub("^[^:]+::", "", x),
    na.translate = FALSE
  ) +
  scale_x_continuous(limits = c(0, 3.18), breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = c(HOA2 = COL_NEG, HOA3 = COL_TEAL), guide = "none") +
  scale_shape_manual(
    values = c(`red rings denote donor-specific FDR<0.05` = 21),
    name = NULL
  ) +
  guides(
    shape = guide_legend(
      title = "Top 10 downstream ranks per donor;",
      title.position = "left", direction = "horizontal",
      override.aes = list(
        alpha = 1, size = 2.35, colour = COL_POS,
        fill = "white", stroke = 0.70
      )
    )
  ) +
  labs(
    title = "Predicted nuclear differential regulation",
    x = "Standardized WT--KO manifold distance (Z)", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y = element_text(face = "italic", size = 5.2),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    strip.text = element_text(size = 5.9, face = "bold"),
    panel.spacing = unit(4.2, "mm"),
    plot.title.position = "plot",
    legend.location = "plot",
    legend.position = "top", legend.justification = "left",
    legend.direction = "horizontal",
    legend.title = element_text(size = 5.4, face = "plain", colour = "#4B5563"),
    legend.text = element_text(size = 5.4, colour = "#4B5563"),
    legend.key.width = unit(3.1, "mm"), legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(0.35, "mm"), legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(t = -2, r = 0, b = -2, l = -5),
    plot.margin = margin(5, 8, 5, -14)
  )

# C: cross-donor ranks, excluding SQSTM1.
rank2 <- fread(file.path(ROOT, "results", "official_r_vko_manuscript_no_mt_encoded", "vko_sqstm1_hoa2_official_r.csv")) %>%
  filter(Gene != "SQSTM1") %>% select(Gene, rank_hoa2 = rank, fdr_hoa2 = adjusted_p_value)
rank3 <- fread(file.path(ROOT, "results", "official_r_vko_manuscript_no_mt_encoded", "vko_sqstm1_hoa3_official_r.csv")) %>%
  filter(Gene != "SQSTM1") %>% select(Gene, rank_hoa3 = rank, fdr_hoa3 = adjusted_p_value)
ranks <- inner_join(rank2, rank3, by = "Gene") %>%
  mutate(
    hoa3_only_fdr = fdr_hoa3 < 0.05 & fdr_hoa2 >= 0.05,
    shared_top50 = rank_hoa2 <= 50 & rank_hoa3 <= 50,
    label = ifelse(hoa3_only_fdr, Gene, NA_character_)
  )
shared_top50_n <- sum(ranks$shared_top50)
audit <- fread(file.path(ROOT, "results", "official_r_vko_no_mt_cross_donor_audit.csv")) %>%
  filter(profile == "manuscript", model_variant == "mt_excluded") %>% slice(1)
panel_c_key <- data.frame(
  x = NA_real_, y = NA_real_,
  role = factor("Red rings: HOA3-only hits", levels = "Red rings: HOA3-only hits")
)
p_c <- ggplot(ranks, aes(rank_hoa2, rank_hoa3)) +
  annotate(
    "rect", xmin = 0, xmax = 50, ymin = 0, ymax = 50,
    fill = "#EAF2F8", colour = NA, alpha = 0.92
  ) +
  geom_abline(
    intercept = 0, slope = 1, linewidth = 0.32,
    colour = "#AAB4BE"
  ) +
  geom_point(size = 0.70, colour = "#B9C1C9", alpha = 0.55) +
  geom_point(
    data = filter(ranks, shared_top50), size = 0.90,
    colour = "#5B7FA3", alpha = 0.82
  ) +
  geom_point(data = filter(ranks, hoa3_only_fdr), size = 2.1, shape = 21,
             fill = "white", colour = COL_POS, stroke = 0.7) +
  geom_point(
    data = panel_c_key, aes(x = x, y = y, shape = role),
    alpha = 0, size = 0.1, inherit.aes = FALSE, show.legend = TRUE,
    na.rm = TRUE
  ) +
  geom_text_repel(
    data = filter(ranks, hoa3_only_fdr), aes(label = label),
    seed = 20260903, size = 1.8, family = FONT, fontface = "italic",
    box.padding = 0.35, point.padding = 0.18, min.segment.length = 0,
    segment.colour = "#70777D", segment.size = 0.25
  ) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "#AEB7BF", linewidth = 0.30) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "#AEB7BF", linewidth = 0.30) +
  scale_x_reverse(expand = expansion(mult = 0.03)) +
  scale_y_reverse(expand = expansion(mult = 0.03)) +
  scale_shape_manual(
    values = c(`Red rings: HOA3-only hits` = 21),
    name = NULL
  ) +
  guides(
    shape = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        alpha = 1, size = 2.25, colour = COL_POS,
        fill = "white", stroke = 0.70
      )
    )
  ) +
  labs(
    title = "Cross-donor rank agreement",
    subtitle = sprintf(
      "Spearman rho=%.3f; shared top-50=%d; replicated FDR hits=0",
      audit$spearman_rho, shared_top50_n
    ),
    x = "HOA2 perturbation rank (1 = strongest)",
    y = "HOA3 perturbation rank (1 = strongest)"
  ) +
  coord_cartesian() + theme_pub() +
  theme(
    legend.position = "top", legend.justification = "left",
    legend.location = "plot", legend.direction = "horizontal",
    legend.text = element_text(size = 5.35, colour = "#4B5563"),
    legend.key.width = unit(3.1, "mm"), legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(0.35, "mm"), legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(t = -2, r = 0, b = -2, l = -5)
  )

# D: SQSTM1 and 20 matched genes.
cal_dir <- file.path(ROOT, "results", "official_r_vko_matched_controls")
cal <- fread(file.path(cal_dir, "sqstm1_matched_control_calibration.csv")) %>%
  filter(metric == "cross_donor_spearman") %>%
  mutate(
    profile_y = ifelse(profile == "manuscript", 2, 1),
    p_label = paste0("exact p=", sprintf("%.3f", exact_empirical_p))
  )
comparators <- fread(file.path(cal_dir, "matched_control_vko_summary.csv")) %>%
  filter(role == "matched comparator") %>%
  group_by(profile) %>%
  arrange(cross_donor_spearman, .by_group = TRUE) %>%
  mutate(
    profile_y = ifelse(profile == "manuscript", 2, 1),
    y_jitter = profile_y + seq(-0.105, 0.105, length.out = n())
  ) %>%
  ungroup()
calibration_key <- data.frame(
  x = rep(0.20, 4),
  y = rep(0.70, 4),
  role = factor(
    c("Comparators", "CALCOCO2", "Median", "SQSTM1"),
    levels = c("Comparators", "CALCOCO2", "Median", "SQSTM1")
  )
)
p_d <- ggplot() +
  geom_point(
    data = filter(comparators, perturbation_gene != "CALCOCO2"),
    aes(x = cross_donor_spearman, y = y_jitter),
    size = 1.25, colour = "#B8C0C8", alpha = 0.88
  ) +
  geom_point(
    data = filter(comparators, perturbation_gene == "CALCOCO2"),
    aes(x = cross_donor_spearman, y = y_jitter),
    shape = 24, size = 2.35, fill = COL_GOLD, colour = "#815B00", stroke = 0.55
  ) +
  geom_segment(
    data = cal,
    aes(x = comparator_median, xend = comparator_median,
        y = profile_y - 0.16, yend = profile_y + 0.16),
    linewidth = 0.55, colour = "#4B5563"
  ) +
  geom_point(
    data = cal, aes(x = comparator_median, y = profile_y),
    shape = 21, size = 2.25, fill = "white", colour = "#4B5563", stroke = 0.55
  ) +
  geom_point(
    data = cal, aes(x = sqstm1_value, y = profile_y),
    shape = 22, size = 2.55, fill = COL_POS, colour = "#B74B46", stroke = 0.52
  ) +
  # In-panel symbol key.
  geom_point(
    data = calibration_key, aes(x = x, y = y, shape = role),
    alpha = 0, size = 0.1, inherit.aes = FALSE, show.legend = TRUE
  ) +
  scale_shape_manual(
    values = c(Comparators = 16, CALCOCO2 = 24, Median = 21, SQSTM1 = 22),
    name = NULL
  ) +
  guides(
    shape = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(
        alpha = 1,
        size = c(2.0, 2.35, 2.20, 2.20),
        colour = c("#B8C0C8", "#815B00", "#4B5563", "#B74B46"),
        fill = c("#B8C0C8", COL_GOLD, "white", COL_POS),
        stroke = c(0.25, 0.55, 0.55, 0.52)
      )
    )
  ) +
  geom_text(
    data = cal, aes(x = 0.705, y = profile_y, label = p_label),
    hjust = 1, size = 1.52, family = FONT
  ) +
  scale_x_continuous(limits = c(0.20, 0.72), breaks = seq(0.2, 0.7, 0.1)) +
  scale_y_continuous(
    limits = c(0.70, 2.30), breaks = c(1, 2),
    labels = c("Package-default sensitivity", "Primary profile")
  ) +
  labs(
    title = "Receptor-context specificity calibration",
    x = "Cross-donor Spearman rho", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "top", legend.justification = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 5.1, colour = "#4B5563"),
    legend.key.width = unit(3.2, "mm"), legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(0.45, "mm"), legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(t = -2, r = 0, b = -2, l = -5)
  )

# E: matched-null pathway ranks from the nuclear-only refit.
pathways <- fread(file.path(ROOT, "results", "official_r_vko_no_mt_pathway_by_donor.csv")) %>%
  filter(profile == "manuscript", model_variant == "mt_excluded") %>%
  mutate(
    donor = factor(donor, levels = c("hoa2", "hoa3"), labels = c("HOA2", "HOA3")),
    pathway_label = recode(
      pathway,
      Mito_fission = "Mito fission", Mitophagy_core = "Selective clearance",
      Mito_proteostasis = "Mito proteostasis", ROS_defense = "ROS defense",
      mtDNA_release = "mtDNA release", cGAS_STING = "cGAS-STING",
      EC_inflammation = "EC inflammation", YAP_mTOR = "YAP-mTOR"
    ),
    pathway_label = factor(pathway_label, levels = rev(c(
      "cGAS-STING", "OXPHOS", "Mito fission", "Selective clearance", "mtDNA release",
      "Mito proteostasis", "YAP-mTOR", "Angiogenesis", "ROS defense", "EC inflammation"
    ))),
    score = -log10(pmax(empirical_p, 1e-6)),
    cell_label = paste0("p=", sprintf("%.3f", empirical_p)),
    text_colour = ifelse(score > 1.5, "white", "#24282C")
  )
panel_e_key <- data.frame(
  x = NA_character_, y = NA_character_,
  role = factor(
    "Red outline: empirical FDR<0.05 in one donor only; no pathway replicated",
    levels = "Red outline: empirical FDR<0.05 in one donor only; no pathway replicated"
  )
)
p_e <- ggplot(pathways, aes(donor, pathway_label, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_tile(
    data = filter(pathways, empirical_fdr < 0.05),
    fill = NA, colour = COL_POS, linewidth = 0.72
  ) +
  geom_point(
    data = panel_e_key, aes(x = x, y = y, shape = role),
    alpha = 0, size = 0.1, inherit.aes = FALSE, show.legend = TRUE,
    na.rm = TRUE
  ) +
  geom_text(aes(label = cell_label, colour = text_colour), size = 1.55, family = FONT) +
  scale_colour_identity() +
  scale_x_discrete(na.translate = FALSE) +
  scale_y_discrete(na.translate = FALSE) +
  scale_shape_manual(
    values = c(`Red outline: empirical FDR<0.05 in one donor only; no pathway replicated` = 22),
    name = NULL
  ) +
  scale_fill_gradient(low = "#EDF4FA", high = "#2166AC", name = "-log10 empirical p") +
  guides(
    shape = guide_legend(
      position = "top", nrow = 1, byrow = TRUE,
      override.aes = list(
        alpha = 1, size = 2.15, colour = COL_POS,
        fill = "white", stroke = 0.72
      )
    ),
    fill = guide_colourbar(position = "right")
  ) +
  labs(
    title = "Matched-null pathway rank recovery",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "right", legend.justification.top = "left",
    legend.location = "plot",
    legend.text = element_text(size = 5.05, colour = "#4B5563"),
    legend.key.width = unit(3.1, "mm"), legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(0.35, "mm"), legend.box.spacing = unit(0, "pt"),
    legend.margin = margin(t = -2, r = 0, b = -2, l = -5)
  )

top_row <- free(p_a, side = "lr")
middle_row <- (p_b | p_c) + plot_layout(widths = c(1.18, 0.82))
bottom_row <- (p_d | p_e) + plot_layout(widths = c(0.98, 1.02))
figure <- top_row / middle_row / bottom_row +
  plot_layout(heights = c(0.52, 1.05, 0.68)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.tag = element_text(size = 8.8, face = "bold")
  )

save_pub(figure, "Figure6", 183, 185)

writeLines(
  c(
    "Genes Figure 6 QA",
    "Backend: R only.",
    "A: donors are modeled separately; labels clear the network outlines.",
    "B: top 10 nuclear ranks per donor; red rings mark donor-specific FDR hits.",
    sprintf("C: shared top-50=%d; red rings mark HOA3-only hits.", shared_top50_n),
    "No nuclear FDR hit replicated across donors.",
    "D: 20 matched genes, CALCOCO2, comparator medians, and SQSTM1 are shown.",
    "E: red outlines mark donor-specific pathway FDR hits; none replicated.",
    "No experimental knockout observations are represented.",
    "Exports: SVG, PDF, 300-dpi PNG, and 600-dpi TIFF."
  ),
  file.path(QA, "Figure6_genes_revision_QA_notes.txt")
)

cat(normalizePath(file.path(OUT, "Figure6.pdf"), winslash = "/"), "\n")
