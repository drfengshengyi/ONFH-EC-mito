#!/usr/bin/env Rscript

# Genes-oriented Figure 6. The primary visualization uses the refit that
# excludes mitochondrially encoded genes. It emphasizes calibration and
# non-specificity rather than a positive MT-ND1/SQSTM1 narrative.

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

# A: explicit audit workflow.
boxes <- data.frame(
  xmin = 0.08, xmax = 0.92,
  ymin = c(0.69, 0.40, 0.11), ymax = c(0.88, 0.59, 0.30),
  fill = c("#DCEAF5", "#E4F2E8", "#FBE6D8"),
  label = c(
    "Control ECs: HOA2 n=1,043; HOA3 n=760",
    "scTenifoldKnk v1.1; 300 shared nuclear genes",
    "Virtual SQSTM1 KO\nzero outgoing row; WT--KO alignment"
  )
)
p_a <- ggplot() +
  geom_rect(
    data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = boxes$fill, colour = "#536273", linewidth = 0.45
  ) +
  geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
            size = 1.85, lineheight = 0.95, family = FONT) +
  geom_segment(aes(x = 0.50, y = 0.685, xend = 0.50, yend = 0.60),
               arrow = arrow(type = "closed", length = unit(2.2, "mm")),
               linewidth = 0.45, colour = "#536273") +
  geom_segment(aes(x = 0.50, y = 0.395, xend = 0.50, yend = 0.31),
               arrow = arrow(type = "closed", length = unit(2.2, "mm")),
               linewidth = 0.45, colour = "#536273") +
  annotate(
    "text", x = 0.50, y = 0.015,
    label = "Primary: 20-network q=0.95 profile; sensitivity: package-default profile",
    size = 1.60, family = FONT, colour = "#665C34"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "Prespecified, donor-separated virtual-knockout design") +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 4)),
        plot.margin = margin(5, 7, 5, 7))

# B: nuclear-only WT neighborhood, averaged only for display.
edge_dir <- file.path(ROOT, "results", "official_r_vko_figure_data_no_mt_encoded")
edge2 <- fread(file.path(edge_dir, "sqstm1_outgoing_edges_hoa2.csv")) %>%
  select(gene, w2 = wt_outgoing_weight)
edge3 <- fread(file.path(edge_dir, "sqstm1_outgoing_edges_hoa3.csv")) %>%
  select(gene, w3 = wt_outgoing_weight)
network <- full_join(edge2, edge3, by = "gene") %>%
  mutate(across(c(w2, w3), ~replace_na(.x, 0)),
         mean_abs = (abs(w2) + abs(w3)) / 2,
         mean_weight = (w2 + w3) / 2) %>%
  arrange(desc(mean_abs)) %>%
  slice_head(n = 18) %>%
  mutate(
    theta = seq(0, 2 * pi - 2 * pi / n(), length.out = n()),
    x = 1.24 * cos(theta), y = sin(theta),
    lx = 1.42 * cos(theta), ly = 1.12 * sin(theta),
    hjust = ifelse(cos(theta) > 0.15, 0, ifelse(cos(theta) < -0.15, 1, 0.5)),
    sign = ifelse(mean_weight >= 0, "Positive", "Negative")
  )
p_b <- ggplot(network) +
  geom_segment(
    aes(x = 0, y = 0, xend = 0.91 * x, yend = 0.91 * y, colour = sign,
        linewidth = mean_abs),
    arrow = arrow(type = "closed", length = unit(1.2, "mm")), alpha = 0.72
  ) +
  geom_point(aes(x, y, fill = sign, size = mean_abs), shape = 21, colour = "white", stroke = 0.35) +
  geom_text(aes(lx, ly, label = gene, hjust = hjust), size = 1.63, family = FONT,
            fontface = "italic") +
  annotate("point", x = 0, y = 0, shape = 21, size = 9.0, fill = "#F5D8D3",
           colour = "#A63D35", stroke = 0.7) +
  annotate("text", x = 0, y = 0, label = "SQSTM1\nWT source", size = 2.1,
           family = FONT, fontface = "bold", lineheight = 0.88) +
  scale_colour_manual(values = c(Positive = COL_POS, Negative = COL_NEG), name = "Edge sign") +
  scale_fill_manual(values = c(Positive = COL_POS, Negative = COL_NEG), guide = "none") +
  scale_linewidth_continuous(range = c(0.25, 0.85), guide = "none") +
  scale_size_continuous(range = c(1.8, 3.2), guide = "none") +
  coord_equal(xlim = c(-1.62, 1.62), ylim = c(-1.25, 1.25), clip = "off") +
  labs(
    title = "Top nuclear WT outgoing edges",
    subtitle = "Mean absolute edge weight across HOA2 and HOA3; topology is descriptive"
  ) +
  theme_void(base_family = FONT) +
  theme(
    plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 5.5, colour = "#4B5563"),
    legend.position = "bottom", legend.title = element_text(size = 5.5, face = "bold"),
    legend.text = element_text(size = 5.2), plot.margin = margin(5, 10, 5, 10)
  )

# C: direct cross-donor rank agreement, excluding the perturbed target itself.
rank2 <- fread(file.path(ROOT, "results", "official_r_vko_manuscript_no_mt_encoded", "vko_sqstm1_hoa2_official_r.csv")) %>%
  filter(Gene != "SQSTM1") %>% select(Gene, rank_hoa2 = rank, fdr_hoa2 = adjusted_p_value)
rank3 <- fread(file.path(ROOT, "results", "official_r_vko_manuscript_no_mt_encoded", "vko_sqstm1_hoa3_official_r.csv")) %>%
  filter(Gene != "SQSTM1") %>% select(Gene, rank_hoa3 = rank, fdr_hoa3 = adjusted_p_value)
ranks <- inner_join(rank2, rank3, by = "Gene") %>%
  mutate(hoa3_only_fdr = fdr_hoa3 < 0.05 & fdr_hoa2 >= 0.05,
         label = ifelse(hoa3_only_fdr, Gene, NA_character_))
audit <- fread(file.path(ROOT, "results", "official_r_vko_no_mt_cross_donor_audit.csv")) %>%
  filter(profile == "manuscript", model_variant == "mt_excluded") %>% slice(1)
p_c <- ggplot(ranks, aes(rank_hoa2, rank_hoa3)) +
  geom_point(size = 0.75, colour = "#AEB7BF", alpha = 0.65) +
  geom_point(data = filter(ranks, hoa3_only_fdr), size = 2.1, shape = 21,
             fill = "white", colour = COL_POS, stroke = 0.7) +
  geom_text_repel(
    data = filter(ranks, hoa3_only_fdr), aes(label = label),
    seed = 20260903, size = 1.8, family = FONT, fontface = "italic",
    box.padding = 0.35, point.padding = 0.18, min.segment.length = 0,
    segment.colour = "#70777D", segment.size = 0.25
  ) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = COL_GREY, linewidth = 0.35) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = COL_GREY, linewidth = 0.35) +
  scale_x_reverse() + scale_y_reverse() +
  annotate(
    "label", x = Inf, y = -Inf,
    label = sprintf("Spearman rho = %.3f\nreplicated nuclear FDR hits = 0", audit$spearman_rho),
    hjust = 1.04, vjust = -0.12, size = 1.7, family = FONT,
    fill = "white", colour = "#333333", linewidth = 0.18
  ) +
  labs(
    title = "Cross-donor nuclear rank agreement",
    subtitle = "Donor-correlated but non-specific; outlines: HOA3 FDR<0.05 only",
    x = "HOA2 perturbation rank (1 = strongest)",
    y = "HOA3 perturbation rank (1 = strongest)"
  ) +
  coord_equal() + theme_pub()

# D: SQSTM1 against matched comparator genes.
cal <- fread(file.path(ROOT, "results", "official_r_vko_matched_controls", "sqstm1_matched_control_calibration.csv")) %>%
  filter(metric == "cross_donor_spearman") %>%
  mutate(
    profile_label = factor(profile, levels = c("official_default", "manuscript"),
                           labels = c("Package-default sensitivity", "Primary profile")),
    p_label = paste0("exact p=", sprintf("%.3f", exact_empirical_p))
  )
p_d <- ggplot(cal, aes(comparator_median, profile_label)) +
  geom_segment(aes(x = comparator_min, xend = comparator_max, yend = profile_label),
               linewidth = 1.8, lineend = "round", colour = "#C4CBD2") +
  geom_point(shape = 21, size = 2.2, fill = "white", colour = "#5C6670", stroke = 0.55) +
  geom_point(aes(x = sqstm1_value), shape = 18, size = 3.0, colour = COL_POS) +
  geom_text(aes(x = comparator_max, label = p_label), hjust = -0.08, size = 1.65, family = FONT) +
  scale_x_continuous(limits = c(0.20, 0.72), breaks = seq(0.2, 0.7, 0.1)) +
  labs(
    title = "Matched-gene calibration",
    subtitle = "Grey range and circle: comparator range and median; diamond: SQSTM1",
    x = "Cross-donor Spearman rho", y = NULL
  ) +
  theme_pub()

# E: matched-null pathway ranks for the nuclear-only primary refit.
pathways <- fread(file.path(ROOT, "results", "official_r_vko_no_mt_pathway_by_donor.csv")) %>%
  filter(profile == "manuscript", model_variant == "mt_excluded") %>%
  mutate(
    donor = factor(donor, levels = c("hoa2", "hoa3"), labels = c("HOA2", "HOA3")),
    pathway_label = recode(
      pathway,
      Mito_fission = "Mito fission", Mitophagy_core = "Mitophagy core",
      Mito_proteostasis = "Mito proteostasis", ROS_defense = "ROS defense",
      mtDNA_release = "mtDNA release", cGAS_STING = "cGAS-STING",
      EC_inflammation = "EC inflammation", YAP_mTOR = "YAP-mTOR"
    ),
    pathway_label = factor(pathway_label, levels = rev(c(
      "cGAS-STING", "OXPHOS", "Mito fission", "Mitophagy core", "mtDNA release",
      "Mito proteostasis", "YAP-mTOR", "Angiogenesis", "ROS defense", "EC inflammation"
    ))),
    score = -log10(pmax(empirical_p, 1e-6)),
    cell_label = ifelse(empirical_fdr < 0.05,
                        paste0("p=", sprintf("%.3f", empirical_p), "*"),
                        paste0("p=", sprintf("%.3f", empirical_p)))
  )
p_e <- ggplot(pathways, aes(donor, pathway_label, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label), size = 1.55, family = FONT,
            colour = ifelse(pathways$score > 1.5, "white", "#24282C")) +
  scale_fill_gradient(low = "#EDF4FA", high = "#2166AC", name = "-log10 p") +
  labs(
    title = "Matched-null pathway rank recovery",
    subtitle = "* empirical FDR<0.05 in one donor only; no pathway replicated",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")

right_lower <- p_d / p_e + plot_layout(heights = c(0.42, 1.18))
figure <- ((p_a | p_b) + plot_layout(widths = c(0.95, 1.20))) /
  ((p_c | right_lower) + plot_layout(widths = c(0.95, 1.20))) +
  plot_layout(heights = c(0.80, 1.20)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_pub(figure, "Figure6", 183, 170)

writeLines(
  c(
    "Genes Figure 6 QA contract",
    "Backend: R only.",
    "Core finding: after excluding mitochondrially encoded features, SQSTM1 perturbation ranks are donor-correlated but not specific relative to matched comparator genes.",
    "The network panel contains nuclear genes only and is explicitly descriptive.",
    "HOA3-only FDR genes (VWF, NFKBIA and C7) are not described as cross-donor hits.",
    "No nuclear downstream gene reached FDR<0.05 in both donors.",
    "One pathway reached empirical FDR<0.05 in HOA3 only; no pathway replicated.",
    "Matched-comparator exact empirical p values are shown for both parameter profiles."
  ),
  file.path(QA, "Figure6_genes_revision_QA_notes.txt")
)

cat(normalizePath(file.path(OUT, "Figure6.pdf"), winslash = "/"), "\n")
