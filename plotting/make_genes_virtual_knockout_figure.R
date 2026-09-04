#!/usr/bin/env Rscript

# Genes-oriented Figure 6. The primary visualization uses the refit that
# excludes mitochondrially encoded genes. Core conclusion: donor-separated
# SQSTM1 in silico perturbation produces moderately concordant nuclear-gene
# ranks, but no replicated nuclear FDR hit and no specificity over matched
# comparator perturbations. CALCOCO2 is displayed as a post-analysis secondary candidate
# within the comparator distribution, not as a validated second target.
# Evidence chain: A, design; B, leading predicted
# nuclear differential regulation; C, global donor reproducibility; D,
# matched-gene specificity calibration; E, pathway-level consequence.
# Archetype: schematic-led quantitative composite. No experimental knockout
# observations are represented in this figure.

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

# A: two explicit donor lanes; networks and perturbations were never pooled.
input_boxes <- data.frame(
  xmin = 0.03, xmax = 0.30,
  ymin = c(0.57, 0.26), ymax = c(0.79, 0.48),
  fill = c("#DCEAF5", "#DDF1EA"),
  label = c(
    "HOA2 control ECs\nn=1,043\nSQSTM1+=70.8%",
    "HOA3 control ECs\nn=760\nSQSTM1+=72.2%"
  )
)
analysis_boxes <- data.frame(
  xmin = 0.40, xmax = 0.98,
  ymin = c(0.57, 0.26), ymax = c(0.79, 0.48),
  label = c(
    "Independent WT GRN\n300 shared nuclear genes\nSQSTM1 row=0  ->  WT--KO alignment\nDonor-specific nuclear-gene ranks",
    "Independent WT GRN\n300 shared nuclear genes\nSQSTM1 row=0  ->  WT--KO alignment\nDonor-specific nuclear-gene ranks"
  )
)
arrows <- data.frame(y = c(0.68, 0.37))
p_a <- ggplot() +
  geom_rect(
    data = input_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = input_boxes$fill, colour = "#536273", linewidth = 0.38
  ) +
  geom_rect(
    data = analysis_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#F7EDDF", colour = "#536273", linewidth = 0.38
  ) +
  geom_text(
    data = input_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 1.55, lineheight = 0.94, family = FONT
  ) +
  geom_text(
    data = analysis_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    size = 1.47, lineheight = 0.94, family = FONT
  ) +
  geom_segment(
    data = arrows,
    aes(x = 0.305, y = y, xend = 0.395, yend = y),
    arrow = arrow(type = "closed", length = unit(1.55, "mm")),
    linewidth = 0.38, colour = "#536273"
  ) +
  annotate(
    "text", x = 0.50, y = 0.105,
    label = paste0(
      "Shared nuclear feature universe; donor matrices and networks were not pooled\n",
      "HOA1 ineligible (26 ECs)  |  two parameter profiles  |  20 matched-gene perturbations  |  no experimental KO"
    ),
    size = 1.32, lineheight = 1.04, family = FONT, colour = "#665C34"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(
    title = "Donor-separated in silico perturbation design",
    subtitle = "HOA2 and HOA3 networks were fitted and perturbed independently"
  ) +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(size = 7.4, face = "bold", margin = margin(b = 2)),
        plot.subtitle = element_text(size = 5.5, colour = "#4B5563", margin = margin(b = 2)),
        plot.margin = margin(5, 7, 5, 7))

# B: the primary perturbation output, not the input WT topology. Distances are
# non-directional; positive/negative expression change is not inferred.
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

p_b <- ggplot(leading, aes(x = Z, y = gene_key)) +
  geom_segment(aes(x = 0, xend = Z, yend = gene_key),
               linewidth = 0.45, colour = "#CBD2D9") +
  geom_point(aes(fill = donor_label), shape = 21, size = 2.15,
             colour = "white", stroke = 0.35) +
  geom_point(data = filter(leading, fdr_hit), shape = 21, size = 3.0,
             fill = NA, colour = COL_POS, stroke = 0.70) +
  facet_wrap(~donor_label, nrow = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^:]+::", "", x)) +
  scale_x_continuous(limits = c(0, 3.18), breaks = c(0, 1, 2, 3),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = c(HOA2 = COL_NEG, HOA3 = COL_TEAL), guide = "none") +
  labs(
    title = "Predicted nuclear differential regulation",
    subtitle = "Top 10 downstream ranks per donor; red rings denote donor-specific FDR<0.05",
    x = "Standardized WT--KO manifold distance (Z)", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y = element_text(face = "italic", size = 5.2),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    strip.text = element_text(size = 5.9, face = "bold"),
    panel.spacing = unit(4.2, "mm"),
    plot.margin = margin(5, 10, 5, 10)
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
  scale_x_reverse(expand = expansion(mult = 0.03)) +
  scale_y_reverse(expand = expansion(mult = 0.03)) +
  labs(
    title = "Cross-donor rank reproducibility",
    subtitle = sprintf(
      "Spearman rho=%.3f; replicated FDR hits=0\nRed rings: HOA3-only hits",
      audit$spearman_rho
    ),
    x = "HOA2 perturbation rank (1 = strongest)",
    y = "HOA3 perturbation rank (1 = strongest)"
  ) +
  coord_cartesian() + theme_pub()

# D: SQSTM1 against the full distribution of 20 matched comparator genes.
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
    shape = 18, size = 3.0, colour = COL_POS
  ) +
  geom_text(
    data = cal, aes(x = 0.655, y = profile_y, label = p_label),
    hjust = 0, size = 1.60, family = FONT
  ) +
  scale_x_continuous(limits = c(0.20, 0.72), breaks = seq(0.2, 0.7, 0.1)) +
  scale_y_continuous(
    limits = c(0.70, 2.30), breaks = c(1, 2),
    labels = c("Package-default sensitivity", "Primary profile")
  ) +
  labs(
    title = "Receptor-context specificity calibration",
    subtitle = "Grey: matched genes; gold triangle: CALCOCO2; open circle: median; red diamond: SQSTM1",
    x = "Cross-donor Spearman rho", y = NULL
  ) +
  theme_pub() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# E: matched-null pathway ranks for the nuclear-only primary refit.
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
p_e <- ggplot(pathways, aes(donor, pathway_label, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_tile(
    data = filter(pathways, empirical_fdr < 0.05),
    fill = NA, colour = COL_POS, linewidth = 0.72
  ) +
  geom_text(aes(label = cell_label, colour = text_colour), size = 1.55, family = FONT) +
  scale_colour_identity() +
  scale_fill_gradient(low = "#EDF4FA", high = "#2166AC", name = "-log10 empirical p") +
  labs(
    title = "Matched-null pathway rank recovery",
    subtitle = "Red outline: empirical FDR<0.05 in one donor only; no pathway replicated",
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
    "Core finding: selective-autophagy receptor evidence is heterogeneous; after excluding mitochondrially encoded features, SQSTM1 perturbation ranks are donor-correlated but not specific relative to matched comparator genes.",
    "Panel B displays the top ten nuclear downstream differential-regulation ranks per donor; the WT topology network is not used as the primary result display.",
    "HOA3-only FDR genes (VWF, NFKBIA and C7) are not described as cross-donor hits.",
    "No nuclear downstream gene reached FDR<0.05 in both donors.",
    "One pathway reached empirical FDR<0.05 in HOA3 only; no pathway replicated.",
    "Matched-comparator exact empirical p values are shown for both parameter profiles, with all 20 comparator values displayed and CALCOCO2 identified as a post-analysis secondary candidate inside the comparator distribution.",
    "Pathway cells with donor-specific empirical FDR<0.05 are outlined; no pathway replicated.",
    "No experimental knockout observations are represented."
  ),
  file.path(QA, "Figure6_genes_revision_QA_notes.txt")
)

cat(normalizePath(file.path(OUT, "Figure6.pdf"), winslash = "/"), "\n")
