#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pROC)
  library(scales)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(flag, default) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value after ", flag)
  args[[index + 1]]
}

root <- normalizePath(
  option_value("--project-root", Sys.getenv("ONFH_ROOT", unset = getwd())),
  winslash = "/", mustWork = TRUE
)
out <- normalizePath(
  option_value("--output-dir", file.path(root, "results", "serum_enhanced_20260910")),
  winslash = "/", mustWork = TRUE
)
formal <- file.path(out, "02_formal")
figdir <- file.path(out, "03_figures")
reportdir <- file.path(out, "04_report")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(reportdir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(name) {
  path <- file.path(formal, name)
  if (!file.exists(path)) stop("Missing formal output: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

agg <- read_required("aggregate_performance.csv")
oof <- read_required("oof_predictions_aggregated.csv")
repeats <- read_required("repeat_performance.csv")
perms <- read_required("label_permutation_null.csv")
perm_sum <- read_required("label_permutation_summary.csv")
random <- read_required("matched_random_space_auc_null.csv")
random_sum <- read_required("matched_random_space_summary.csv")
stability <- read_required("feature_stability_summary.csv")
coef <- read_required("feature_coefficients_long.csv")
fixed <- read_required("fixed_oof_paired_bootstrap_summary.csv")
refit <- read_required("refit_oob_bootstrap_summary.csv")
calibration <- read_required("calibration_bins.csv")
cal_sum <- read_required("calibration_summary.csv")
boot <- read_required("refit_oob_bootstrap_all.csv")
features <- read_required("features_per_outer_fit.csv")
random_qc <- read_required("matched_random_space_qc.csv")

candidate_auc <- agg$auc[agg$model == "candidate"]
ma_auc <- agg$auc[agg$model == "ma"]
roc_candidate <- roc(oof$label, oof$candidate_probability, quiet = TRUE, direction = "<")
roc_ma <- roc(oof$label, oof$ma_probability, quiet = TRUE, direction = "<")
delong <- roc.test(roc_candidate, roc_ma, paired = TRUE, method = "delong")
delong_out <- data.frame(
  comparison = "candidate_vs_ma",
  candidate_auc = as.numeric(auc(roc_candidate)),
  ma_auc = as.numeric(auc(roc_ma)),
  delta_auc = as.numeric(auc(roc_candidate) - auc(roc_ma)),
  paired_delong_p = as.numeric(delong$p.value),
  scope = "paired DeLong test of fixed aggregate cross-fitted predictions"
)
write.csv(delong_out, file.path(formal, "delong_comparison.csv"), row.names = FALSE)

candidate_col <- "#C95B5B"
ma_col <- "#4C78A8"
null_col <- "#B9C7D5"
gold_col <- "#D2A73A"
text_col <- "#26313D"
grid_col <- "#E7ECF0"

theme_review <- theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11, colour = text_col, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 8.2, colour = "#5D6975", margin = margin(b = 6)),
    axis.title = element_text(size = 8.6, colour = text_col),
    axis.text = element_text(size = 7.8, colour = text_col),
    legend.position = "top", legend.title = element_blank(),
    legend.text = element_text(size = 7.6), panel.grid.major.y = element_line(colour = grid_col, linewidth = .3),
    plot.margin = margin(7, 9, 7, 7)
  )

panel_a_data <- perms[, c("task_id", "auc")]
panel_a_data$type <- "label permutation"
write.csv(panel_a_data, file.path(figdir, "Figure7A_label_permutation_plot_data.csv"), row.names = FALSE)
pA <- ggplot(panel_a_data, aes(auc)) +
  geom_histogram(aes(y = after_stat(density)), bins = 32, fill = null_col, colour = "white", linewidth = .25) +
  geom_density(colour = "#708399", linewidth = .65, adjust = 1.1) +
  geom_vline(xintercept = candidate_auc, colour = candidate_col, linewidth = .9) +
  annotate("label", x = candidate_auc, y = Inf, label = sprintf("Observed %.3f\nempirical p=%.4f", candidate_auc, perm_sum$empirical_p[1]),
           hjust = 1.04, vjust = 1.18, size = 2.55, linewidth = .2, colour = text_col, fill = alpha("white", .9)) +
  labs(title = "A  Full-pipeline label permutation", subtitle = sprintf("n=%d; nested CV, probe selection and tuning repeated", nrow(perms)),
       x = "Aggregate repeated-OOF AUC", y = "Density") + theme_review

panel_b_data <- random[, c("task_id", "auc", "average_precision", "brier")]
write.csv(panel_b_data, file.path(figdir, "Figure7B_matched_random_plot_data.csv"), row.names = FALSE)
pB <- ggplot(panel_b_data, aes(auc)) +
  geom_histogram(aes(y = after_stat(density)), bins = 32, fill = "#E3D6AB", colour = "white", linewidth = .25) +
  geom_density(colour = gold_col, linewidth = .7, adjust = 1.1) +
  geom_vline(xintercept = candidate_auc, colour = candidate_col, linewidth = .9) +
  annotate("label", x = candidate_auc, y = Inf, label = sprintf("Candidate %.3f\nspecificity p=%.4f", candidate_auc, random_sum$empirical_p_specificity[1]),
           hjust = 1.04, vjust = 1.18, size = 2.55, linewidth = .2, colour = text_col, fill = alpha("white", .9)) +
  labs(title = "B  Expression/variance-matched spaces", subtitle = "1,000 matched 77-gene spaces;\nouter-training matching only",
       x = "Aggregate repeated-OOF AUC", y = "Density") + theme_review

top15 <- stability[order(-stability$nonzero_frequency, -stability$direction_consistency_given_selected, stability$gene), ][1:15, ]
top15$gene <- factor(top15$gene, levels = rev(top15$gene))
top15$direction <- factor(top15$dominant_direction, levels = c("positive", "negative", "not_selected"))
write.csv(top15, file.path(figdir, "Figure7C_feature_stability_plot_data.csv"), row.names = FALSE)
pC <- ggplot(top15, aes(nonzero_frequency, gene, colour = direction)) +
  geom_segment(aes(x = 0, xend = nonzero_frequency, yend = gene), colour = "#D9E0E6", linewidth = 1.1) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%+.2f", coefficient_median_nonzero)), hjust = -.18, size = 2.25, colour = text_col, na.rm = TRUE) +
  scale_colour_manual(values = c(positive = candidate_col, negative = ma_col, not_selected = "grey70"), labels = c("Positive", "Negative", "Not selected")) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, min(1.12, max(top15$nonzero_frequency) + .18)), expand = c(0, 0)) +
  labs(title = "C  Selection and coefficient direction", subtitle = "Top 15 of 77 genes; label = median non-zero\nstandardized coefficient (red +; blue -)",
       x = "Non-zero frequency across 25 outer fits", y = NULL, colour = NULL) + theme_review +
  theme(panel.grid.major.y = element_blank(), legend.position = "none")

repeat_wide <- reshape(repeats[, c("repeat", "model", "auc", "average_precision", "brier")], idvar = "repeat", timevar = "model", direction = "wide")
metric_long <- rbind(
  data.frame(`repeat` = repeat_wide[["repeat"]], metric = "AUC", Candidate = repeat_wide$auc.candidate, Ma = repeat_wide$auc.ma, check.names = FALSE),
  data.frame(`repeat` = repeat_wide[["repeat"]], metric = "AP", Candidate = repeat_wide$average_precision.candidate, Ma = repeat_wide$average_precision.ma, check.names = FALSE),
  data.frame(`repeat` = repeat_wide[["repeat"]], metric = "Brier score", Candidate = repeat_wide$brier.candidate, Ma = repeat_wide$brier.ma, check.names = FALSE)
)
metric_plot <- rbind(
  data.frame(`repeat` = metric_long[["repeat"]], metric = metric_long$metric, model = "Candidate", value = metric_long$Candidate, check.names = FALSE),
  data.frame(`repeat` = metric_long[["repeat"]], metric = metric_long$metric, model = "Ma comparator", value = metric_long$Ma, check.names = FALSE)
)
write.csv(metric_plot, file.path(figdir, "Figure7D_paired_metrics_plot_data.csv"), row.names = FALSE)
fixed_auc <- fixed[fixed$metric == "delta_auc", ]
refit_auc <- refit[refit$metric == "delta_auc", ]
annotation_d <- sprintf("Aggregate Delta AUC %.3f; paired DeLong p=%.3f\nFixed-OOF interval %.3f to %.3f\nRefit-aware OOB interval %.3f to %.3f",
                        candidate_auc - ma_auc, delong_out$paired_delong_p, fixed_auc$q025, fixed_auc$q975, refit_auc$q025, refit_auc$q975)
pD <- ggplot(metric_plot, aes(model, value, group = .data[["repeat"]])) +
  geom_line(colour = "#AEB8C2", linewidth = .45) +
  geom_point(aes(colour = model), size = 2.1) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(Candidate = candidate_col, `Ma comparator` = ma_col)) +
  labs(title = "D  Paired internal performance", subtitle = annotation_d, x = NULL, y = NULL) + theme_review +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), panel.spacing = unit(10, "pt"), legend.position = "none")

roc_points <- rbind(
  data.frame(model = "Candidate", specificity = roc_candidate$specificities, sensitivity = roc_candidate$sensitivities),
  data.frame(model = "Ma comparator", specificity = roc_ma$specificities, sensitivity = roc_ma$sensitivities)
)
write.csv(roc_points, file.path(figdir, "Figure7E_ROC_plot_data.csv"), row.names = FALSE)
write.csv(calibration, file.path(figdir, "Figure7E_calibration_plot_data.csv"), row.names = FALSE)
pE1 <- ggplot(roc_points, aes(1 - specificity, sensitivity, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") + geom_step(linewidth = .8) +
  coord_equal() + scale_colour_manual(values = c(Candidate = candidate_col, `Ma comparator` = ma_col)) +
  labs(title = "E  ROC and exploratory calibration", subtitle = "Internal cross-fitted predictions; n=40", x = "1 - specificity", y = "Sensitivity") + theme_review +
  theme(legend.position = "bottom")
pE2 <- ggplot(calibration, aes(mean_predicted, observed_fraction, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey70") +
  geom_line(linewidth = .7) + geom_point(aes(size = n), stroke = .2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c(candidate = candidate_col, ma = ma_col), labels = c(candidate = "Candidate", ma = "Ma comparator")) +
  scale_size_continuous(range = c(2.2, 3.5), guide = "none") +
  labs(title = "Four quantile bins", subtitle = sprintf("Brier: candidate %.3f; Ma %.3f", cal_sum$brier[cal_sum$model == "candidate"], cal_sum$brier[cal_sum$model == "ma"]),
       x = "Mean predicted probability", y = "Observed SONFH fraction") + theme_review +
  theme(legend.position = "none")
pE <- pE1 + pE2 + plot_layout(widths = c(1, 1))

combined <- (pA + pB) / (pC + pD + plot_layout(widths = c(1, 1.14))) / pE +
  plot_layout(heights = c(1, 1.20, 1.08))

pdf_path <- file.path(figdir, "Figure7_enhanced_review.pdf")
png_path <- file.path(figdir, "Figure7_enhanced_review.png")
ggsave(pdf_path, combined, width = 190, height = 230, units = "mm", device = cairo_pdf)
agg_png(png_path, width = 190, height = 230, units = "mm", res = 300, background = "white")
print(combined)
dev.off()

boot_plot <- rbind(
  data.frame(metric = "Candidate OOB AUC", value = boot$candidate_auc),
  data.frame(metric = "Candidate - Ma OOB AUC", value = boot$delta_auc)
)
pS1 <- ggplot(boot_plot, aes(value, fill = metric)) + geom_histogram(bins = 32, colour = "white") +
  facet_wrap(~metric, scales = "free_x", ncol = 1) + guides(fill = "none") + labs(title = "Refit-aware OOB distributions", x = NULL, y = "Count") + theme_review
pS2 <- ggplot(features, aes(n_selected_genes)) + geom_histogram(binwidth = 1, fill = candidate_col, colour = "white") +
  labs(title = "Selected genes per outer fit", x = "Number of non-zero coefficients", y = "Count") + theme_review
pS3 <- ggplot(random_qc, aes(expanded_fraction)) + geom_histogram(bins = 25, fill = gold_col, colour = "white") +
  labs(title = "Matched-space expansion", x = "Fraction requiring Manhattan expansion", y = "Count") + theme_review
pS4 <- ggplot(calibration, aes(factor(bin), n, fill = model)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(candidate = candidate_col, ma = ma_col)) +
  labs(title = "Samples per calibration bin", x = "Quantile bin", y = "n") + theme_review
supp <- (pS1 | pS2) / (pS3 | pS4) + plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave(file.path(figdir, "Figure7_enhanced_supplemental_diagnostics.pdf"), supp, width = 190, height = 150, units = "mm", device = cairo_pdf)
agg_png(file.path(figdir, "Figure7_enhanced_supplemental_diagnostics.png"), width = 190, height = 150, units = "mm", res = 300, background = "white")
print(supp)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(reportdir, "R_sessionInfo.txt"))
cat("Figure 7 enhanced review and DeLong output generated.\n")
