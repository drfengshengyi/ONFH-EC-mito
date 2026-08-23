# GSEA on pseudobulk DE ranks: fgsea + msigdbr (Hallmark, GO BP, Reactome)
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
})
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "analysis/pathway_enrichment.R"
default_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- Sys.getenv("ONFH_ROOT", unset = default_root)
logf <- file.path(root, "analysis/r_gsea.log")
if (file.exists(logf)) file.remove(logf)
lg <- function(...) { s <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = " ")); cat(s, "\n"); cat(s, "\n", file = logf, append = TRUE); flush.console() }

lg("fetching msigdbr gene sets")
lg("enrichment null: preranked competitive gene-set null; participant labels are not permuted")
msig <- tryCatch(msigdbr(species = "Homo sapiens"),
                 error = function(e) { lg("msigdbr error:", conditionMessage(e)); NULL })
if (is.null(msig)) stop("msigdbr unavailable")

hallmark <- msig[msig$gs_collection == "H", ]
gobp <- msig[msig$gs_collection == "C5" & msig$gs_subcollection == "GO:BP", ]
react <- msig[msig$gs_collection == "C2" & msig$gs_subcollection == "CP:REACTOME", ]
mk_list <- function(df) split(df$gene_symbol, df$gs_name)
gs_all <- list(H = mk_list(hallmark), GOBP = mk_list(gobp), REACTOME = mk_list(react))
lg("genesets:", sapply(gs_all, length))

run_fgsea <- function(de_csv, tag) {
  de <- fread(de_csv)
  de <- de[!is.na(de$stat) & !is.na(de$pvalue), ]
  ranks <- de$stat
  names(ranks) <- de$V1
  if ("Row.names" %in% names(de)) names(ranks) <- de$Row.names
  if (all(is.na(names(ranks)))) names(ranks) <- de[[1]]
  ranks <- ranks[!is.na(names(ranks))]
  ranks <- ranks[!duplicated(names(ranks))]
  for (db in names(gs_all)) {
    set.seed(42)
    fg <- fgseaMultilevel(pathways = gs_all[[db]], stats = ranks,
                         minSize = 10, maxSize = 500, eps = 0)
    fg <- fg[order(fg$padj), ]
    out <- file.path(root, "analysis", paste0("gsea_", tag, "_", db, ".csv"))
    fwrite(as.data.frame(fg)[, c("pathway", "pval", "padj", "NES", "size")], out)
    n_sig <- sum(fg$padj < 0.05, na.rm = TRUE)
    lg(tag, db, "sig pathways:", n_sig)
    top <- head(as.data.frame(fg)[, c("pathway", "NES", "padj")], 8)
    print(paste("==", tag, db, "=="))
    print(top)
  }
}

# Only the participant-mapped within-Liao contrast is used for enrichment.
# The cross-cohort SONFH effect table is descriptive and is not submitted to
# participant-level pathway inference.
run_fgsea(file.path(root, "analysis/de_ec_ONFH_3A_vs_HOA_v4.csv"), "ONFH3A_vs_HOA")
lg("DONE")

fwrite(data.frame(
  method = "fgseaMultilevel preranked competitive enrichment",
  ranking_statistic = "participant-pseudobulk DESeq2 Wald statistic",
  null_model = "random gene-set membership conditional on the ranked gene list",
  participant_label_permutation = FALSE,
  interpretive_scope = "hypothesis-generating gene-set prioritization; not independent participant-level validation",
  fgsea_version = as.character(packageVersion("fgsea")),
  msigdbr_version = as.character(packageVersion("msigdbr"))
), file.path(root, "analysis/gsea_method_v7.csv"))
