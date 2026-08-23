# Stage 1a v4: explicit library/participant metadata and inference guardrails.
# Run with: Rscript analysis/prepare_matrices.R from any working directory.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "analysis/prepare_matrices.R"
default_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
root <- Sys.getenv("ONFH_ROOT", unset = default_root)
logf <- file.path(root, "analysis", "stage1a_v4.log")
if (file.exists(logf)) file.remove(logf)
lg <- function(...) {
  s <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(s, "\n"); cat(s, "\n", file = logf, append = TRUE); flush.console()
}

sample_meta <- read.csv(
  file.path(root, "analysis", "sample_metadata_v4.csv"),
  stringsAsFactors = FALSE,
  na.strings = c("")
)
stopifnot(!anyDuplicated(sample_meta$sample))
stopifnot(sum(sample_meta$dataset == "sonfh_cystic") == 4)
stopifnot(all(!sample_meta$independent_for_inference[sample_meta$dataset == "sonfh_cystic"]))

ALIAS <- c("STING1" = "TMEM173")
harmonize <- function(m) {
  rn <- rownames(m)
  hit <- rn %in% names(ALIAS)
  if (any(hit)) {
    rn[hit] <- ALIAS[rn[hit]]
    rownames(m) <- rn
    if (anyDuplicated(rownames(m))) {
      m <- as(m, "dgTMatrix")
      m <- Matrix.utils::aggregate.Matrix(as(m, "dgCMatrix"), groupings = rownames(m), fun = "sum")
    }
  }
  m
}

mats <- list()
for (s in sample_meta$sample[sample_meta$dataset == "liao_alcohol"]) {
  m <- readRDS(file.path(root, "data", "liao2022", paste0(s, ".rds")))
  colnames(m) <- paste0(s, "_", colnames(m))
  mats[[s]] <- harmonize(m); lg("rds", s, ncol(m))
}
gse1 <- file.path(root, "data", "gse169396")
gsm1 <- c(S1 = "GSM5201883_S1", S2 = "GSM5201884_S2", S3 = "GSM5201885_S3", S4 = "GSM5201886_S4")
for (s in names(gsm1)) {
  g <- gsm1[[s]]
  m <- ReadMtx(
    mtx = file.path(gse1, paste0(g, "_matrix.mtx.gz")),
    cells = file.path(gse1, paste0(g, "_barcodes.tsv.gz")),
    features = file.path(gse1, paste0(g, "_features.tsv.gz")), feature.column = 2
  )
  colnames(m) <- paste0(s, "_", colnames(m))
  mats[[s]] <- harmonize(m); lg("mtx", s, ncol(m))
}
gse2 <- file.path(root, "data", "gse290411")
gsm2 <- c(AL = "GSM8812280_AL", bone = "GSM8812281_bone", DWR = "GSM8812282_DWR", WXJ = "GSM8812283_WXJ")
for (s in names(gsm2)) {
  g <- gsm2[[s]]
  m <- ReadMtx(
    mtx = file.path(gse2, paste0(g, "matrix.mtx.gz")),
    cells = file.path(gse2, paste0(g, "barcodes.tsv.gz")),
    features = file.path(gse2, paste0(g, "features.tsv.gz")), feature.column = 2
  )
  colnames(m) <- paste0(s, "_", colnames(m))
  mats[[s]] <- harmonize(m); lg("mtx", s, ncol(m))
}

genes <- Reduce(intersect, lapply(mats, rownames))
mats <- lapply(mats, function(m) m[genes, , drop = FALSE])
big <- do.call(cbind, mats)
rm(mats); gc(verbose = FALSE)

nFeature <- diff(big@p)
nCount <- Matrix::colSums(big)
mt_rows <- grep("^MT-", rownames(big))
hb_rows <- match(c("HBA1", "HBA2", "HBB"), rownames(big)); hb_rows <- hb_rows[!is.na(hb_rows)]
pmt <- Matrix::colSums(big[mt_rows, , drop = FALSE]) / nCount * 100
phb <- Matrix::colSums(big[hb_rows, , drop = FALSE]) / nCount * 100
keep <- nFeature >= 500 & nFeature <= 5000 & pmt < 10 & phb < 1
big2 <- big[, keep, drop = FALSE]

cell_meta <- data.frame(
  cell = colnames(big2),
  sample = sub("_.*$", "", colnames(big2)),
  stringsAsFactors = FALSE
)
idx <- match(cell_meta$sample, sample_meta$sample)
for (v in c("dataset", "group", "participant_id", "inferential_unit", "independent_for_inference", "source_note")) {
  cell_meta[[v]] <- sample_meta[[v]][idx]
}
cell_meta$nFeature <- nFeature[keep]
cell_meta$nCount <- nCount[keep]
cell_meta$percent.mt <- pmt[keep]
cell_meta$percent.hb <- phb[keep]

saveRDS(big2, file.path(root, "analysis", "counts_qc_v4.rds"), compress = FALSE)
saveRDS(cell_meta, file.path(root, "analysis", "cell_meta_v4.rds"), compress = FALSE)
lg("saved v4 matrices; 19 libraries, SONFH inferential units disabled pending participant mapping")
