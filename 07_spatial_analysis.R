# 07_spatial_analysis.R
# Purpose : Spatial transcriptomics analysis of DRI (GSE250521)
# Input   : data/GSE250521/ (Visium files: *_visium_matrix.mtx.gz,
#           *_visium_barcodes.tsv.gz, *_visium_features.tsv.gz,
#           *_visium_tissue_positions_list.csv.gz)
#           output/datExpr_TCGA.rds, output/TDS_score.rds
#           (used to recompute the DRI signature gene list)
# Output  : output/TableS8_spatial_summary.csv
#           output/*_spatial_scores.rds
# Depends : Seurat, ape, Matrix, dplyr

suppressPackageStartupMessages({
  library(Seurat); library(ape); library(Matrix); library(dplyr)
})

data_dir <- "data/GSE250521"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

# Recompute signature gene list (consistent with 01_compute_DRI_TCGA.R)
expr    <- readRDS("output/datExpr_TCGA.rds")
tds_all <- readRDS("output/TDS_score.rds")
keep    <- apply(expr, 2, median) > 1
expr    <- expr[, keep]
nsun2_cor <- cor(expr, expr[, "NSUN2"], method = "spearman")[, 1]
tds_cor   <- cor(expr, tds_all, method = "spearman")[, 1]
core <- names(which(nsun2_cor > 0.15 & tds_cor < -0.10))
supp <- setdiff(names(which(tds_cor < -0.15)), core)
dri_genes <- union(core, supp)
cat("DRI signature genes:", length(dri_genes), "\n")

tds_genes <- c("TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1")

mat_files  <- list.files(data_dir, pattern = "_visium_matrix.mtx.gz$",
                         full.names = TRUE)
sample_ids <- sub("_visium_matrix.mtx.gz$", "", basename(mat_files))
cat("Visium samples found:", length(sample_ids), "\n")

results <- data.frame()

for (sid in sample_ids) {
  cat("\n==========", sid, "==========\n")

  bc <- readLines(gzfile(file.path(data_dir,
              paste0(sid, "_visium_barcodes.tsv.gz"))))

  counts <- ReadMtx(
    mtx      = file.path(data_dir, paste0(sid, "_visium_matrix.mtx.gz")),
    cells    = file.path(data_dir, paste0(sid, "_visium_barcodes.tsv.gz")),
    features = file.path(data_dir, paste0(sid, "_visium_features.tsv.gz")),
    feature.column = 2, cell.column = 1
  )
  colnames(counts) <- bc

  pos <- read.csv(gzfile(file.path(data_dir,
                paste0(sid, "_visium_tissue_positions_list.csv.gz"))),
                header = FALSE)
  colnames(pos) <- c("barcode", "in_tissue", "array_row", "array_col",
                     "pxl_row", "pxl_col")
  pos <- pos[pos$barcode %in% bc, ]

  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj, verbose = FALSE)
  colnames(obj) <- bc

  nCount   <- colSums(counts)
  nFeature <- colSums(counts > 0)

  dri_avail <- intersect(dri_genes, rownames(obj))
  tds_avail <- intersect(tds_genes, rownames(obj))

  if (length(dri_avail) < 500 || length(tds_avail) < 3) {
    cat("Insufficient gene coverage, skipping\n")
    next
  }

  # AddModuleScore to remove sequencing-depth confounding
  set.seed(123)
  obj <- AddModuleScore(obj, features = list(dri = dri_avail),
                        name = "DRI_score", nbin = 24, ctrl = 100)
  obj <- AddModuleScore(obj, features = list(tds = tds_avail),
                        name = "TDS_score", nbin = 24, ctrl = 100)
  dri <- obj$DRI_score1
  tds <- obj$TDS_score1
  names(dri) <- colnames(obj)
  names(tds) <- colnames(obj)

  common <- intersect(names(dri), pos$barcode)
  common <- common[is.finite(dri[common]) & is.finite(tds[common])]
  if (length(common) < 50) next

  # Raw correlation
  ct_raw <- cor.test(dri[common], tds[common], method = "spearman")

  # Partial correlation controlling nCount + nFeature
  d <- data.frame(dri = dri[common], tds = tds[common],
                  nCount = nCount[common], nFeature = nFeature[common])
  res_dri <- residuals(lm(dri ~ nCount + nFeature, data = d))
  res_tds <- residuals(lm(tds ~ nCount + nFeature, data = d))
  ct_par  <- cor.test(res_dri, res_tds, method = "spearman")

  # Moran's I (kNN, k = 6)
  pos_sub <- pos[match(common, pos$barcode), ]
  coords  <- as.matrix(pos_sub[, c("pxl_row", "pxl_col")])
  k <- min(6, nrow(coords) - 1)
  W <- matrix(0, nrow(coords), nrow(coords))
  for (i in seq_len(nrow(coords))) {
    dd <- sqrt(rowSums(sweep(coords, 2, coords[i, ])^2))
    dd[i] <- Inf
    W[i, order(dd)[1:k]] <- 1
  }
  W <- (W + t(W)) / 2
  moran <- Moran.I(res_dri, W)

  cat("raw rho =", round(ct_raw$estimate, 3),
      "| partial rho =", round(ct_par$estimate, 3),
      "| Moran's I =", round(moran$observed, 3), "\n")

  results <- rbind(results, data.frame(
    Sample = sid, n_spot = length(common),
    DRI_genes = length(dri_avail), TDS_genes = length(tds_avail),
    rho_raw = round(ct_raw$estimate, 3), p_raw = ct_raw$p.value,
    rho_partial = round(ct_par$estimate, 3), p_partial = ct_par$p.value,
    moran_I = round(moran$observed, 3), moran_p = moran$p.value))

  saveRDS(list(dri = dri[common], tds = tds[common],
               dri_resid = res_dri, tds_resid = res_tds,
               pos = pos_sub, nCount = nCount[common],
               nFeature = nFeature[common], sample = sid),
          file.path(out_dir, paste0(sid, "_spatial_scores.rds")))
}

write.csv(results, file.path(out_dir, "TableS8_spatial_summary.csv"),
          row.names = FALSE)
cat("\nDone. TableS8_spatial_summary.csv saved.\n")

