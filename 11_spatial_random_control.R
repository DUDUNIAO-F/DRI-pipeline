# 11_spatial_random_control.R
# Purpose: Spatial DRI vs expression-matched random gene sets (GSE250521)
# Input  : GSE250521_RAW/ Visium files, DRI_genes.txt
# Output : output/Table_S14_spatial_null.csv, output/Fig_S9_spatial_null.pdf
suppressPackageStartupMessages({
  library(Seurat); library(ape); library(Matrix); library(ggplot2); library(data.table)
})

# ==== User-configurable paths ====
data_dir <- Sys.getenv("DRI_SPATIAL_DIR",
  "/Users/fengweiming/Desktop/tsRNA生信文章/data/GSE250521_RAW")
dri_file <- Sys.getenv("DRI_GENE_LIST",
  "/Users/fengweiming/Desktop/tsRNA生信文章/data/未命名文件夹/DRI_genes.txt")
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dri_genes <- readLines(dri_file); dri_genes <- dri_genes[dri_genes != ""]
tds_genes <- c("TG","TSHR","PAX8","SLC5A5","NKX2-1")
pos_files  <- list.files(data_dir, pattern = "_visium_tissue_positions_list.csv.gz$")
sample_ids <- sub("_visium_tissue_positions_list.csv.gz$", "", pos_files)
cat("Samples with positions:", length(sample_ids), "\n")

results <- data.frame()
for (sid in sample_ids) {
  bc <- readLines(gzfile(file.path(data_dir, paste0(sid, "_visium_barcodes.tsv.gz"))))
  counts <- ReadMtx(
    mtx = file.path(data_dir, paste0(sid, "_visium_matrix.mtx.gz")),
    cells = file.path(data_dir, paste0(sid, "_visium_barcodes.tsv.gz")),
    features = file.path(data_dir, paste0(sid, "_visium_features.tsv.gz")),
    feature.column = 2, cell.column = 1)
  colnames(counts) <- bc
  pos <- read.csv(gzfile(file.path(data_dir,
                paste0(sid, "_visium_tissue_positions_list.csv.gz"))), header = FALSE)
  colnames(pos) <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
  pos <- pos[pos$barcode %in% bc, ]
  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj, verbose = FALSE)
  nCount <- colSums(counts); nFeature <- colSums(counts > 0)
  dri_avail <- intersect(dri_genes, rownames(obj))
  tds_avail <- intersect(tds_genes, rownames(obj))
  if (length(dri_avail) < 500 || length(tds_avail) < 3) next
  set.seed(123)
  obj <- AddModuleScore(obj, features = list(dri = dri_avail),
                        name = "DRI_score", nbin = 24, ctrl = 100)
  obj <- AddModuleScore(obj, features = list(tds = tds_avail),
                        name = "TDS_score", nbin = 24, ctrl = 100)
  dri <- obj$DRI_score1; tds <- obj$TDS_score1
  names(dri) <- colnames(obj); names(tds) <- colnames(obj)
  common <- intersect(names(dri), pos$barcode)
  common <- common[is.finite(dri[common]) & is.finite(tds[common])]
  if (length(common) < 50) next
  d <- data.frame(dri = dri[common], tds = tds[common],
                  nCount = nCount[common], nFeature = nFeature[common])
  res_dri <- residuals(lm(dri ~ nCount + nFeature, data = d))
  res_tds <- residuals(lm(tds ~ nCount + nFeature, data = d))
  rho_obs <- cor(res_dri, res_tds, method = "spearman")

  n_random <- 50; rho_null <- numeric(n_random); set.seed(123)
  pool <- setdiff(rownames(obj), tds_avail)
  for (j in 1:n_random) {
    rg <- sample(pool, length(dri_avail), replace = FALSE)
    obj_tmp <- AddModuleScore(obj, features = list(r = rg),
                              name = "Rand_score", nbin = 24, ctrl = 100)
    rs <- obj_tmp$Rand_score1; names(rs) <- colnames(obj_tmp)
    d2 <- data.frame(r = rs[common], nCount = nCount[common], nFeature = nFeature[common])
    res_r <- residuals(lm(r ~ nCount + nFeature, data = d2))
    rho_null[j] <- cor(res_r, res_tds, method = "spearman")
  }
  p_emp <- mean(abs(rho_null) >= abs(rho_obs))
  ci <- quantile(rho_null, c(.025, .975), na.rm = TRUE)
  cat(sid, ": rho =", round(rho_obs, 3), "| null CI =",
      round(ci[1], 3), "-", round(ci[2], 3), "| p =", round(p_emp, 3), "\n")
  results <- rbind(results, data.frame(
    Sample = sid, n_spot = length(common),
    rho_partial = round(rho_obs, 3),
    rho_null_median = round(median(rho_null), 3),
    rho_null_lower = round(ci[1], 3), rho_null_upper = round(ci[2], 3),
    p_empirical = round(p_emp, 3)))
}
fwrite(results, file.path(out_dir, "Table_S14_spatial_null.csv"))
print(results)

if (nrow(results) > 0) {
  p <- ggplot(results, aes(x = Sample)) +
    geom_point(aes(y = rho_partial, color = "Observed DRI"), size = 3) +
    geom_errorbar(aes(ymin = rho_null_lower, ymax = rho_null_upper,
                      color = "Random null 95% CI"), width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "", y = "Partial Spearman rho (DRI vs TDS)",
         title = "Spatial DRI vs expression-matched random gene sets", color = "") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "Fig_S9_spatial_null.pdf"), p, width = 9, height = 5)
}
cat("Done: Table S14 and Fig S9 saved to", out_dir, "\n")
