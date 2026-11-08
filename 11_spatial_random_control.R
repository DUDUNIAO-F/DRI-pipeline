# 11_spatial_random_control.R
# Purpose: Spatial DRI vs permutation null (GSE250521)
# Input  : data/GSE250521_RAW/ Visium files, data/DRI_genes.txt
# Output : output/Table_S14_spatial_null.csv, output/Fig_S8_spatial_null.pdf
suppressPackageStartupMessages({
  library(Seurat); library(ape); library(Matrix); library(ggplot2); library(data.table)
})

# ---- Portable project root ----
project_root <- Sys.getenv("DRI_PROJECT_ROOT", unset = normalizePath("."))
data_root <- Sys.getenv("DRI_DATA_DIR", unset = file.path(project_root, "data"))
out_dir   <- Sys.getenv("DRI_OUT_DIR",  unset = file.path(project_root, "output"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==== User-configurable paths ====
data_dir <- Sys.getenv("DRI_SPATIAL_DIR", unset = file.path(data_root, "GSE250521_RAW"))
dri_file <- Sys.getenv("DRI_GENE_LIST", unset = file.path(data_root, "DRI_genes.txt"))
stopifnot(dir.exists(data_dir))
stopifnot(file.exists(dri_file))

dri_genes <- readLines(dri_file); dri_genes <- dri_genes[dri_genes != ""]
tds_genes <- c("TG","TSHR","PAX8","SLC5A5","NKX2-1")
pos_files  <- list.files(data_dir, pattern = "_visium_tissue_positions_list.csv.gz$")
if (length(pos_files) == 0) stop("No Visium position files found in: ", data_dir)
sample_ids <- sub("_visium_tissue_positions_list.csv.gz$", "", pos_files)
cat("Samples with positions:", length(sample_ids), "\n")

results <- data.frame(
  Sample = character(), n_spot = integer(),
  rho_partial = numeric(), rho_null_median = numeric(),
  rho_null_lower = numeric(), rho_null_upper = numeric(),
  p_empirical = numeric(), p_empirical_str = character(),
  stringsAsFactors = FALSE
)

n_perm <- 1000
set.seed(123)

for (sid in sample_ids) {
  mtx_file  <- file.path(data_dir, paste0(sid, "_visium_matrix.mtx.gz"))
  bc_file   <- file.path(data_dir, paste0(sid, "_visium_barcodes.tsv.gz"))
  feat_file <- file.path(data_dir, paste0(sid, "_visium_features.tsv.gz"))
  pos_file  <- file.path(data_dir, paste0(sid, "_visium_tissue_positions_list.csv.gz"))
  if (!all(file.exists(c(mtx_file, bc_file, feat_file, pos_file)))) {
    warning("Missing files for sample ", sid, "; skipped.")
    next
  }
  counts <- ReadMtx(
    mtx = mtx_file,
    cells = bc_file,
    features = feat_file,
    feature.column = 2, cell.column = 1)
  pos <- read.csv(gzfile(pos_file), header = FALSE)
  colnames(pos) <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
  pos <- pos[pos$in_tissue == 1, ]
  pos <- pos[pos$barcode %in% colnames(counts), ]
  if (nrow(pos) < 50) next

  obj <- CreateSeuratObject(counts = counts)
  obj <- NormalizeData(obj, verbose = FALSE)
  nCount <- colSums(counts); nFeature <- colSums(counts > 0)
  dri_avail <- intersect(dri_genes, rownames(obj))
  tds_avail <- intersect(tds_genes, rownames(obj))
  if (length(dri_avail) < 500 || length(tds_avail) < 3) next

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
  rho_obs <- cor(res_dri, res_tds, method = "spearman", use = "complete.obs")

  # Permutation null: shuffle depth-corrected DRI residuals across spots
  rho_null <- numeric(n_perm)
  for (j in 1:n_perm) {
    res_r <- sample(res_dri)
    rho_null[j] <- cor(res_r, res_tds, method = "spearman", use = "complete.obs")
  }
  p_emp <- mean(abs(rho_null) >= abs(rho_obs), na.rm = TRUE)
  ci <- quantile(rho_null, c(.025, .975), na.rm = TRUE)
  p_str <- if (p_emp == 0) paste0("<", 1/(n_perm+1)) else sprintf("%.4f", p_emp)
  cat(sid, ": rho =", round(rho_obs, 3), "| perm null CI =",
      round(ci[1], 3), "-", round(ci[2], 3), "| p =", p_str, "\n")
  results <- rbind(results, data.frame(
    Sample = sid, n_spot = length(common),
    rho_partial = round(rho_obs, 3),
    rho_null_median = round(median(rho_null, na.rm = TRUE), 3),
    rho_null_lower = round(ci[1], 3), rho_null_upper = round(ci[2], 3),
    p_empirical = round(p_emp, 4), p_empirical_str = p_str,
    stringsAsFactors = FALSE))
}

if (nrow(results) > 0) {
  fwrite(results, file.path(out_dir, "Table_S14_spatial_null.csv"))
  print(results)

  p <- ggplot(results, aes(x = Sample)) +
    geom_point(aes(y = rho_partial, color = "Observed DRI"), size = 3) +
    geom_errorbar(aes(ymin = rho_null_lower, ymax = rho_null_upper,
                      color = "Permutation null 95% CI"), width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "", y = "Partial Spearman rho (DRI vs TDS)",
         title = "Spatial DRI vs permutation null",
         subtitle = paste0(n_perm, " spatial permutations of depth-corrected DRI residuals"),
         color = "") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "Fig_S8_spatial_null.pdf"), p, width = 9, height = 5)
  cat("Done: Table S14 and Fig S8 saved to", out_dir, "\n")
} else {
  warning("No samples passed filters; no output written.")
}
