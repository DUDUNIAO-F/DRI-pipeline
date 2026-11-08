# 10_downsampling.R
# Purpose: Gene-coverage-matched downsampling in bulk TCGA (FPKM)
# Input  : data/HiSeqV2.tsv (UCSC Xena TCGA-THCA FPKM), data/DRI_genes.txt
# Output : output/Table_S13_downsampling.csv, output/Fig_S7_downsampling.pdf
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

# ---- Portable project root ----
project_root <- Sys.getenv("DRI_PROJECT_ROOT", unset = normalizePath("."))
data_dir <- Sys.getenv("DRI_DATA_DIR", unset = file.path(project_root, "data"))
out_dir  <- Sys.getenv("DRI_OUT_DIR",  unset = file.path(project_root, "output"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==== User-configurable paths ====
tcga_file <- Sys.getenv("DRI_TCGA_FPKM", unset = file.path(data_dir, "HiSeqV2.tsv"))
dri_file  <- Sys.getenv("DRI_GENE_LIST", unset = file.path(data_dir, "DRI_genes.txt"))
stopifnot(file.exists(tcga_file))
stopifnot(file.exists(dri_file))

expr_raw <- read.delim(tcga_file, header = TRUE, row.names = 1, check.names = FALSE)
datExpr <- t(expr_raw)
datExpr <- datExpr[substr(rownames(datExpr), 14, 15) == "01", , drop = FALSE]
datExpr <- datExpr[, !duplicated(colnames(datExpr))]
keep <- apply(datExpr, 2, median) > 1
datExpr <- datExpr[, keep]

tds_genes <- c("TG","TSHR","PAX8","SLC5A5","NKX2-1")
tds_avail <- intersect(tds_genes, colnames(datExpr))
if (length(tds_avail) < 3) stop("Too few TDS genes available: ", paste(tds_avail, collapse = ", "))
tds <- rowMeans(datExpr[, tds_avail, drop = FALSE], na.rm = TRUE)

dri_genes <- readLines(dri_file); dri_genes <- dri_genes[dri_genes != ""]
dri_avail <- intersect(dri_genes, colnames(datExpr))
cat("DRI genes available:", length(dri_avail), "/", length(dri_genes), "\n")
if (length(dri_avail) < 100) stop("Too few DRI genes available: ", length(dri_avail))

n_single <- min(2964, length(dri_avail))
n_iter <- 1000
set.seed(42); rho_down <- numeric(n_iter)
for (i in 1:n_iter) {
  g <- sample(dri_avail, n_single, replace = FALSE)
  rho_down[i] <- cor(rowMeans(scale(datExpr[, g, drop = FALSE]), na.rm = TRUE),
                     tds, method = "spearman", use = "complete.obs")
}
rho_full <- cor(rowMeans(scale(datExpr[, dri_avail, drop = FALSE]), na.rm = TRUE),
                tds, method = "spearman", use = "complete.obs")
cat("Full DRI rho:", round(rho_full, 3),
    "| Downsampled median:", round(median(rho_down), 3),
    "| 95% CI:", round(quantile(rho_down, c(.025, .975)), 3), "\n")

fwrite(data.frame(Iteration = 1:n_iter, Rho = rho_down,
                  N_genes = n_single,
                  Method = paste0("Downsampled DRI (", format(n_single, big.mark = ","), " genes)")),
       file.path(out_dir, "Table_S13_downsampling.csv"))
p <- ggplot(data.frame(Rho = rho_down), aes(x = Rho)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = median(rho_down), color = "red", linetype = "dashed") +
  geom_vline(xintercept = 0.101, color = "black", linetype = "dotted") +
  labs(x = "Spearman rho (Downsampled DRI vs TDS)", y = "Frequency",
       title = "Gene coverage-matched downsampling (TCGA-THCA FPKM)",
       subtitle = sprintf("%s of %s DRI genes; n = %d iterations",
                          format(n_single, big.mark = ","),
                          format(length(dri_avail), big.mark = ","), n_iter)) +
  theme_bw()
ggsave(file.path(out_dir, "Fig_S7_downsampling.pdf"), p, width = 7, height = 5)
cat("Done: Table S13 and Fig S7 saved to", out_dir, "\n")
