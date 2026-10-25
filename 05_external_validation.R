# 05_external_validation.R
# Purpose: DRI-TDS in GSE65144, GSE33630, GSE29265
# Input  : data/GSE65144_series_matrix.txt.gz, etc.
# Output : output/TableS4_external_validation.csv

library(GEOquery); library(hgu133plus2.db); library(dplyr)
dri_genes <- readRDS('output/DRI_score.rds')
tds_genes <- c('TG','TSHR','PAX8','SLC5A5','NKX2-1')
validate_cohort <- function(f, id) {
  if (!file.exists(f)) return(data.frame(Cohort = id, status = 'missing'))
  gset <- getGEO(filename = f, GSEMatrix = TRUE)
  if (is.list(gset) && !is(gset, 'ExpressionSet')) gset <- gset[[1]]
  expr <- exprs(gset)
  sym <- mapIds(hgu133plus2.db, rownames(expr), 'SYMBOL', 'PROBEID', multiVals = 'first')
  d <- as.data.frame(expr); d$Gene <- sym; d <- d[!is.na(d$Gene), ]
  g <- d |> group_by(Gene) |> summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')
  m <- as.matrix(g[, -1]); rownames(m) <- g$Gene
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)
  da <- intersect(names(dri_genes), rownames(m)); ta <- intersect(tds_genes, rownames(m))
  dri <- colMeans(t(scale(t(m[da, , drop = FALSE]))), na.rm = TRUE)
  tds <- colMeans(m[ta, , drop = FALSE], na.rm = TRUE)
  ct <- cor.test(dri, tds, method = 'spearman')
  data.frame(Cohort = id, n = ncol(m), DRI_genes = length(da),
             TDS_genes = length(ta), rho = round(ct$estimate, 3),
             p = ct$p.value, status = 'ok')
}
res <- bind_rows(
  validate_cohort('data/GSE65144_series_matrix.txt.gz', 'GSE65144'),
  validate_cohort('data/GSE33630_series_matrix.txt.gz', 'GSE33630'),
  validate_cohort('data/GSE29265_series_matrix.txt.gz', 'GSE29265'))
print(res)
write.csv(res, 'output/TableS4_external_validation.csv', row.names = FALSE)
