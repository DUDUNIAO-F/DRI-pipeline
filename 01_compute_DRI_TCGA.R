# 01_compute_DRI_TCGA.R
# Purpose: Construct DRI in TCGA-THCA and compute DRI-TDS correlation
# Input  : data/TCGA-THCA.htseq_fpkm.tsv.gz
# Output : output/DRI_score.rds, output/TDS_score.rds, output/TableS1_top50_core_genes.csv

expr_raw <- read.delim(gzfile('data/TCGA-THCA.htseq_fpkm.tsv.gz'),
                       header = TRUE, row.names = 1, check.names = FALSE)
datExpr <- t(expr_raw)
datExpr <- datExpr[substr(rownames(datExpr), 14, 15) == '01', , drop = FALSE]
datExpr <- datExpr[, !duplicated(colnames(datExpr))]
keep <- apply(datExpr, 2, median) > 1
expr <- datExpr[, keep]
tds <- rowMeans(expr[, c('TG','TSHR','PAX8','SLC5A5','NKX2-1')], na.rm = TRUE)
nsun2_cor <- cor(expr, expr[, 'NSUN2'], method = 'spearman')[, 1]
tds_cor <- cor(expr, tds, method = 'spearman')[, 1]
core <- names(which(nsun2_cor > 0.15 & tds_cor < -0.10))
supp <- setdiff(names(which(tds_cor < -0.15)), core)
dri_genes <- union(core, supp)
dri <- rowMeans(scale(expr[, dri_genes]), na.rm = TRUE)
saveRDS(dri, 'output/DRI_score.rds')
saveRDS(tds, 'output/TDS_score.rds')
saveRDS(datExpr, 'output/datExpr_TCGA.rds')
cat('DRI-TDS rho:', round(cor(dri, tds, method = 'spearman'), 3), '\n')
