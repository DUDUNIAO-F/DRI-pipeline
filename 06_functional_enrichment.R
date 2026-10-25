# 06_functional_enrichment.R
# Purpose: GO BP and KEGG enrichment
# Input  : output/datExpr_TCGA.rds, output/TDS_score.rds
# Output : output/TableS2_GO_BP_DRI.csv, output/TableS3_KEGG_DRI.csv

library(clusterProfiler); library(org.Hs.eg.db)
expr <- readRDS('output/datExpr_TCGA.rds')
tds <- readRDS('output/TDS_score.rds')
expr <- expr[, apply(expr, 2, median) > 1]
nsun2_cor <- cor(expr, expr[, 'NSUN2'], method = 'spearman')[, 1]
tds_cor <- cor(expr, tds, method = 'spearman')[, 1]
core <- names(which(nsun2_cor > 0.15 & tds_cor < -0.10))
supp <- setdiff(names(which(tds_cor < -0.15)), core)
dri_genes <- union(core, supp)
ego <- enrichGO(dri_genes, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                ont = 'BP', pAdjustMethod = 'BH', qvalueCutoff = 0.05, readable = TRUE)
write.csv(ego@result, 'output/TableS2_GO_BP_DRI.csv', row.names = FALSE)
gene_entrez <- bitr(dri_genes, 'SYMBOL', 'ENTREZID', org.Hs.eg.db)
ekegg <- enrichKEGG(gene_entrez$ENTREZID, organism = 'hsa',
                    pAdjustMethod = 'BH', qvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, org.Hs.eg.db, keyType = 'ENTREZID')
write.csv(ekegg@result, 'output/TableS3_KEGG_DRI.csv', row.names = FALSE)
