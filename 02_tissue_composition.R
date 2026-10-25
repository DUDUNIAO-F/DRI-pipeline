# 02_tissue_composition.R
# Purpose: ESTIMATE, partial correlation, purity-stratified analysis
# Input  : output/DRI_score.rds, output/TDS_score.rds, output/datExpr_TCGA.rds
# Output : output/TableS5_tissue_composition.csv

library(tidyestimate)
datExpr <- readRDS('output/datExpr_TCGA.rds')
dri <- readRDS('output/DRI_score.rds')
tds <- readRDS('output/TDS_score.rds')
est <- as.data.frame(t(datExpr)) |>
  filter_common_genes(id = 'hgnc_symbol', tell_missing = FALSE, find_alias = TRUE) |>
  estimate_score(is_affymetrix = FALSE)
est$purity <- cos(0.6049872018 + 0.0001467884 * est$estimate)
common <- intersect(est$sample, names(dri))
df <- data.frame(DRI = dri[common], TDS = tds[common],
                 Purity = est$purity[match(common, est$sample)],
                 Immune = est$immune[match(common, est$sample)])
res1 <- residuals(lm(DRI ~ Purity + Immune, data = df))
res2 <- residuals(lm(TDS ~ Purity + Immune, data = df))
cat('Partial rho:', round(cor(res1, res2, method = 'spearman'), 3), '\n')
write.csv(df, 'output/TableS5_tissue_composition.csv', row.names = FALSE)
