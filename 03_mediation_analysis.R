# 03_mediation_analysis.R
# Purpose: Mediation BRAF/RAS -> DRI -> TDS
# Input  : output/DRI_score.rds, output/TDS_score.rds, data/THCA_mc3_gene_level.txt.gz
# Output : output/mediation_result.rds, output/TableS6_mutation_mediation.csv

library(mediation)
dri <- readRDS('output/DRI_score.rds')
tds <- readRDS('output/TDS_score.rds')
mut <- read.table(gzfile('data/THCA_mc3_gene_level.txt.gz'),
                  header = TRUE, row.names = 1, check.names = FALSE,
                  sep = '\t', fill = TRUE, comment.char = '', quote = '')
drivers <- intersect(c('BRAF','NRAS','HRAS','KRAS'), rownames(mut))
mut_any <- colSums(mut[drivers, , drop = FALSE] > 0) > 0
df <- data.frame(id12 = substr(names(dri), 1, 12),
                 DRI = as.numeric(dri), TDS = as.numeric(tds[names(dri)]))
df$Mutant <- mut_any[match(df$id12, substr(gsub('\\.', '-', names(mut_any)), 1, 12))]
df <- df[!is.na(df$Mutant), ]
df$MutantBinary <- as.numeric(df$Mutant)
set.seed(123)
res <- mediate(lm(DRI ~ MutantBinary, data = df),
               lm(TDS ~ MutantBinary + DRI, data = df),
               treat = 'MutantBinary', mediator = 'DRI',
               boot = 1000, sims = 1000)
print(summary(res))
saveRDS(res, 'output/mediation_result.rds')
