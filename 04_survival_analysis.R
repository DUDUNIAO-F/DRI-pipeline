# 04_survival_analysis.R
# Purpose: PFI survival, Cox, ROC
# Input  : output/DRI_score.rds, data/TCGA-THCA_FULL_merged_with_CDR.tsv
# Output : output/survival_results.rds, output/TableS7_survival.csv

library(survival); library(timeROC)
dri <- readRDS('output/DRI_score.rds')
clin <- read.delim('data/TCGA-THCA_FULL_merged_with_CDR.tsv', check.names = FALSE)
id_col <- grep('sample|barcode|submitter', colnames(clin), ignore.case = TRUE, value = TRUE)[1]
df <- data.frame(id12 = substr(gsub('\\.', '-', clin[[id_col]]), 1, 12),
                 PFI = as.numeric(clin$cdr_PFI),
                 PFI_time = as.numeric(clin$cdr_PFI.time),
                 stage = clin$cdr_ajcc_pathologic_tumor_stage)
df <- df[!is.na(df$PFI) & !is.na(df$PFI_time) & df$PFI_time > 0, ]
df$stage_group <- ifelse(grepl('Stage I$|Stage II$', df$stage), 'Early',
                  ifelse(grepl('Stage III|Stage IV', df$stage), 'Advanced', NA))
df <- merge(data.frame(id12 = substr(names(dri), 1, 12), DRI = as.numeric(dri)),
            df, by = 'id12')
df <- df[!is.na(df$DRI) & !is.na(df$stage_group), ]
df$DRI_group <- factor(ifelse(df$DRI > median(df$DRI), 'High DRI', 'Low DRI'),
                       levels = c('Low DRI', 'High DRI'))
df$time_y <- df$PFI_time / 365.25
cox_uni <- coxph(Surv(time_y, PFI) ~ DRI_group, data = df)
cox_adj <- coxph(Surv(time_y, PFI) ~ DRI_group + stage_group, data = df)
print(summary(cox_uni)$conf.int)
print(summary(cox_adj)$conf.int)
roc <- timeROC(T = df$time_y, delta = df$PFI, marker = df$DRI,
               cause = 1, times = c(3, 5), iid = TRUE)
cat('3y AUC:', round(roc$AUC[1], 3), ' 5y AUC:', round(roc$AUC[2], 3), '\n')
saveRDS(list(cox_uni = cox_uni, cox_adj = cox_adj, roc = roc),
        'output/survival_results.rds')
