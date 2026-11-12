# 09_smallRNA_validation.R
# Purpose: Independent small RNA-seq validation (9 metastatic vs 9 primary)
# Input  : M_vs_Primary_PTC.xlsx
# Output : output/Table_S11_tsRNA_activity_score.csv, Table_S12_mt_signature_score.csv
#          output/Fig_2A_tsRNA_validation.pdf, Fig_2E_mt_signature.pdf
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(effsize); library(data.table)
})

# ==== User-configurable paths ====
data_dir <- Sys.getenv("DRI_DATA_DIR",
  "/Users/fengweiming/Desktop/tsRNA生信文章/data/未命名文件夹")
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_excel(file.path(data_dir, "M_vs_Primary_PTC.xlsx"), skip = 4)
expr_df  <- df[, 10:27]
expr_mat <- matrix(NA_real_, nrow = nrow(df), ncol = 18)
for (i in 1:18) expr_mat[, i] <- suppressWarnings(as.numeric(as.character(expr_df[[i]])))
rownames(expr_mat) <- make.unique(as.character(df[[2]]))
sample_names <- c("BoM-1","BoM-2","BoM-3","BrM-1","BrM-2","BrM-3",
                  "LuM-1","LuM-2","LuM-3", paste0("Primary_PTC-", 1:9))
colnames(expr_mat) <- sample_names

expr_mat <- expr_mat[rowSums(!is.na(expr_mat)) > 0, ]
for (i in 1:nrow(expr_mat)) {
  if (any(is.na(expr_mat[i, ]))) {
    expr_mat[i, is.na(expr_mat[i, ])] <- median(expr_mat[i, ], na.rm = TRUE)
  }
}
group <- factor(c(rep("Metastasis", 9), rep("Primary", 9)),
                levels = c("Primary", "Metastasis"))

# ---- Global activity score ----
z <- t(scale(t(expr_mat)))
vars <- apply(z, 1, var, na.rm = TRUE)
top_var <- names(sort(vars, decreasing = TRUE))[1:round(length(vars) * 0.5)]
score <- colMeans(z[top_var, , drop = FALSE], na.rm = TRUE)
w  <- wilcox.test(score ~ group); cd <- cliff.delta(score ~ group)
med <- tapply(score, group, median)
cat("== Global activity score ==\nWilcoxon p:", w$p.value,
    "| Cliff's delta:", cd$estimate, "\n")
fwrite(data.frame(Comparison = "M_vs_P", Wilcoxon_p = w$p.value,
                  Cliff_delta = cd$estimate,
                  CI_lower = cd$conf.int[1], CI_upper = cd$conf.int[2],
                  Median_Primary = med["Primary"],
                  Median_Metastasis = med["Metastasis"]),
       file.path(out_dir, "Table_S11_tsRNA_activity_score.csv"))
p1 <- ggplot(data.frame(Sample = sample_names, Group = group, Score = score),
             aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = 0.4, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, aes(color = Group)) + theme_bw() +
  labs(title = "tsRNA activity score (top 50% variable)", y = "Score (mean z-score)")
ggsave(file.path(out_dir, "Fig_2A_tsRNA_validation.pdf"), p1, width = 6, height = 5)

# ---- Focused mt-tsRNA signature ----
sig_ids <- c("mt-tRF3-27-LysTTT","mt-tRF3-21-LysTTT","mt-tRF3-37-LysTTT",
             "mt-tRF3-29-LysTTT_2","mt-tRF3-28-LysTTT","mt-tRF3-29-LysTTT_1",
             "mt-3'tiRNA-44-LysTTT","mt-tRF3-23-LysTTT",
             "mt-3'tiRNA-38-TrpTCA","mt-tRF3-32-TrpTCA")
matched <- character(0)
for (id in sig_ids) {
  if (id %in% rownames(expr_mat)) matched <- c(matched, id)
  else {
    base <- sub("_[12]$", "", id)
    hit <- rownames(expr_mat)[grepl(base, rownames(expr_mat), fixed = TRUE)]
    matched <- c(matched, hit[1])
  }
}
matched <- unique(matched[!is.na(matched)])
cat("Matched mt-tsRNAs:", length(matched), "/ 10\n")
sig_z <- t(scale(t(expr_mat[matched, , drop = FALSE])))
sig_score <- colMeans(sig_z, na.rm = TRUE)
w2 <- wilcox.test(sig_score ~ group); cd2 <- cliff.delta(sig_score ~ group)
med2 <- tapply(sig_score, group, median)
cat("== mt-tsRNA signature ==\nWilcoxon p:", w2$p.value,
    "| Cliff's delta:", cd2$estimate, "\n")
fwrite(data.frame(Signature = "10_mt_tsRNAs", N = length(matched),
                  Wilcoxon_p = w2$p.value, Cliff_delta = cd2$estimate,
                  CI_lower = cd2$conf.int[1], CI_upper = cd2$conf.int[2],
                  Median_Primary = med2["Primary"],
                  Median_Metastasis = med2["Metastasis"]),
       file.path(out_dir, "Table_S12_mt_signature_score.csv"))
p2 <- ggplot(data.frame(Sample = sample_names, Group = group, Score = sig_score),
             aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(width = 0.4, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, aes(color = Group)) + theme_bw() +
  labs(title = "mt-tsRNA signature score (10 pan-metastasis mt-tsRNAs)",
       y = "Score (mean z-score)")
ggsave(file.path(out_dir, "Fig_2E_mt_signature.pdf"), p2, width = 6, height = 5)
cat("Done: Tables S11/S12 and Figs S6/S7 saved to", out_dir, "\n")
