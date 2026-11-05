# 08_framework_necessity.R
# Purpose: Compare DRI against single-gene proxies across cohorts
# Input  : GSE65144/GSE33630/GSE29265 series matrices, GPL570.soft.gz, DRI_genes.txt
# Output : output/Table_S9_framework_necessity.csv, output/Fig_1F_framework_necessity.pdf
suppressPackageStartupMessages({
  library(GEOquery); library(ggplot2); library(data.table)
})

# ==== User-configurable paths ====
data_dir <- Sys.getenv("DRI_DATA_DIR",
  "/Users/fengweiming/Desktop/tsRNA生信文章/data/未命名文件夹")
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gse_files <- c("GSE65144" = "GSE65144_series_matrix.txt.gz",
               "GSE33630" = "GSE33630_series_matrix.txt.gz",
               "GSE29265" = "GSE29265_series_matrix.txt.gz")

gpl <- getGEO(filename = file.path(data_dir, "GPL570.soft.gz"))
gpl_tab <- Table(gpl)
sym_col <- grep("symbol", colnames(gpl_tab), ignore.case = TRUE, value = TRUE)[1]
probe2gene <- gpl_tab[, c("ID", sym_col)]
colnames(probe2gene) <- c("ID", "SYMBOL")
probe2gene <- probe2gene[probe2gene$SYMBOL != "" & !is.na(probe2gene$SYMBOL), ]
probe2gene$SYMBOL <- sapply(strsplit(probe2gene$SYMBOL, " /// "), `[`, 1)

dri_genes <- readLines(file.path(data_dir, "DRI_genes.txt"))
dri_genes <- dri_genes[dri_genes != ""]

process_eset <- function(eset, probe2gene) {
  expr <- exprs(eset)
  qx <- as.numeric(quantile(expr, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  if ((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
    expr[expr <= 0] <- NaN; expr <- log2(expr + 1)
  }
  df <- as.data.frame(expr); df$ID <- rownames(df)
  merged <- merge(probe2gene, df, by = "ID")
  agg <- aggregate(merged[, -c(1, 2)], by = list(SYMBOL = merged$SYMBOL),
                   FUN = mean, na.rm = TRUE)
  rownames(agg) <- agg$SYMBOL; agg$SYMBOL <- NULL
  as.matrix(agg)
}

calc_tds <- function(m) {
  genes <- c("TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1")
  colMeans(m[intersect(genes, rownames(m)), , drop = FALSE], na.rm = TRUE)
}

calc_rho_ci <- function(x, y, R = 1000) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  if (length(x) < 10) return(data.frame(rho = NA, lower = NA, upper = NA, n = length(x)))
  rho <- cor(x, y, method = "spearman")
  br <- replicate(R, { i <- sample.int(length(x), length(x), TRUE)
                       cor(x[i], y[i], method = "spearman") })
  ci <- quantile(br, c(.025, .975), na.rm = TRUE)
  data.frame(rho = rho, lower = ci[1], upper = ci[2], n = length(x))
}

res_list <- list()
for (cohort in names(gse_files)) {
  message("Processing ", cohort)
  eset <- getGEO(filename = file.path(data_dir, gse_files[cohort]))
  expr <- process_eset(eset, probe2gene)
  tds  <- calc_tds(expr)
  z    <- t(scale(t(expr)))
  s <- list(NSUN2 = z["NSUN2", ], METTL1 = z["METTL1", ], WDR4 = z["WDR4", ],
            ThreeGeneMean = colMeans(z[c("NSUN2","METTL1","WDR4"), ], na.rm = TRUE))
  avail <- intersect(dri_genes, rownames(z))
  s[["DRI"]] <- colMeans(z[avail, , drop = FALSE], na.rm = TRUE)
  message("  DRI genes used: ", length(avail), "/", length(dri_genes))
  for (nm in names(s)) {
    r <- calc_rho_ci(s[[nm]], tds)
    res_list[[length(res_list) + 1]] <- data.frame(Cohort = cohort, Proxy = nm, r)
  }
}
res_df <- do.call(rbind, res_list)
fwrite(res_df, file.path(out_dir, "Table_S9_framework_necessity.csv"))
print(res_df)

lv <- c("DRI","ThreeGeneMean","NSUN2","METTL1","WDR4")
res_df$Proxy <- factor(res_df$Proxy, levels = intersect(lv, unique(res_df$Proxy)))
p <- ggplot(res_df, aes(x = rho, y = Proxy, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.2, position = position_dodge(width = 0.6)) +
  labs(x = "Spearman rho (Proxy vs TDS)", y = "",
       title = "Framework necessity: DRI vs single-gene proxies") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "Fig_1F_framework_necessity.pdf"), p, width = 8, height = 5)
cat("Done: Table_S9 and Fig_1F saved to", out_dir, "\n")
