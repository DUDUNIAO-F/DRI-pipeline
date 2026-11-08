# make_template.R
# Generate www/DRI_input_template.csv from REAL TCGA-THCA FPKM data.
# Why real data: calculate_DRI() requires (1) NSUN2 as anchor gene,
# (2) genes passing nsun2_cor > 0.15 & tds_cor < -0.10, which random data cannot satisfy.

tcga_file <- Sys.getenv("DRI_TCGA_FPKM", unset = file.path("data", "HiSeqV2.tsv"))
dri_file  <- Sys.getenv("DRI_GENE_LIST", unset = file.path("data", "DRI_genes.txt"))
if (!file.exists(tcga_file)) stop("HiSeqV2.tsv not found. Set DRI_TCGA_FPKM.")
if (!file.exists(dri_file))  stop("DRI_genes.txt not found. Set DRI_GENE_LIST.")

expr_raw <- read.delim(tcga_file, header = TRUE, row.names = 1, check.names = FALSE)
dri_genes <- readLines(dri_file); dri_genes <- dri_genes[dri_genes != ""]
tds_genes <- c("TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1")
anchor    <- "NSUN2"

required <- c(anchor, tds_genes)
missing_req <- setdiff(required, rownames(expr_raw))
if (length(missing_req) > 0) stop("Required genes missing: ", paste(missing_req, collapse = ", "))

dri_avail <- intersect(dri_genes, rownames(expr_raw))
if (length(dri_avail) < 100) stop("Too few DRI genes in expression matrix: ", length(dri_avail))
n_dri <- min(500, length(dri_avail))

genes <- unique(c(required, head(dri_avail, n_dri)))
genes <- intersect(genes, rownames(expr_raw))

set.seed(1)
n_samples <- 20
if (ncol(expr_raw) < n_samples) n_samples <- ncol(expr_raw)
sample_cols <- sample(colnames(expr_raw), n_samples)

mat <- expr_raw[genes, sample_cols, drop = FALSE]
mat <- round(as.matrix(mat), 3)

dir.create("www", showWarnings = FALSE)
out <- file.path("www", "DRI_input_template.csv")
write.csv(mat, out, quote = FALSE)

cat("Template written:", out, "\n")
cat("Genes:", length(genes), "(", n_dri, "DRI +", length(tds_genes), "TDS + NSUN2 )\n")
cat("Samples:", n_samples, "\n")

# Self-test: must pass calculate_DRI before we ship the template
source("DRI_core.R")
res <- calculate_DRI(mat, tds_genes)
cat("Self-test OK: DRI genes selected =", res$N_genes, "\n")
if (res$N_genes == 0) stop("Template self-test failed: no DRI genes passed filter.")
