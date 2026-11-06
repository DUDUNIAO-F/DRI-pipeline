# Minimal test script: Simulates 3 samples and 100 genes to run the pipeline in 5 seconds
source("DRI_core.R")

cat(">>> Generating toy data...\n")
set.seed(123)

# Simulate expression matrix (3 samples, 100 genes) including NSUN2 and TDS genes
genes <- c("NSUN2", "TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1", paste0("Gene", 1:94))
expr_toy <- matrix(rnorm(300), nrow = 100, ncol = 3, dimnames = list(genes, paste0("Sample", 1:3)))

# Introduce artificial correlations to ensure the DRI calculation succeeds
expr_toy["Gene1", ] <- expr_toy["NSUN2", ] * 0.8 + rnorm(3, sd = 0.2)
expr_toy["Gene2", ] <- -expr_toy["TG", ] * 0.8 + rnorm(3, sd = 0.2)

cat(">>> Running DRI framework...\n")
res <- calculate_DRI(
  expr_mat = expr_toy, 
  tds_genes = c("TG", "TSHR", "PAX8", "SLC5A5", "NKX2-1")
)

cat(">>> Success! DRI scores calculated:\n")
print(res$DRI_score)
cat(">>> Number of DRI genes identified:", res$N_genes, "\n")

# Save results
output_file <- "DRI_example_results.csv"
write.csv(data.frame(Sample = names(res$DRI_score), DRI = as.numeric(res$DRI_score)), 
          output_file, row.names = FALSE)
cat(">>> Output saved to", output_file, "\n")