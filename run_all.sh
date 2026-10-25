#!/bin/bash
set -e
mkdir -p output
echo "Step 1"; Rscript 01_compute_DRI_TCGA.R
echo "Step 2"; Rscript 02_tissue_composition.R
echo "Step 3"; Rscript 03_mediation_analysis.R
echo "Step 4"; Rscript 04_survival_analysis.R
echo "Step 5"; Rscript 05_external_validation.R
echo "Step 6"; Rscript 06_functional_enrichment.R
echo "Done"
