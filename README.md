# DRI-pipeline

Code accompanying the manuscript:

**A tRNA Regulatory Index (DRI) quantifies tsRNA network activity and positions it as a downstream mediator of MAPK-driven dedifferentiation in thyroid cancer**
---
## 🚀 Quick Start (5-second test)
Clone this repository and run the following command in R. It will generate a simulated dataset, compute the DRI score, and output the results instantly. No external data is required.

```R
source("run_example.R")

🧩 Core Function

The DRI_core.R file contains the parameterized calculate_DRI() function. It allows users to apply the framework to any cancer type by providing a user-defined differentiation gene list and an anchor gene (default: NSUN2).

🌐 Web Interface (Shiny App)

For non-computational users, a Shiny app is provided in app.R. You can run it locally by executing shiny::runApp("app.R") or use it to calculate DRI by uploading a standard bulk RNA-seq CSV matrix via a user-friendly interface.

## Directory structure

```
DRI-pipeline/
├── data/          # input data (not included)
├── output/        # generated tables and figures
├── 01_compute_DRI_TCGA.R
├── 02_tissue_composition.R
├── 03_mediation_analysis.R
├── 04_survival_analysis.R
├── 05_external_validation.R
├── 06_functional_enrichment.R
├── 07_spatial_analysis.R
└── run_all.sh
```

## Input data

Place the following in `data/`:

| File | Source |
|:---|:---|
| TCGA-THCA.htseq_fpkm.tsv.gz | UCSC Xena |
| THCA_mc3_gene_level.txt.gz | UCSC Xena |
| TCGA-THCA_FULL_merged_with_CDR.tsv | tcga-gdc-downloader |
| GSE65144/GSE33630/GSE29265_series_matrix.txt.gz | GEO |
| GSE250521_* | GEO (spatial transcriptomics, Visium) |

## Software

R 4.5.2; WGCNA, GSVA, tidyestimate, mediation, survival, survminer, timeROC, rms, clusterProfiler, GEOquery, hgu133plus2.db, ape

## Usage

```bash
bash run_all.sh
```

## License

MIT
## Analysis scripts

| Script | Purpose | Key outputs |
|---|---|---|
| `01_compute_DRI_TCGA.R` | Construct DRI in TCGA-THCA | DRI score, TDS, DRI gene list |
| `02_tissue_composition.R` | ESTIMATE-based purity/stromal/immune adjustment | Partial correlations |
| `03_mediation_analysis.R` | Formal mediation of BRAF/RAS → DRI → TDS | ACME, ADE, proportion mediated |
| `04_survival_analysis.R` | PFI survival analysis and nomogram | KM curves, Cox models, ROC |
| `05_external_validation.R` | DRI–TDS correlation in three GEO cohorts | Cross-cohort rho |
| `06_functional_enrichment.R` | GO/KEGG enrichment of DRI signature | Enrichment tables |
| `07_spatial_analysis.R` | Visium spatial analysis of DRI | Per-sample partial rho, Moran's I |
| `08_framework_necessity.R` | DRI vs single-gene proxies across cohorts | Table S9, Fig. 1F |
| `09_smallRNA_validation.R` | Independent small RNA-seq validation (9 metastatic vs 9 primary) | Tables S11–S12, Figs S6–S7 |
| `10_downsampling.R` | Gene-coverage-matched downsampling in TCGA FPKM | Table S13, Fig. S8 |
| `11_spatial_random_control.R` | Spatial DRI vs expression-matched random gene sets | Table S14, Fig. S9 |

## Input data

All scripts read from a user-specified data directory. The default path is set via
environment variables (e.g., `DRI_DATA_DIR`, `DRI_TCGA_FPKM`, `DRI_GENE_LIST`,
`DRI_SPATIAL_DIR`) and falls back to a local path. Public datasets used:

- TCGA-THCA FPKM: UCSC Xena `HiSeqV2`
- GSE65144, GSE33630, GSE29265 (Affymetrix U133 Plus 2.0)
- GSE184362 (single-cell RNA-seq)
- GSE250521 (Visium spatial transcriptomics)
- 18-sample small RNA-seq validation cohort (available from the corresponding author upon reasonable request)

No raw data are deposited in this repository. `.gitignore` excludes all data
files (`*.xlsx`, `*.tsv`, `*.mtx`, `*.soft.gz`, etc.).

## Reproducibility

All analyses were performed in R 4.5.2. Each script writes its outputs to
`output/` and prints key statistics to the console. Software versions and
package dependencies are documented in the header of each script.
