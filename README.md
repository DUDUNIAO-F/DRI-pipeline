# DRI-pipeline

Code accompanying the manuscript:

**A tRNA Regulatory Index (DRI) quantifies tsRNA network activity and positions it as a downstream mediator of MAPK-driven dedifferentiation in thyroid cancer**

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

## Software

R 4.5.2; WGCNA, GSVA, tidyestimate, mediation, survival, survminer, timeROC, rms, clusterProfiler, GEOquery, hgu133plus2.db

## Usage

```bash
bash run_all.sh
```

## License

MIT
