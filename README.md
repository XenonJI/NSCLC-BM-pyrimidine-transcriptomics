# NSCLC brain metastasis pyrimidine-transcription analysis

This repository contains the R scripts supporting the manuscript:

**An apparent pyrimidine transcriptional signature in NSCLC brain metastases is proliferation-coupled and sensitive to cohort structure**

## Data sources

The analysis uses publicly available datasets from NCBI GEO:

- **GSE131907** — exploratory single-cell dataset
- **GSE223499** — primary single-nucleus RNA-seq dataset
- **GSE223501** — same-study Slide-seq spatial dataset

Raw data are not redistributed in this repository.

## Repository structure

```text
01_exploratory/           Exploratory GSE131907 analysis
02_primary_GSE223499/     Primary single-nucleus analysis
03_spatial_GSE223501/     Spatial analysis
04_integrated/            CR2/common-support and specificity analyses
05_reporting/             Final figure and table generation
SCRIPT_MANIFEST.csv       Script descriptions and order
required_packages.txt     Main R package requirements
```

## Local project root

The analysis scripts were developed with the project root:

```r
ROOT <- "D:/LUAD_LM_PYRIMIDINE"
```

Before running the scripts, replace `ROOT` with the local project directory or reproduce the same folder structure.

## Suggested data layout

```text
project_root/
├── data/
│   ├── raw/
│   │   ├── GSE131907/
│   │   ├── GSE223499/
│   │   └── GSE223501/
│   ├── processed/
│   └── reference/
│       └── GSE223499/
│           └── Supplementary_Table_1.xlsx
├── results/
└── scripts/
```

## R packages

Main packages used across the selected scripts include:

`data.table`, `ggplot2`, `readxl`, `MASS`, `clubSandwich`, `patchwork`, `svglite`, and `scales`.

Some scripts additionally use `fwildclusterboot`, `dqrng`, `generics`, or `Seurat`. See `required_packages.txt`.

## Analysis order

### Exploratory analysis

1. `01_exploratory/01_GSE131907_exploratory_analysis.R`

### Primary GSE223499 analysis

1. `02_primary_GSE223499/01_GSE223499_extract_fixed_module.R`
2. `02_primary_GSE223499/02_KRAS17_MTX_reconstruction.R`
3. rerun `02_primary_GSE223499/01_GSE223499_extract_fixed_module.R`
4. `02_primary_GSE223499/03_GSE223499_fixed_robustness.R`
5. `02_primary_GSE223499/04_build_primary_analysis_table.R`
6. `02_primary_GSE223499/05_match_official_metadata.R`
7. `02_primary_GSE223499/06_primary_patient_aware_analysis.R`
8. `02_primary_GSE223499/07_cell_state_sensitivity.R`
9. `02_primary_GSE223499/08_count_pseudobulk_NB.R`
10. `02_primary_GSE223499/09_cohort_provenance_covariates.R`
11. `02_primary_GSE223499/10_covariate_sensitivity.R`

### Spatial GSE223501 analysis

1. `03_spatial_GSE223501/01_spatial_file_audit.R`
2. `03_spatial_GSE223501/02_spatial_matrix_gene_audit.R`
3. `03_spatial_GSE223501/03_extract_spatial_targets.R`
4. `03_spatial_GSE223501/04_compute_spatial_library_depth.R`
5. `03_spatial_GSE223501/05_coordinate_grid_audit.R`
6. `03_spatial_GSE223501/06_spatial_coupling_analysis.R`
7. `03_spatial_GSE223501/07_spatial_scale_sensitivity.R`
8. `03_spatial_GSE223501/08_spatial_duplicate_audit.R`
9. `03_spatial_GSE223501/09_spatial_robustness.R`

### Integrated analyses

1. `04_integrated/01_CR2_common_support_attenuation.R`
2. `04_integrated/02_target_panel_specificity.R`

### Final reporting

1. `05_reporting/01_final_figures_tables.R`

## Notes

- Patient/specimen-level inference is used for the primary comparisons.
- GSE223501 is used as same-study spatial context rather than an independent validation cohort.
- The five-gene panel is treated as discovery-derived, not prospectively preregistered.
- Source datasets should be obtained directly from GEO.

## Authors

Xiang Ji, Zicheng Peng, Peiyu Zhao, Zhuoyan Jiang, Jiqing Chen, Tongyan Liu, and Yin Rong.
