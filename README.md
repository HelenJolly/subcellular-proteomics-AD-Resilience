# Subcellular-proteomics-AD-Resilience: Analysis code
Analysis code accompanying our study of subcellular proteomics in ROSMAP samples
This repository contains R Markdown workflow for a study of subcellular proteomics in ROSMAP samples. Scripts are being prepared for manuscript submission.

## Analysis index

| File | Purpose | Manuscript output |
| --- | --- | --- |
| `protein_matrix_QC.Rmd` | Protein-level filtering, log2 transformation, MAI imputation, SummarizedExperiment (SE) storage and run-level PCA | Processed protein matrix, SE object, and Supplementary Figure S1 |
| `Figure2_subcell_profiles.Rmd` | Protein fractionation profile clustering, GO:CC enrichment, curated marker profile visualisation, t-SNE and within-group subject correlations | Figure 2A, 2C, 2D and Supplementary Figure S2 |

## Data availability

This repository contains analysis code, not restricted participant-level data.
The code can be inspected without the input data; some elements of execution require authorised access to the restricted metadata:
ROSMAP clinical and neuropathological data access is subject to the applicable
data-use agreements, through [Rush Alzheimer's Disease Center](https://www.radc.rush.edu) and [Synapse](https://www.synapse.org).

## Protein matrix QC workflow

The script imports DIA-NN protein-group intensities, maps runs to the experimental layout, removes excluded samples and contaminants, filters proteins with at least 70% missingness, removes entries listing multiple accessions or genes, transforms intensities to log2 scale, and applies MAI imputation.

Configure `data_dir` and `output_dir` in the script. The current code reads:

- `data/report.pg_matrix.tsv`: DIA-NN protein-group intensities.
- `data/mapnamesreport.csv`: experimental sample mapping.
- `data/contaminants.csv`: contaminant exclusion list.

The PCA also requires an externally prepared `demogs` object with `Subject`, `Group`, `age_death`, `pmi`, `study`, and `msex` columns.

The processed SummarizedExperiment contains `log2Quant` and `MAI_imputed` assays. The manuscript object has 7,271 proteins and 524 runs. The script's `saveRDS` line is currently commented out; the PCA section reads an existing `outputs_local/MAI_prot_model_new.rds` which can be generated locally using the raw pg_matrix from PRIDE.

## Subcellular profiles workflow

`Figure2_subcell_profiles.Rmd` uses the MAI-imputed protein matrix to assess subcellular fractionation and reproducibility across individuals. 

Analyses include:
- Hierarchical clustering of profiles and GO:CC enrichment (Figure 2A).
- Curated marker profiles (Figure 2Ci–ii and Supplementary Figure S2).
- t-SNE of scaled mean proportional profiles, coloured by marker compartment (Figure 2Ciii).
- Pairwise subject correlations of log2-imputed protein intensities across fractions, by disease group (Figure 2D).

Configure `data_dir` to point to the directory containing:

- `data/MAI_prot_model_new.rds`: processed SummarizedExperiment from the protein matrix QC workflow.
- `data/merged_markers_filtered.csv`: final curated marker assignments, with `Protein.Group` and `markers` columns.
- `data/sample_category.csv`: subject-to-disease-group mapping, with `Subject` and `Group` columns (CON`,`RES`,`AD-DEM`.

Automatic code execution during knitting is disabled (`eval = FALSE`). Individual chunks can be run locally once the required inputs and packages are available.

## Software

Analyses use R. Packages referenced include: tidyverse, data.table, readr, SummarizedExperiment, S4Vectors, MAI, dplyr, tidyr, tibble, purrr, stringr, forcats, ggplot2, patchwork, viridis, scales, paletteer, ComplexHeatmap, circlize, clusterProfiler, org.Hs.eg.db, tidytext, and Rtsne. Generating R Markdown also requires rmarkdown and knitr.

## Manuscript and citation

The full manuscript title, author list, citation, and contact details are detailed in the submission. Further scripts will be listed in the analysis index as they are prepared.
