# Subcellular-proteomics-AD-Resilience: Analysis code
Analysis code accompanying our study of subcellular proteomics in ROSMAP samples
This repository contains R Markdown workflow for a study of subcellular proteomics in ROSMAP samples. Scripts are being prepared for manuscript submission.

## Analysis index

| File | Purpose | Manuscript output |
| --- | --- | --- |
| `protein_matrix_QC.Rmd` | Protein-level filtering, log2 transformation, MAI imputation, and run-level PCA | Processed protein matrix and Supplementary Figure S1 |

## Data availability

This repository contains analysis code, not restricted participant-level data.
The code can be inspected without the input data; some elements of execution require authorised access to the restricted metadata:
ROSMAP clinical and neuropathological data access is subject to the applicable
data-use agreements, through [Rush Alzheimer's Disease Center](https://www.radc.rush.edu) and [Synapse](https://www.synapse.org).

## Protein matrix QC workflow

The script imports DIA-NN protein-group intensities, maps runs to the experimental
layout, removes excluded samples and contaminants, filters proteins with at least
70% missingness, removes entries listing multiple accessions or genes, transforms
intensities to log2 scale, and applies MAI imputation. The supplementary PCA uses
imputed intensities and is centred but not scaled to unit variance.

Configure `data_dir` and `output_dir` in the script. The current code reads:

- `data/report.pg_matrix.tsv`: DIA-NN protein-group intensities.
- `data/mapnamesreport.csv`: experimental sample mapping.
- `data/contaminants.csv`: contaminant exclusion list.

The PCA also requires an externally prepared `demogs` object with `Subject`,
`Group`, `age_death`, `pmi`, `study`, and `msex` columns.

The processed SummarizedExperiment contains `log2Quant` and `MAI_imputed`
assays. The manuscript object has 7,271 proteins and 524 runs. The script's
`saveRDS` line is currently commented out; the PCA section reads an existing
`outputs_local/MAI_prot_model_new.rds`. To regenerate that file locally, enable
the save line after configuring paths and obtaining the required inputs.

## Software

Analyses use R. Packages referenced include tidyverse, data.table, readr,
SummarizedExperiment, S4Vectors, MAI, dplyr, tibble, ggplot2, patchwork, viridis,
and scales. Rendering an R Markdown document additionally requires rmarkdown
and knitr.

R and package versions from the original analysis environment should be recorded
in `session_info.txt` before final submission. This repository does not currently include a verified environment lockfile.

## Manuscript and citation

The full manuscript title, author list, citation, and contact details are detailed in the submission. Further scripts will be listed in the analysis index as they are prepared.
