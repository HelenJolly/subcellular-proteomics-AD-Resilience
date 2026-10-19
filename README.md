# Subcellular-proteomics-AD-Resilience: Analysis code
Analysis code accompanying our study of subcellular proteomics in ROSMAP samples
This repository contains R Markdown workflow for a study of subcellular proteomics in ROSMAP samples. Scripts are being prepared for manuscript submission.

## Analysis index

| File | Purpose | Manuscript output |
| --- | --- | --- |
| `protein_matrix_QC.Rmd` | Protein-level filtering, log2 transformation, MAI imputation, SummarizedExperiment (SE) storage and run-level PCA | Processed protein matrix, SE object, and Supplementary Figure S1 |
| `Figure2_subcell_profiles.Rmd` | Protein fractionation profile clustering, GO:CC enrichment, curated marker profile visualisation, t-SNE and within-group subject correlations | Figure 2A, 2C, 2D and Supplementary Figure S2 |
| `Figure3_Modelling_categorical.Rmd` | Relative predictor importance, disease-by-fraction mixed models, comparison of imputed and unimputed results, protein profiles and GO:CC enrichment | Figure 3A–E |
| `Figure3_Modelling_Resilience.Rmd` | Longitudinal cognition modelling, derivation of resilience scores and resilience-by-fraction mixed models | Figure 3F–H |
| `Figure5_peptide_analysis.Rmd` | Peptide-level filtering, MAPT peptide inspection, mapping, and plotting | Figure 5 |

## Data availability

This repository contains analysis code, not restricted participant-level data.
The code can be inspected without the input data; some elements of execution require authorised access to the restricted metadata:
ROSMAP clinical and neuropathological data access is subject to the applicable data-use agreements, through [Rush Alzheimer's Disease Center](https://www.radc.rush.edu) and [Synapse](https://www.synapse.org).

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

## Disease and resilience modelling workflows

### Disease-associated protein fractionation — Figure 3A–E

`Figure3_Modelling_categorical.Rmd` prepares protein profiles from the imputed and unimputed assays and identifies those eligible for localisation modelling. Relative predictor importance is compared using drop-one-out changes in marginal R2. Disease-associated (CON/RES/AD-DEM) differences in fractionation profiles are identified using mixed models with disease-group x fraction interactions.

Results are compared between imputed and unimputed analyses to identify robust associations. Example protein profiles and GO Cellular Component enrichment analysis is also included.

Required inputs are:

- `MAI_prot_model_new.rds` 
- `demogs.xlsx` 


### Resilience-associated protein fractionation — Figure 3F–H

`Figure3_Modelling_Resilience.Rmd` fits a longitudinal mixed model of cognition using the full available ROSMAP cohort, retaining individuals with at least three lognitudinal cognitive observations and the required covariates.
These resilience scores are linked to the proteomics cohort and used in the localisation mixed model to identify a resilience x fraction interaction effect. Results from imputed and unimputed protein intensities are compared to identify robust associations. The script also includes cognitive trajectory and resilience visualisations.

Additional required inputs are:

- `dataset_longitudinal.xlsx`: restricted longitudinal cognition data.
- `dataset_cross_sectional.xlsx`: restricted pathology and demographic covariates.
- `demogs.xlsx`

The protein preparation and eligibility sections of `Figure3_Modelling_categorical.Rmd` need to be run before the resilience workflow. Intermediate data and model outputs are not included in this repository. 

## MAPT peptides workflow

`Figure5_peptide_analysis.Rmd` imports DIA-NN precursor intensities (These can be downloaded from ProteomeXchange), maps runs to the experimental layout, removes excluded samples and contaminants, and filters by missingness to obtain a clean PR (precursor) matrix. Unmodified precursors are median-aggregated. MAPT peptides with <5% missingness are mapped to the canonical 2N4R sequence and plotted by disease-profile. The accompanying reference cytoskeletal and cytosolic protein-level compartment profiles can be obtained in the figure 2Ci section of `Figure2_subcell_profiles.Rmd`. 

## Software

Analyses use R. Packages referenced include: tidyverse, data.table, readr, SummarizedExperiment, S4Vectors, MAI, dplyr, tidyr, tibble, purrr, stringr, forcats, ggplot2, patchwork, viridis, scales, paletteer, ComplexHeatmap, circlize, clusterProfiler, org.Hs.eg.db, tidytext, and Rtsne. Generating R Markdown also requires rmarkdown and knitr.

## Manuscript and citation

The full manuscript title, author list, citation, and contact details are detailed in the submission. Further scripts will be listed in the analysis index as they are prepared.
