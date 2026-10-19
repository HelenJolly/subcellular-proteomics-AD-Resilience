# Subcellular-proteomics-AD-Resilience: Analysis code
Analysis code accompanying the study of subcellular proteomics in ROSMAP samples. This repository contains R Markdown scripts used for the relevant analyses and figure generation. Automatic code execution during knitting is disabled (`eval = FALSE`). Individual chunks can be run locally where required inputs and packages are available.

## Index

| File | Purpose | Manuscript output |
| --- | --- | --- |
| `protein_matrix_QC.Rmd` | Protein-level filtering, log2 transformation, MAI imputation, SummarizedExperiment (SE) storage and run-level PCA | Processed protein matrix, SE object, and Supplementary Figure S1 |
| `Figure2_subcell_profiles.Rmd` | Protein fractionation profile clustering, GO:CC enrichment, curated marker profile visualisation, t-SNE and within-group subject correlations | Figure 2A, 2C, 2D and Supplementary Figure S2 |
| `Figure3_Modelling_categorical.Rmd` | Relative predictor importance, disease-by-fraction mixed models, comparison of imputed and unimputed results, protein profiles and GO:CC enrichment | Figure 3A–E |
| `Figure3_Modelling_Resilience.Rmd` | Longitudinal cognition modelling, derivation of resilience scores and resilience-by-fraction mixed models | Figure 3F–H |
| `Figure5_peptide_analysis.Rmd` | Peptide-level filtering, MAPT peptide inspection, mapping, and plotting | Figure 5 |

## Data availability

This repository contains analysis code, not restricted participant-level data.
The code can be inspected without the input data; several elements of execution require authorised access to the restricted metadata. This is subject to the applicable data-use agreements, through [Rush Alzheimer's Disease Center](https://www.radc.rush.edu) and [Synapse](https://www.synapse.org).

## Protein matrix QC workflow

The script imports DIA-NN protein-group intensities, maps runs to the experimental layout, removes excluded samples and contaminants, filters proteins with at least 70% missingness, removes entries listing multiple accessions or genes, transforms intensities to log2 scale, and applies MAI imputation.

Inputs:
- `report.pg_matrix.tsv`: raw DIA-NN output of protein-group intensities.
- `mapnamesreport.csv`: experimental sample mapping.
- `contaminants.csv`: contaminant exclusion list.

Note that the supplementary PCA requires some restricted external variables to plot all panels.
The raw protein group matrix, and resulting processed SummarizedExperiment containing `log2Quant` and `MAI_imputed` assays can be downloaded from proteomeXchange.

## Subcellular profiles workflow

`Figure2_subcell_profiles.Rmd` uses the MAI-imputed protein matrix stored in the SummarizedExperiment object to assess subcellular fractionation and reproducibility across individuals. 

Analyses include:
- Hierarchical clustering of scaled profiles and GO:CC enrichment (Figure 2A).
- Curated marker profiles (Figure 2Ci–ii and Supplementary Figure S2).
- t-SNE of scaled mean proportional profiles, coloured by marker compartment (Figure 2Ciii).
- Pairwise subject correlations of log2-imputed protein intensities across fractions, by disease group (Figure 2D).

Inputs:
- `MAI_prot_model_new.rds`: processed SummarizedExperiment from the protein matrix QC workflow.
- `merged_markers_filtered.csv`: final curated marker assignments, with `Protein.Group` and `markers` columns.
- `sample_category.csv`: sample-to-disease-group mapping, with `Subject` and `Group` columns (CON`,`RES`,`AD-DEM`. This allows execution of chunk to generate Figure 2D without any restricted external variables.

## Disease and resilience modelling workflows

### Disease-associated protein fractionation — Figure 3A–E

`Figure3_Modelling_categorical.Rmd` prepares protein profiles from the imputed and unimputed assays and identifies those eligible for localisation modelling. Relative predictor importance is compared using drop-one-out changes in marginal R2. Disease-associated (CON/RES/AD-DEM) differences in fractionation profiles are identified using mixed models with disease-group x fraction interactions.

Results are compared between imputed and unimputed analyses to identify robust associations. Example protein profiles and GO Cellular Component enrichment analysis is also included.

Required inputs are:
- `MAI_prot_model_new.rds` 
- `demogs.xlsx`
- `sample_category.csv` can be used in place of demogs for the plotting of fractionation profiles. This allows execution of chunk to plot figures 3B and 3D without requiring restricted external variables.
 
### Resilience-associated protein fractionation — Figure 3F–H

`Figure3_Modelling_Resilience.Rmd` fits a longitudinal mixed model of cognition using the full available ROSMAP cohort, retaining individuals with at least three lognitudinal cognitive observations and the required covariates.
These resilience scores are linked to the proteomics cohort and used in the localisation mixed model to identify a resilience x fraction interaction effect. Results from imputed and unimputed protein intensities are compared to identify robust associations. The script also includes cognitive trajectory and resilience visualisations.

Additional required input

- `dataset_longitudinal.xlsx`: restricted longitudinal cognition data.
- `dataset_cross_sectional.xlsx`: restricted pathology and demographic covariates.

These can be requested through [Rush Alzheimer's Disease Center](https://www.radc.rush.edu) and [Synapse](https://www.synapse.org), subject to the applicable data-use agreements.

Note that the protein preparation and eligibility sections of `Figure3_Modelling_categorical.Rmd` need to be run before the resilience workflow. Intermediate data and model outputs are not included in this repository as they contain restricted participant-level variables.

## MAPT peptides workflow

`Figure5_peptide_analysis.Rmd` imports DIA-NN precursor intensities (These can be downloaded from ProteomeXchange), maps runs to the experimental layout, removes excluded samples and contaminants, and filters by missingness to obtain a clean PR (precursor) matrix. Unmodified precursors are median-aggregated. MAPT peptides with <5% missingness are mapped to the canonical 2N4R sequence and plotted by disease-profile. The accompanying reference cytoskeletal and cytosolic protein-level compartment profiles can be obtained in the figure 2Ci section of `Figure2_subcell_profiles.Rmd`. 

Inputs:
- `report.pr_matrix.tsv`: raw DIA-NN output of precursor intensities.
- `sample_category.csv`: sample-to-disease-group mapping, with `Subject` and `Group` columns (CON`,`RES`,`AD-DEM`. This allows execution of chunk to generate Figure  5 without any restricted external variables.

## Software

Analyses use R. Packages referenced include: tidyverse, data.table, readr, SummarizedExperiment, S4Vectors, MAI, dplyr, tidyr, tibble, purrr, stringr, forcats, ggplot2, patchwork, viridis, scales, paletteer, ComplexHeatmap, circlize, clusterProfiler, org.Hs.eg.db, tidytext, and Rtsne. Generating R Markdown also requires rmarkdown and knitr.

## Manuscript and citation

The full manuscript title, author list, citation, and contact details are detailed in the submission. Further scripts will be listed in the analysis index as they are prepared.
