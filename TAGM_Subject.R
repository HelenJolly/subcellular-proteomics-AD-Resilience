
library("QFeatures")
library("pRoloc")
library("bandle")
library("tidyverse")
library("dbscan")
library("tidyr")
library("dplyr")
library(ggplot2)
library(rlang)
library(readxl)
require(mvtnorm)
require(dplyr)
require(ggrepel)
library(SummarizedExperiment)

# Step 1: Load the SE object
se_log2_07 <- readRDS("MAI_prot_model_new.rds")
#Summarised experiment, global filter of <0.7 missingness, MAI_imputed

#extract the imputed assay data and backtransform
MAI_backtransformed <- 2^(assay(se_log2_07, "MAI_imputed")) #linear scale

#Extract row (Gene) metadata and col (subject/run) metadata
rowmeta <- as.data.frame(rowData(se_log2_07)) %>% dplyr::select(Genes) %>% tibble::rownames_to_column("Protein.Group")
colmeta <- as.data.frame(colData(se_log2_07)) %>% dplyr::select(Subject, Fraction) %>% tibble::rownames_to_column("RunID")
subjects <- colmeta %>% filter(Subject !="G2") %>% pull(Subject) %>% unique() #extract 74 subjects (G2 excluded due to incomplete fractionation profile)

######################### markers and proteins of interest ###########################


markers_to_keep <- read_csv(input_file("merged_markers_filtered.csv"))
marker_vec <- setNames(markers_to_keep$markers, markers_to_keep$Protein.Group); marker_vec[1:10] #protein group

# Prepare list to store predictions per subject
all_subject_predictions <- list()

# Loop through each subject
for (i in subjects) { 
  cat("\n=========== Processing subject:", i, "===========\n")
  cat("Subsetting and ordering expression data...\n")
  assay_df <- as.data.frame(MAI_backtransformed) %>%
    dplyr::select(matches(paste0("-", i, "-"))) %>%
    {
      fraction_index <- str_extract(colnames(.), "(?<=-)([1-7])(?=p$)")
      col_info <- data.frame(colname = colnames(.), fraction = as.numeric(fraction_index), stringsAsFactors = FALSE)
      col_order <- col_info$colname[order(col_info$fraction)]
      .[, col_order]
    } %>% tibble::rownames_to_column("Protein.Group") #or Genes
  
  cat("Creating QFeatures object...\n")
  qf <- readQFeatures(assayData = assay_df, quantCols = 2:8, name = "MAI_prots")
  
  cat("Applying sum normalisation...\n")
  qf <- normalize(qf, i = "MAI_prots", name = "MAI_prots_norm", method = "sum")
  
  cat("Adding markers...\n")
  msn <- as(qf[["MAI_prots_norm"]], "MSnSet") %>%
    addMarkers(markers = marker_vec, fcol = "Protein.Group", mcol = "markers")
  featureNames(msn) <- fData(msn)$Protein.Group 
  fData(msn)$markers <- fData(msn)$markers
  
  cat("Training TAGM-MAP model...\n")
  set.seed(2)
  map_params <- tagmMapTrain(msn, fcol = "markers", numIter = 100)
  
  ## Optional segment (run-time much longer)
  #cat("Plotting log-posterior convergence...\n")
  #plot(logPosteriors(map_params), type = "b", col = "blue", cex = 0.3, ylab = "log-posterior", xlab = "Iteration")
  #title(main = paste("TAGM-MAP training - Subject", i))
  
  cat("Predicting localisations using MAP parameters...\n")
  msn <- tagmMapPredict(
    object      = msn,
    params      = map_params,
    fcol        = "markers", probJoint   = TRUE, probOutlier = TRUE) # Original prediction options retained.
  
  cat("Classifications (before thresholding):\n")
  print(table(fData(msn)$tagm.map.allocation))
  
  cat("Summary of compartment localisation probabilities:\n")
  print(summary(fData(msn)$tagm.map.probability))
  
  cat("Summary of outlier probabilities:\n")
  print(summary(fData(msn)$tagm.map.outlier))
  
  # Compute overall probability
  cat("Computing overall probability...\n")
  tagm_prob <- fData(msn)[, "tagm.map.probability"]
  tagm_out  <- 1 - fData(msn)[, "tagm.map.outlier"]
  fData(msn)[, "overall_prob"] <- tagm_prob * tagm_out
  
  cat("Applying 0.95 threshold to predictions...\n")
  msn <- getPredictions(msn, 
                        fcol = "tagm.map.allocation", 
                        scol = "overall_prob", 
                        mcol  = "markers",
                        t     = 0.95)
  
  cat("turning into dataframe...\n")
  fdat <- fData(msn) %>%
    as.data.frame() %>%
    mutate(Subject = i)
  
  ##check
  cat("checking range is 1...\n")
  fdat %>%
    mutate(total_prob = rowSums(dplyr::select(., matches("^tagm\\.map\\.joint")))) %>%
    summarise(
      min_prob = min(total_prob),
      max_prob = max(total_prob))
  
  cat("Storing predictions for subject", i, "\n")
  all_subject_predictions[[i]] <- fdat
}

#Save as an RDS 
#saveRDS(all_subject_predictions, "all_subject_predictions_raw.rds")

#Save as a dataframe
all_subject_predictions_expanded <- lapply(all_subject_predictions, function(df) {
  
  joint_mat <- df$tagm.map.joint #extract joint (outlier weghted probability) matrix
  
  #convert matrix to a dataframe with proper names
  joint_df <- as.data.frame(joint_mat)
  colnames(joint_df) <- colnames(joint_mat)  #already correct
  
  #bind to the rest of the columns *except* the matrix (so there's no name clash)
  df_out <- cbind( df[, !names(df) %in% "tagm.map.joint"], joint_df )
  
  return(df_out)
})

localisation_predictions <- bind_rows(all_subject_predictions_expanded)
## Write out concatenated dataset
#write.csv(localisation_predictions, "localisation_Apr26.csv"))
