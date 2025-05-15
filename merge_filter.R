# Combine the AWS csv files into one csv file for each biomarker

# Directory containing results
path <- "Results/"

# Function to merge all chunk files for a given biomarker
merge_filter_biomarker_files <- function(biomarker, threshold) {
  
  files <- list.files(path = paste0(path, biomarker), pattern = paste0(biomarker, "_EWAS_.*\\.csv"), full.names = TRUE)
  
  biomarker_data <- do.call(rbind, lapply(files, read.csv))
  
  # Write merged file
  output_file <- paste0(path, biomarker, "_EWAS_Results_Full.csv")
  write.csv(biomarker_data, file = output_file, row.names = FALSE)
  
  
  # Apply FDR correction using the Benjamini-Hochberg method
  biomarker_data$FDR_pval <- p.adjust(biomarker_data$pval, method = "BH")
  
  # Filter significant results based on FDR threshold
  significant_results <- biomarker_data[biomarker_data$FDR_pval < threshold, ]
  
  # Save significant results file
  significant_output_file <- paste0(path, biomarker, "_EWAS_Results_Significant_10.csv")
  write.csv(significant_results, file = significant_output_file, row.names = FALSE)
  

}
threshold <- 0.10
biomarkers <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
lapply(biomarkers, merge_filter_biomarker_files, threshold = threshold)