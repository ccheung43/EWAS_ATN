# Load required packages 
library(tidyverse)
library(dplyr)
library(survey)
library(future)
library(furrr)
library(rlang)
system("echo $LD_LIBRARY_PATH") 
library(arrow)

# Source Meth Parquet Functions 
#source("Code/cpg_extract.R")
#source("Code/create_methlist.R")
source("Code/runsvy.R")

args <- commandArgs(trailingOnly = TRUE)
blk <- as.numeric(args[1])

cpg_annotation <- read.csv("Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7)

# Split data into chunks
chunk_size <- 20000
chunk_labels  <- rep(1:ceiling(nrow(cpg_annotation)/chunk_size),each=chunk_size)[1:nrow(cpg_annotation)]
chunks <- split(cpg_annotation$Name,chunk_labels)
selected_chunk <- chunks[[blk]]

# # Running not in Parellel
# # Run runsvy function on the selected chunk
# biomarkers = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
# res <- runsvy(chunks, atn_phen, biomarkers) # when not running in parallel
# lapply(names(res), function(biomarker) {
#   path = "/home/ec2-user/EBS4T/Projects/2024_EWAS_ATN/Results/"
#   file_name = paste0(path, biomarker, "/", biomarker, "_EWAS_", blk, ".csv")
#   write.csv(res[[biomarker]], file = file_name, row.names = FALSE)
# })

# Running in Parallel
# Split data into subchunks (for parallelization)
subchunk_size <- 4000
subchunk_labels <- rep(1:ceiling(length(selected_chunk)/subchunk_size),each=subchunk_size)[1:length(selected_chunk)]
subchunks <- split(selected_chunk, subchunk_labels)

# Run runsvy function on the selected chunk
biomarkers = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
ppath = "Data/20250224_ATN_outliers_removed.csv"
mpath <-"Data/Combat_V02/"
options(future.globals.maxSize = 1e9) 
plan(multicore, workers=5) 
res <- future_map(subchunks, ~ runsvy(.x, ppath, mpath, biomarkers), .options = furrr_options(seed = TRUE), .progress = TRUE)


lapply(names(res[[1]]), function(biomarker) { 
  biomarker_res <- lapply(names(res), function(subchunk) { 
    res[[subchunk]][[biomarker]]
  }) 
  biomarker_res <- do.call(rbind, biomarker_res)
  path = "Results/"
  file_name = paste0(path, biomarker, "/", biomarker, "_EWAS_", blk, ".csv")
  write.csv(biomarker_res, file = file_name, row.names = FALSE)
}) 

