# Function to run survey-weighted GLM on CpG methylation data

runsvy <- function(cpg_list, ppath, mpath, biomarkers) {
  # Read in phenotype data
  pdata <- read.csv(ppath) %>% 
    mutate(CENTER = factor(CENTER,levels = c("B","C","M","S"),labels=c("Bronx","Chicago", "Miami", "San Diego")),
           SEX = factor(SEX,levels = c("F","M"),labels=c("Female","Male")))
  
  # Read in methylation data for the specified CpG sites
  # subset columns to match phenotypic data, transpose, convert to dataframe
  col_keep <- c('CpG',pdata$Row_names) # select the CpG column and samples you want to extract
  mdata <- open_dataset(mpath) %>% # open connection to parquet database
    select(all_of(col_keep)) %>% 
    filter(CpG %in% cpg_list) %>% 
    as.data.frame() %>% 
    remove_rownames() %>% 
    column_to_rownames('CpG') %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Row_names") 
  
  # Merge methylation and phenotype data 
  sdata <- mdata %>% left_join(pdata, by = "Row_names") 
    
  
  results_list <- setNames(lapply(biomarkers, function(biomarker) {
    
    sdata_biomarker <- sdata %>% filter(!is.na(.data[[biomarker]]))  # Remove rows with missing biomarker values
    
    # Define survey design with stratification, clustering, and weights
    survey_obj <- svydesign(id = ~PSU_ID, strata = ~STRAT, weights = ~WEIGHT_NORM_OVERALL_INCA, 
                            nest = TRUE, data = sdata_biomarker)
    
    # Identify CpG sites in the dataset that match the provided CpG list
    cpg <- colnames(sdata_biomarker)[colnames(sdata_biomarker) %in% cpg_list]

    # Construct regression formulas for each CpG
    mod_formulas <- lapply(cpg, function(x) { 
      as.formula(paste0(biomarker, " ~ ", x, " + AGE + SEX + BMI + CENTER + PC1 + PC2 + PC3 + PC4 + PC5 + CD4T + NK + Bcell + Mono + Gran")) 
    }) 
    
    # Fit survey-weighted GLM models
    mods <- lapply(mod_formulas, function(x) svyglm(x, design = survey_obj, family = "gaussian"))
    
    # Extract coefficients, SEs, p-values, and confidence intervals
    coef.m <- lapply(mods, function(x) summary(x)$coefficients)
    beta <- sapply(coef.m, function(x) x[2, 1])  # Beta coefficients
    p.val <- sapply(coef.m, function(x) x[2, 4])  # P-values
    sd <- sapply(coef.m, function(x) x[2, 2])  # Standard errors
    
    ci <- lapply(mods, function(x) confint(x))
    lower <- sapply(ci, function(x) x[2, 1])  # Lower CI bound
    upper <- sapply(ci, function(x) x[2, 2])  # Upper CI bound
    
    # Store results in data frame
    data.frame(
      CpG = cpg, 
      coef = beta, 
      sd = sd,
      pval = p.val, 
      lower = lower,
      upper = upper
    )
  }), biomarkers)  # **Set names to biomarkers**
  
  return(results_list)
}
