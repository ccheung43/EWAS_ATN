# Fitting a regression model with basic covariates, and subtracting the estimated effect of plate number and plate position
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(lme4)

# Read in data from csv file 
data_file = "2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_pheno.csv"
dat <- read.csv(data_file)

# Removing ATN biomarker batch effects

# Loop through each ATN biomarkers
biomarker_vars = c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")
for (biomarker in biomarker_vars) {
  # Set the current ATN biomarker as outcome and create a temporary data frame
  temp_dat <- dat %>% filter(!is.na(dat[[biomarker]]))
  
  fit.warn <- tryCatch(
    {list(lmer(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER + (1 | PLATE_NUM) + (1 | PLATE_POS)")), data = temp_dat), "MixedOK")},
    warning = function(Warn) {
      print(paste("MY_WARNING:  ", Warn))
      fit <- lm(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER")), data = temp_dat)
      return(list(fit, "Warn"))
    },
    error = function(err) {
      print(paste("MY_ERROR:  ", err))
      fit <- lm(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER")), data = temp_dat)
      return(list(fit, "err"))
    }
  )
  
  fit <- fit.warn[[1]]
  
  if (fit.warn[[2]] == "MixedOK") {
    ranef_list <- ranef(fit)
    
    # Subtract random effect of PLATE_NUM
    dat$PLATE_RAND <- ranef_list$PLATE_NUM[match(dat$PLATE_NUM, rownames(ranef_list$PLATE_NUM)), 1]
    
    # Subtract random effect of PLATE_POS
    dat$PLATE_POS_RAND <- ranef_list$PLATE_POS[match(dat$PLATE_POS, rownames(ranef_list$PLATE_POS)), 1]
    
    dat[[biomarker]][!is.na(dat[[biomarker]])] <- 
      dat[[biomarker]][!is.na(dat[[biomarker]])] - dat$PLATE_RAND[!is.na(dat[[biomarker]])] - dat$PLATE_POS_RAND[!is.na(dat[[biomarker]])]
  }
  
  # Remove temporary columns
  dat <- dat %>% select(-PLATE_RAND, -PLATE_POS_RAND)
}

dat_batch_adj <- dat %>% 
  select(-c(PLATE_NUM, PLATE_POS)) %>% 
  mutate(ABETA4240 = ABETA42/ABETA40)

write.csv(dat_batch_adj, "2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_batch_adj.csv", row.names = FALSE)


