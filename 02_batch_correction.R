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
# Convert plate position to linear position 
dat <- dat %>%
  mutate(PLATE_POS_LINEAR = as.numeric(factor(PLATE_POS, levels = paste0(rep(LETTERS[1:8], each = 12), sprintf("%02d", 1:12)))))

# Loop through each ATN biomarkers
biomarker_vars = c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")
for (biomarker in biomarker_vars) {
  # Set the current ATN biomarker as outcome and create a temporary data frame
  temp_dat <- dat %>% filter(!is.na(dat[[biomarker]]))
  
  fit.warn <- tryCatch(
    {list(lmer(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER + (1 | PLATE_NUM) + PLATE_POS_LINEAR")), data = temp_dat), "MixedOK")},
    warning = function(Warn) {
      print(paste("MY_WARNING:  ", Warn))
      fit <- lm(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER + PLATE_POS_LINEAR")), data = temp_dat)
      return(list(fit, "Warn"))
    },
    error = function(err) {
      print(paste("MY_ERROR:  ", err))
      fit <- lm(as.formula(paste0(biomarker, " ~ AGE + SEX + BMI + CENTER + PLATE_POS_LINEAR")), data = temp_dat)
      return(list(fit, "err"))
    }
  )
  
  fit <- fit.warn[[1]]
  
  # Remove plate position effect
  model.mat <- model.matrix(fit)
  PLATE_POS_LINEAR.mat <- model.mat[, "PLATE_POS_LINEAR"]
  PLATE_POS_LINEAR.effects <- summary(fit)$coefficients["PLATE_POS_LINEAR", "Estimate"]
  
  # Create a vector of adjusted biomarker values for the original dataset
  dat[[biomarker]][!is.na(dat[[biomarker]])] <- dat[[biomarker]][!is.na(dat[[biomarker]])] - PLATE_POS_LINEAR.mat * PLATE_POS_LINEAR.effects
  
  # If the mixed model worked, remove the random plate number effect
  if (fit.warn[[2]] == "MixedOK") {
    dat$PLATE_RAND <- ranef(fit)$PLATE_NUM[match(dat$PLATE_NUM, rownames(ranef(fit)$PLATE_NUM)), 1]
    dat[[biomarker]] <- dat[[biomarker]] - dat$PLATE_RAND
  }
  
  dat <- dat %>% select(-PLATE_RAND)
}

dat_batch_adj <- dat %>% 
  select(-c(PLATE_NUM, PLATE_POS, PLATE_POS_LINEAR)) %>% 
  mutate(ABETA4240 = ABETA42/ABETA40)

write.csv(dat_batch_adj, "2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_batch_adj.csv", row.names = FALSE)


