library(haven)
library(corrplot)
library(RColorBrewer)
library(survey)
library(tidyverse)


# Read in the ATN + phenotype data
dat <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_outliers_removed.csv")
# Read in the MRI structure data
struc <- read_sas("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Datasets/SOL-INCA-MRI/Datasets/inca_mri_struc_inv2.sas7bdat") %>% 
  rename(HCHS_ID = ID) %>% 
  unique() %>% 
  group_by(HCHS_ID) %>% # remove any rows with duplicated HCHS IDs
  filter(n() == 1) %>%
  mutate(HCHS_ID = as.integer(HCHS_ID)) 
mri <- merge(dat, struc, by = "HCHS_ID")
# Read in the sample weights
dt <- read_sas('//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Datasets/SOL-INCA-MRI/Datasets/inca_mri_part_derv_inv2.sas7bdat') %>% 
  rename(HCHS_ID = ID) %>%
  unique() %>% 
  group_by(HCHS_ID) %>% # remove any rows with duplicated HCHS IDs
  filter(n() == 1) %>%
  mutate(HCHS_ID = as.integer(HCHS_ID)) %>% 
  select(c("HCHS_ID", "WEIGHT_NORM_OVERALL_INCAMRI"))
mri <- merge(dt, mri, by = "HCHS_ID")
mri$CORTICAL_TGV <- mri$FRONTAL_CORTICAL_GRAY + mri$OCCIPITAL_CORTICAL_GRAY + mri$PARIETAL_CORTICAL_GRAY + mri$TEMPORAL_CORTICAL_GRAY

write.csv(mri, "2024_EWAS_ATN/ATN_EWAS/Data/20250508_MRI.csv", row.names = FALSE)


mri <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250508_MRI.csv")

# Correlation plots of MRI data
cor.mri <- mri %>% 
  select(c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP", 
                            "CEREBRUMV2_TCV", "CEREBRUMV2_TCB", "CEREBRUMV2_WMH", "HIPPO", "CORTICAL_TGV", 
                            "BMI", "AGE")) %>% 
  rename(CEREB_TCV = CEREBRUMV2_TCV, CEREB_TCB = CEREBRUMV2_TCB, CEREB_WMH = CEREBRUMV2_WMH)

M <-cor(cor.mri,use = "complete.obs")
corrplot(M, type="lower",mar=c(0,0,1,0),addCoef.col = 'black',
         col=brewer.pal(n=8, name="RdYlBu"))

# # Normalize to TCV using linear regression
# cc <- complete.cases(mri$CEREBRUMV2_TCB, mri$CEREBRUMV2_TCV) 
# tcb <- lm(CEREBRUMV2_TCB ~ CEREBRUMV2_TCV, data=mri)
# mri$tcb.res <- NA 
# mri$tcb.res[cc] <- residuals(tcb)
# 
# cc <- complete.cases(mri$HIPPO, mri$CEREBRUMV2_TCV) 
# hip <- lm(HIPPO ~ CEREBRUMV2_TCV, data=mri)
# mri$hip.res <- NA 
# mri$hip.res[cc] <- residuals(hip)
# 
# mri$CEREBRUMV2_WMH[mri$CEREBRUMV2_WMH == 0] <- NA
# mri$log_WMH <- log(mri$CEREBRUMV2_WMH)
# cc <- complete.cases(mri$log_WMH, mri$CEREBRUMV2_TCV) 
# wmh <- lm(log_WMH ~ CEREBRUMV2_TCV, data=mri)
# mri$wmh.res <- NA
# mri$wmh.res[cc] <- residuals(wmh)


# Normalize to TCV using survey-weighted linear regression
mri$CEREBRUMV2_WMH[mri$CEREBRUMV2_WMH == 0] <- NA
mri$log_WMH <- log(mri$CEREBRUMV2_WMH)
survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~mri$WEIGHT_NORM_OVERALL_INCAMRI, data=mri)
options(survey.lonely.psu="adjust")

cc <- complete.cases(mri$CEREBRUMV2_TCB, mri$CEREBRUMV2_TCV) 
tcb <- svyglm(CEREBRUMV2_TCB ~ CEREBRUMV2_TCV, design = subset(survey_obj, cc), family = "gaussian")
mri$tcb.res <- NA 
mri$tcb.res[cc] <- residuals(tcb)

cc <- complete.cases(mri$HIPPO, mri$CEREBRUMV2_TCV) 
hip <- svyglm(HIPPO ~ CEREBRUMV2_TCV, design = subset(survey_obj, cc), family = "gaussian")
mri$hip.res <- NA 
mri$hip.res[cc] <- residuals(hip)

cc <- complete.cases(mri$log_WMH, mri$CEREBRUMV2_TCV) 
wmh <- svyglm(log_WMH ~ CEREBRUMV2_TCV, design = subset(survey_obj, cc), family = "gaussian")
mri$wmh.res <- NA 
mri$wmh.res[cc] <- residuals(wmh)

cc <- complete.cases(mri$log_WMH, mri$CEREBRUMV2_TCV) 
wmh <- svyglm(log_WMH ~ CEREBRUMV2_TCV, design = subset(survey_obj, cc), family = "gaussian")
mri$wmh.res <- NA 
mri$wmh.res[cc] <- residuals(wmh)

cc <- complete.cases(mri$CORTICAL_TGV, mri$CEREBRUMV2_TCV) 
tgv <- svyglm(CORTICAL_TGV ~ CEREBRUMV2_TCV, design = subset(survey_obj, cc), family = "gaussian")
mri$tgv.res <- NA 
mri$tgv.res[cc] <- residuals(tgv)

# Save normalized data 
write.csv(mri, "2024_EWAS_ATN/ATN_EWAS/Data/20250508_MRI_normalized.csv", row.names = FALSE)

