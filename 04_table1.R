# Creating Table 1 for descriptive statistics 
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries 
library(tidyverse)
library(tableone)
library(survey)
library(forcats)
library(openxlsx)

# Load the raw data 
dat <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_pheno.csv")


dat <- dat %>%
  mutate(AGE_CAT = case_when(
    AGE < 60 ~ "<60",
    AGE >= 60 & AGE <= 70 ~ "60−70",
    AGE > 70 ~ ">70"
  )) %>% 
  mutate(AGE_CAT = factor(AGE_CAT, levels = c("<60", "60−70", ">70"))) %>%
  mutate(BMI_CAT = case_when(
    BMI < 18.5 ~ "<18.5",
    BMI >= 18.5 & BMI <= 25 ~ "18.5−25",
    BMI > 25 ~ ">25"
  )) %>% 
  mutate(APOE_GENOTYPE = as.character(APOE_GENOTYPE)) %>% 
  mutate(E2 = grepl("2", APOE_GENOTYPE), E3 = grepl("3", APOE_GENOTYPE), E4 = grepl("4", APOE_GENOTYPE)) %>% 
  mutate(E2 = factor(E2, levels = c("TRUE", "FALSE")), 
         E3 = factor(E3, levels = c("TRUE", "FALSE")), 
         E4 = factor(E4, levels = c("TRUE", "FALSE"))) %>% 
  mutate(APOE_GENOTYPE = fct_recode(APOE_GENOTYPE,
                                    "ε2/ε2" = "22",
                                    "ε2/ε3" = "23",
                                    "ε2/ε4" = "24", 
                                    "ε3/ε3" = "33",
                                    "ε3/ε4" = "34",
                                    "ε4/ε4" = "44", 
                                    "NA" = "-9", 
                                    "NA" = "99"))

vars <-  c("AGE",
           "AGE_CAT",
           "SEX", 
           "BMI",
           "BMI_CAT",
           "CENTER", 
           "ABETA40", 
           "ABETA42", 
           "PTAU181", 
           "NFLIGHT", 
           "GFAP", 
           "APOE_GENOTYPE", 
           "E2", 
           "E3", 
           "E4")

survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_INCA , data=dat)

tbl1_unweighted <- CreateTableOne(vars = vars, data = dat, strata = "MCI", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

tbl1_weighted <- svyCreateTableOne(vars = vars, data = survey_obj, strata = "MCI", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

tbl1 <- tbl1_unweighted

# keep only "true" levels for APOE allele carrier rows
allele_inds_to_remove <- which(tbl1[,"level"] == "FALSE")
tbl1 <- tbl1[-c(allele_inds_to_remove),]
tbl1_unweighted <- tbl1_unweighted[-c(allele_inds_to_remove),]
tbl1_weighted <- tbl1_weighted[-c(allele_inds_to_remove),]
allele_inds_to_update <- which(tbl1[,"level"] == "TRUE")
tbl1[allele_inds_to_update,"level"] = ""


col_inds_to_update <- 2:4 # all columns except "level"
row_inds_to_update <- c(1, 3:5, 6:7, 9:11, 12:15, 21:27, 28:30) 
# rows corresponding to N, Age (categorical), Sex, BMI (categorical), Center, APOE genotype, allele carriers

# create function to remove any leading/trailing spaces around the parentheses
fix_spacing <- function(x) {
  gsub("\\s*\\(", " (", gsub("\\s*\\)", ") ", gsub("\\s*\\(\\s*", " (", trimws(x)))) 
}

# update tbl1_comb with the percentages from the weighted table
for (i in col_inds_to_update){
  counts <- sapply(tbl1_unweighted[,i], function(x){
    strsplit(x, split = "(", fixed = TRUE)[[1]][1]
  })
  percent <- sapply(tbl1_weighted[,i], function(x){
    paste0("(",  strsplit(x, split = "(", fixed = TRUE)[[1]][2])
  })
  percent[which(names(percent) == "n")] = ""
  tbl1[,i] <- paste0(counts, percent)
  tbl1[,i] <- fix_spacing(tbl1[, i])
}
rownames(tbl1) = c("N", 
                        "Age (mean (SD))", 
                        "Age (%)", "", "", 
                        "Sex (%)", "",
                        "BMI (mean (SD))", 
                        "BMI (%)", "", "", 
                        "Center", "", "", "",
                        "Aꞵ40 (mean (SD))", "Aꞵ42 (mean (SD))", "pTau-181 (mean (SD))", "NfL (mean (SD))", "GFAP (mean (SD))", 
                        "APOE genotype (%)","",  "", "", "", "", "", 
                        "ε2 carrier (%)", "ε3 carrier (%)", "ε4 carrier (%)"
                        )

tbl1 = cbind(rownames(tbl1), tbl1)


write.xlsx(tbl1, file = "2024_EWAS_ATN/ATN_EWAS/Results/20250515_table1.xlsx", overwrite = TRUE)

