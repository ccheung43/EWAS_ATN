# Handling ATN biomarker outliers
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(survey)
library(ggplot2)
library(patchwork)

# Read in data from csv file 
dat = read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_batch_adj.csv") %>% 
  mutate(E4 = if_else(grepl("4", APOE_GENOTYPE), 1, 0))

# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)


biomarker_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
residuals_df <- data.frame(matrix(ncol = 4, nrow = nrow(dat)))
colnames(residuals_df) <- biomarker_vars
for (var in biomarker_vars){
  mod <- svyglm(formula = as.formula(paste0(var, "~SEX+AGE+BMI+CENTER+PC1+PC2+PC3+PC4+PC5+CD4T+NK+Bcell+Mono+Gran")), 
                design = svy_obj, family = gaussian())
  
  # Store residuals in a df
  residuals_df[var][which(!is.na(dat[var])),] <- residuals(mod)
} 

# Get residual stats (including quartiles and IQR) to assess upper and lower bounds of extreme outliers 
residuals_stats <- residuals_df %>% 
  pivot_longer(1:4, names_to = "outcome", values_to = "residuals") %>%
  group_by(outcome) %>% 
  summarize(Q1 = quantile(residuals, 0.25, na.rm=TRUE), 
            M = quantile(residuals, 0.5, na.rm = TRUE),
            Q3 = quantile(residuals, 0.75, na.rm=TRUE), 
            Avg = mean(residuals, na.rm = TRUE), 
            SD= sd(residuals, na.rm = TRUE)) %>% 
  mutate(IQR = Q3-Q1) %>% 
  mutate(lower_bound = Q1 - 3*IQR, upper_bound = Q3 + 3*IQR) %>% # extreme outliers defined as outside of 3xIQR
  mutate(outcome = factor(outcome, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) 

# Histogram of residuals before removal of outliers 
residuals_hist_outliers <- residuals_df %>% pivot_longer(1:4, names_to = "outcome", values_to = "residuals") %>%
  mutate(outcome = factor(outcome, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) %>% 
  ggplot(aes(x = residuals)) +
  geom_histogram(bins = 30) +
  geom_vline(data = residuals_stats, aes(xintercept = M), color = "red", linetype = "solid", linewidth = 1) +
  geom_vline(data = residuals_stats, aes(xintercept = lower_bound), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(data = residuals_stats, aes(xintercept = upper_bound), color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ outcome, scales = "free", nrow = 1, 
             labeller = labeller(outcome = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP"))) +
  theme_minimal(base_size = 14) +
  labs(title = "Histogram of Residuals for Each Biomarker", subtitle = "Before outlier removal", x = "", y = "Count")

# Removal of outlier datapoints from dat 
for (var in biomarker_vars){
  # Define and remove outliers 
  lower_bound <- residuals_stats$lower_bound[which(residuals_stats$outcome == var)]
  upper_bound <- residuals_stats$upper_bound[which(residuals_stats$outcome == var)]
  outliers = which(residuals_df[var] > upper_bound | residuals_df[var] < lower_bound)
  dat[outliers, var] = NA 
  residuals_df[outliers, var] = NA
}

# Histogram of residuals after removal of outliers 
residuals_hist_removed <- residuals_df %>% pivot_longer(1:4, names_to = "outcome", values_to = "residuals") %>%
  mutate(outcome = factor(outcome, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) %>% 
  ggplot(aes(x = residuals)) +
  geom_histogram(bins = 30) +
  geom_vline(data = residuals_stats, aes(xintercept = M), color = "red", linetype = "solid", linewidth = 1) +
  facet_wrap(~ outcome, scales = "free", nrow = 1, 
             labeller = labeller(outcome = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP"))) +
  theme_minimal(base_size = 14) +
  labs(subtitle = "After outlier removal", x = "Residuals", y = "Count")



# QQ plot of residuals after removal of outliers
residuals_qq <- residuals_df %>% pivot_longer(1:4, names_to = "outcome", values_to = "residuals") %>%
  mutate(outcome = factor(outcome, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) %>% 
  ggplot(aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line(color = "red", linewidth = 1) +
  facet_wrap(~ outcome, scales = "free", nrow = 1, 
             labeller = labeller(outcome = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP"))) +
  theme_minimal(base_size = 14) +
  labs(title = "Q-Q Plot of Residuals for Each Model", subtitle = "After outlier removal", x = "Theoretical Quantiles", y = "Sample Quantiles")


residuals_plots <- residuals_hist_outliers / residuals_hist_removed / residuals_qq

# Include only indiviudals with at least one ATN biomarker (after batch correction and outlier removal) 
outlier_removed_dat <- dat %>% filter(!if_all(all_of(biomarker_vars), is.na)) # removed 8 individuals from the dataset

# 2381 individuals
write.csv(outlier_removed_dat, file = "2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_outliers_removed.csv", row.names = FALSE)

