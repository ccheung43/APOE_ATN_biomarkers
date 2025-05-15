# Handling outliers (plus supplementary figure 1?)
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(ggplot2)
library(survey)
library(gridExtra)
library(Cairo)
library(patchwork)

# Read in data from csv file 
dat = read.csv("2024_APOE_ATN_biomarkers/Data/20241031_data_batch_adj.csv")
plot_file = "2024_APOE_ATN_biomarkers/Results/20241031_normality_plots.png"
outlier_removed_data_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"

# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)

outcome_vars <- c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")
residuals_df <- data.frame(matrix(ncol = 4, nrow = nrow(dat)))
colnames(residuals_df) <- outcome_vars
for (var in outcome_vars){
  mod <- svyglm(formula = as.formula(paste0(var, "~E2_add+E4_add+SEX+AGE+CENTER")), 
                design = svy_obj, family = gaussian())
  
  # Store residuals in a df
  residuals_df[var][which(!is.na(dat[var])),] <- residuals(mod)
} 

residuals_stats <- residuals_df %>% 
  pivot_longer(1:4, names_to = "outcome", values_to = "residuals") %>%
  group_by(outcome) %>% 
  summarize(Q1 = quantile(residuals, 0.25, na.rm=TRUE), 
            M = quantile(residuals, 0.5, na.rm = TRUE),
            Q3 = quantile(residuals, 0.75, na.rm=TRUE), 
            Avg = mean(residuals, na.rm = TRUE), 
            SD= sd(residuals, na.rm = TRUE)) %>% 
  mutate(IQR = Q3-Q1) %>% 
  mutate(lower_bound = Q1 - 3*IQR, upper_bound = Q3 + 3*IQR) %>% 
  mutate(outcome = factor(outcome, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) 

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


normality_pvals <- c()
for (var in outcome_vars){
  # Define and remove outliers 
  lower_bound <- residuals_stats$lower_bound[which(residuals_stats$outcome == var)]
  upper_bound <- residuals_stats$upper_bound[which(residuals_stats$outcome == var)]
  outliers = which(residuals_df[var] > upper_bound | residuals_df[var] < lower_bound)
  dat[outliers, var] = NA 
  residuals_df[outliers, var] = NA
  
  # Perform Kolmogorov-Smirnov Test for Normality and test for significance with Bonferroni Correction 
  resids = na.omit(residuals_df[,var])
  ks_test <- ks.test(resids, "pnorm", mean=mean(resids), sd = sd(resids)) 
  normality_pvals <- c(normality_pvals, ks_test$p.value)
}

# Plot residuals: Histogram
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



# Plot residuals: Q-Q Plot
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

CairoPNG(plot_file, width = 10, height = 10, units = "in", res = 100)
print(residuals_plots)
dev.off()

complete_dat <- dat %>% filter(!if_all(all_of(outcome_vars), is.na))
#write.csv(complete_dat, file = outlier_removed_data_file, row.names = FALSE)

