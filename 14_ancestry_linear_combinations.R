# Creating df for Table 4: Associations between ATN biomarkers and APOE allele by genetic ancestry group
# by testing linear combination of regression coefficients for proportion genetic ancestry
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(survey)
library(openxlsx)
library(readxl)
source("2024_APOE_ATN_biomarkers/Code/test_coef_comb.r")

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file) %>%
  mutate(AFRICAN = MEAN_AFR, EUROPEAN = MEAN_EUR, AMERINDIAN = MEAN_AMER)
ancestry_associations_file = "2024_APOE_ATN_biomarkers/Data/20250124_genetic_ancestry_linear_combinations.xlsx"


# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)

APOE_ATN_interaction_by_ANCESTRY <- function(E2_exposure, E4_exposure, covarsFormula) { 
  results_E2 <- c()
  results_E4 <- c()
  for (ancestry_group in genetic_ancestry_groups) { 
    print(ancestry_group)
    modelFormula <- paste0("~", E2_exposure, "+", E4_exposure, "+", ancestry_group, "+", 
                           E2_exposure, "*", ancestry_group, "+", E4_exposure, "*", ancestry_group, "+", covarsFormula) 
    
    for (var in outcome_vars){
      print(var)
      mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj, family = gaussian())
      p_intervals <- seq(0, 1, by = 0.2)
      for (p in p_intervals) { 
        weights <- c()
        weights[E2_exposure] = 1
        weights[paste0(E2_exposure, ":", ancestry_group)] = p
        results_E2 <- rbind(results_E2, c(allele = "E2", ancestry = ancestry_group, biomarker = var, proportion = p, test_coef_comb(mod, weights))) 
        weights <- c()
        weights[E4_exposure] = 1
        weights[paste0(E4_exposure, ":", ancestry_group)] = p
        results_E4 <- rbind(results_E4, c(allele = "E4", ancestry = ancestry_group, biomarker = var, proportion = p, test_coef_comb(mod, weights))) 
      } 
    } 
  } 
  results <- rbind(results_E2, results_E4) 
  return(results)
} 

genetic_ancestry_groups <- c("AFRICAN", "EUROPEAN", "AMERINDIAN") 
outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")

# Additive inheritance mode (partially adjusted)
E2_exposure = "E2_add"
E4_exposure = "E4_add"
covarsFormula = "SEX+AGE+CENTER+EV1+EV2+EV3+EV4+EV5"
ancestry_results <- APOE_ATN_interaction_by_ANCESTRY(E2_exposure, E4_exposure, covarsFormula)

format_est_ci <- function(est, se) {
  est = as.numeric(est)
  se = as.numeric(se)
  # compute upper and lower CI bounds
  CI_lower = signif(est - 2*se, 3) 
  CI_upper = signif(est + 2*se, 3) 
  # Combine into single string
  return(sprintf("%s[%s,%s]", signif(est,3) , CI_lower, CI_upper))
}
ancestry_results <- as.data.frame(ancestry_results) %>% mutate(est = format_est_ci(combined_est, combined_est_se)) 

# Write to xlsx sheet with page for each ancestry group
write.xlsx(ancestry_results, file = ancestry_associations_file, overwrite = TRUE) 
