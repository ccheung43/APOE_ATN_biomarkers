# Heterogeneity testing for effect modification of genetic background group on APOE allele 
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(survey)
source("2024_APOE_ATN_biomarkers/Code/cochran.Q.cor.1.R")
library(CompQuadForm)

# Read in data from csv file 
heterogeneity_file = "2024_APOE_ATN_biomarkers/Data/20241118_gengroup_heterogeneity.csv"
effect_modification_file = "2024_APOE_ATN_biomarkers/Data/20241118_gengroup_effect_modification.csv"
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file) %>% 
  mutate(GENGROUP = factor(GENGROUP)) %>% 
  filter(!is.na(GENGROUP))

# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)
# Heterogeneity Testing 
heterogeneity_df <- data.frame(Biomarker = character(0), Heterogeneity_pval = numeric(0))
# Survey Model with Genetic Analysis Group
modelFormula <- "~ -1 + GENGROUP + E2_add + E4_add + AGE + SEX + CENTER + EV1 + EV2 + EV3 + EV4 + EV5"
outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
for (var in outcome_vars){
  mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj, family = gaussian())
  
  model_betas <- coef(mod)
  sigma_betas <- vcov(mod)
  
  # Extract coefficients
  inds <- grep("GENGROUP", names(model_betas))
  gengrp_betas <- model_betas[inds]
  gengrp_sigma <- sigma_betas[inds, inds]
  
  # Extract vars and covs
  gengrp_variances <- diag(gengrp_sigma)
  gengrp_covariances <- gengrp_sigma[lower.tri(gengrp_sigma)]
  
  # Run Cochran's Q test
  Qtest_pval <- cochran.Q.cor.1(gengrp_betas, gengrp_variances, gengrp_covariances,"MetaGLS")
  heterogeneity_df <- rbind(heterogeneity_df, data.frame(Biomarker = var, Heterogeneity_pval = Qtest_pval))
}

# write.csv(heterogeneity_df, heterogeneity_file, row.names = FALSE)


# Effect Modification Testing 
EM <- vector(mode = "list")
# Survey Model with Genetic Analysis Group
modelFormula <- "~  GENGROUP*(E2_add + E4_add) + AGE + SEX + CENTER"
outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
for (var in outcome_vars){
  mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj, family = gaussian())
  
  model_betas <- coef(mod)
  model_pvals <- summary(mod)$coefficients[, "Pr(>|t|)"]
  
  # Extract coefficients
  interaction_term_inds <- grep(":", names(model_betas))
  interaction_terms <- names(model_betas)[interaction_term_inds]
  interaction_term_betas <- model_betas[interaction_term_inds]
  interaction_term_pvals <- model_pvals[interaction_term_inds]
  
  EM[[var]] <- data.frame(interaction = interaction_terms, est = interaction_term_betas, pval = interaction_term_pvals) %>% 
    mutate(biomarker = var) %>% 
    select(biomarker, interaction, est, pval)
}

EM <- do.call(rbind, EM)
rownames(EM) <- NULL
#write.csv(EM, effect_modification_file, row.names = FALSE)
