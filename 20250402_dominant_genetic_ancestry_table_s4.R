# Creating df for Table S4: Associations between ATN biomarkers and APOE allele by genetic ancestry group (dominant mode)
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(survey)
source("2024_APOE_ATN_biomarkers/Code/format_pvalue.r")
library(openxlsx)
library(readxl)

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file) %>%
  mutate(AFRICAN = MEAN_AFR, EUROPEAN = MEAN_EUR, AMERINDIAN = MEAN_AMER)
table_s4_file = "2024_APOE_ATN_biomarkers/Results/20250402_table_s4.xlsx"
ancestry_associations_file = "2024_APOE_ATN_biomarkers/Data/20250402_dominant_associations_by_genetic_ancestry.xlsx"

# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)

# from Gendered_indices_and_insomnia 
extract_one_out <- function(mod, exposure, out_name, exponentiate = FALSE, sigfigs = 2){
  est <- summary(mod)$coef
  confint <- confint(mod)
  
  ind <- which(rownames(est) == exposure) # modified to find only exact case 
  
  if (exponentiate){
    df_out <- data.frame(outcome = out_name, 
                         est = signif(exp(est[ind, "Estimate"]),sigfigs), 
                         CI = paste0("(", signif(exp(confint[ind, 1]), sigfigs), ",",
                                     signif(exp(confint[ind, 2]), sigfigs), ")"),
                         pval = formatC(est[ind, "Pr(>|t|)"], digits = sigfigs, format = "E"),
                         paper_pval= format_pvalue(est[ind, "Pr(>|t|)"]))
  } else{
    df_out <- data.frame(outcome = out_name, 
                         est = signif(est[ind, "Estimate"],sigfigs), 
                         CI = paste0("(", signif(confint[ind, 1], sigfigs), ",",
                                     signif(confint[ind, 2], sigfigs), ")"),
                         pval = formatC(est[ind, "Pr(>|t|)"], digits = sigfigs, format = "E"),
                         paper_pval = format_pvalue(est[ind, "Pr(>|t|)"]))
  }
  
  return(df_out)
}


APOE_ATN_interaction_by_ANCESTRY <- function(E2_exposure, E4_exposure, covarsFormula, num_permutations) { 
  all_tables <- list()
  for (ancestry_group in genetic_ancestry_groups) { 
    print(ancestry_group)
    modelFormula <- paste0("~", E2_exposure, "*", ancestry_group, "+", E4_exposure, "*", ancestry_group, "+", covarsFormula) 
    
    E2 <- vector(mode = "list")
    E4 <- vector(mode = "list")
    E2_interaction <- vector(mode = "list")
    E4_interaction <- vector(mode = "list")
    ANCESTRY <- vector(mode = "list")
    for (var in outcome_vars){
      print(var)
      mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj, family = gaussian())
      
      E2[[var]] <- extract_one_out(mod, 
                                   exposure = E2_exposure, 
                                   out_name  = var, 
                                   exponentiate = FALSE)
      E4[[var]] <- extract_one_out(mod, 
                                   exposure = E4_exposure, 
                                   out_name  = var, 
                                   exponentiate = FALSE)
      
      E2_interaction[[var]] <- extract_one_out(mod, 
                                               exposure = paste0(E2_exposure, ":", ancestry_group), 
                                               out_name  = var, 
                                               exponentiate = FALSE)
      E4_interaction[[var]] <- extract_one_out(mod, 
                                               exposure = paste0(ancestry_group, ":", E4_exposure), 
                                               out_name  = var, 
                                               exponentiate = FALSE)
      
      ANCESTRY[[var]] <- extract_one_out(mod, 
                                         exposure = ancestry_group, 
                                         out_name  = var, 
                                         exponentiate = FALSE)
      
      # Permutation testing
      e2_permuted_stats <- numeric(num_permutations)
      e4_permuted_stats <- numeric(num_permutations)
      e2_interaction_permuted_stats <- numeric(num_permutations)
      e4_interaction_permuted_stats <- numeric(num_permutations)
      ancestry_permuted_stats <- numeric(num_permutations)
      
      for (i in 1:num_permutations) {
        permuted_dat <- dat %>%
          mutate(across(all_of(outcome_vars), ~ sample(.))) 
        
        # Fit the model with permuted data
        permuted_svy_obj <- svydesign(id = ~PSU_ID, strata = ~STRAT, weights = ~WEIGHT_NORM_OVERALL_INCA, 
                                      nest = TRUE, data = permuted_dat)
        permuted_mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = permuted_svy_obj, family = gaussian())
        
        # Store permuted estimates for E2 and E4 exposures
        e2_permuted_stats[i] <- summary(permuted_mod)$coef[E2_exposure, "Estimate"]
        e4_permuted_stats[i] <- summary(permuted_mod)$coef[E4_exposure, "Estimate"]
        e2_interaction_permuted_stats[i] <- summary(permuted_mod)$coef[paste0(E2_exposure, ":", ancestry_group), "Estimate"]
        e4_interaction_permuted_stats[i] <- summary(permuted_mod)$coef[paste0(ancestry_group, ":", E4_exposure), "Estimate"]
        ancestry_permuted_stats[i] <- summary(permuted_mod)$coef[ancestry_group, "Estimate"]
      }
      
      # Calculate permutation-based p-values
      if (E2[[var]]$est < 0) {
        E2[[var]]$perm_pval <- mean(e2_permuted_stats <= E2[[var]]$est)
      } else {
        E2[[var]]$perm_pval <- mean(e2_permuted_stats >= E2[[var]]$est)
      }
      
      if (E4[[var]]$est < 0) {
        E4[[var]]$perm_pval <- mean(e4_permuted_stats <= E4[[var]]$est)
      } else {
        E4[[var]]$perm_pval <- mean(e4_permuted_stats >= E4[[var]]$est)
      }
      
      if (E2_interaction[[var]]$est < 0) {
        E2_interaction[[var]]$perm_pval <- mean(e2_interaction_permuted_stats <= E2_interaction[[var]]$est)
      } else {
        E2_interaction[[var]]$perm_pval <- mean(e2_interaction_permuted_stats >= E2_interaction[[var]]$est)
      }
      
      if (E4_interaction[[var]]$est < 0) {
        E4_interaction[[var]]$perm_pval <- mean(e4_interaction_permuted_stats <= E4_interaction[[var]]$est)
      } else {
        E4_interaction[[var]]$perm_pval <- mean(e4_interaction_permuted_stats >= E4_interaction[[var]]$est)
      }
      
      if (ANCESTRY[[var]]$est < 0) {
        ANCESTRY[[var]]$perm_pval <- mean(ancestry_permuted_stats <= ANCESTRY[[var]]$est)
      } else {
        ANCESTRY[[var]]$perm_pval <- mean(ancestry_permuted_stats >= ANCESTRY[[var]]$est)
      }
    } 
    
    E2_tab <- do.call(rbind, E2) %>% mutate(term = "E2")
    E4_tab <- do.call(rbind, E4) %>% mutate(term = "E4")
    E2_interaction_tab <- do.call(rbind, E2_interaction) %>% mutate(term = "E2_interaction")
    E4_interaction_tab <- do.call(rbind, E4_interaction) %>% mutate(term = "E4_interaction")
    ANCESTRY_tab <- do.call(rbind, ANCESTRY) %>% mutate(term = ancestry_group)
    
    ancestry_table <- rbind(E2_tab, E2_interaction_tab, E4_tab, E4_interaction_tab, ANCESTRY_tab) %>% 
      mutate(biomarker = factor(outcome, levels = c("ABETA4240", "PTAU181","NFLIGHT", "GFAP"))) %>% 
      arrange(biomarker) %>% 
      select(biomarker, term, est, CI, perm_pval) %>% 
      mutate(paper_pval = format_pvalue(perm_pval))
    
    all_tables[[ancestry_group]] <- ancestry_table
  } 
  return(all_tables)
} 

genetic_ancestry_groups <- c("AFRICAN", "EUROPEAN", "AMERINDIAN") 
outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")

# Additive inheritance mode (partially adjusted)
e2_exposure = "E2_dom"
e4_exposure = "E4_dom"
covarsFormula = "SEX+AGE+CENTER+EV1+EV2+EV3+EV4+EV5"
tbl_by_ancestry <- APOE_ATN_interaction_by_ANCESTRY(e2_exposure, e4_exposure, covarsFormula, 10000)

# Write to xlsx sheet with page for each ancestry group
write.xlsx(tbl_by_ancestry, file = ancestry_associations_file, overwrite = TRUE) 


# Create table S4

# write function to format table2 estimate as est[CI]
format_est_ci <- function(est, CI) {
  # Extract lower and upper bounds from CI
  bounds <- as.numeric(unlist(strsplit(gsub("[()]", "", CI), ",")))
  # Combine into single string
  return(sprintf("%s[%s,%s]", est, bounds[1], bounds[2]))
}

dat_file <- ancestry_associations_file
sheet_names <- excel_sheets(dat_file)
associations_tables <- lapply(sheet_names, function(sheet) read_excel(dat_file, sheet = sheet))
names(associations_tables) <- sheet_names 

tbls4 <- data.frame("ATN biomarker" = c("Aꞵ42/40", "", "pTau-181", "", "NfL", "", "GFAP", ""), 
                   "APOE allele" = c("ε2", "ε4", "ε2", "ε4", "ε2", "ε4", "ε2", "ε4"))
for(ancestry in genetic_ancestry_groups) { 
  ancestry_df <- associations_tables[[ancestry]] %>% 
    mutate("est[CI]" = mapply(format_est_ci, est, CI)) %>% 
    mutate("p-value" = paper_pval)
  ancestry_allele_df <- ancestry_df %>% 
    filter(term %in% c("E2", "E4")) %>% 
    select(c("est[CI]", "p-value"))
  ancestry_interaction_df <- ancestry_df %>% 
    filter(term %in% c("E2_interaction", "E4_interaction"))%>% 
    select(c("est[CI]", "p-value"))
  tbls4 <- cbind(tbls4, ancestry_interaction_df, ancestry_allele_df)
} 

row1 <- c("", "", "African", "", "", "", "European", "", "", "", "Amerindian", "", "", "")
row2 <- c("", "", "Interaction", "", "APOE allele", "", "Interaction", "", "APOE allele", "", "Interaction", "", "APOE allele", "")
row3 <- c("ATN biomarker", "APOE allele", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value")
tbls4 <- rbind(row1, row2, row3, tbls4)

write.xlsx(tbls4, file = table_s4_file, overwrite = TRUE, colNames=FALSE)  

