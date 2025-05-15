
# function to compute a linear combination of effects from regression
# compute value, SE, and test the combination

source("2024_APOE_ATN_biomarkers/Code/format_pvalue.r")


# suppose we fit a model, called mod

test_coef_comb <- function(mod, named_vector_vals){
  betas <- summary(mod)$coef[,"Estimate"]
  sigma_betas <- vcov(mod)
  
  stopifnot(all(names(named_vector_vals) %in% names(betas)))
  
  combined_est <- sum(named_vector_vals*betas[names(named_vector_vals)])
  combined_est_variance <- t(named_vector_vals) %*% 
    sigma_betas[names(named_vector_vals),names(named_vector_vals)] %*% 
    named_vector_vals
  
  
  combined_est_se <- sqrt(combined_est_variance)
  
  chi_sq_test_stat <- combined_est^2/combined_est_variance
  
  pval <- pchisq(chi_sq_test_stat, df =1, lower.tail = FALSE)
  
  return(c(combined_est = combined_est, 
           combined_est_se = combined_est_se, 
           chi_sq_test_stat = chi_sq_test_stat,
           pval = pval, 
           paper_pval = format_pvalue(pval)))
  
}
