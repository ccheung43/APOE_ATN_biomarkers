# code from Gendered_indices_and_insomnia
# function to format p-value with stars following a reviewer's request
# gets a vector of p-values named pval



format_pvalue <- function(pval, max_digits = 3, stars = c(0.05, 0.01, 0.001)){
  if (length(stars) > 3) stop("only allow three stars")
  stopifnot(all(pval>= 0 ))
  
  pval_round <- formatC(pval, digits = max_digits, format = "f")
  
  stars_vec <- rep("", length = length(pval))
  
  inds_3 <- which(pval < stars[3])
  if (length(inds_3) > 0){
    stars_vec[inds_3] <- c("***")
    pval_round[inds_3] <- paste0("<", stars[3]) 
  }
  
  inds <- which(pval < stars[2] & pval>= stars[3]) 
  if (length(inds) > 0){
    stars_vec[inds] <- c("**")
  }
  
  inds <- which(pval < stars[1] & pval>= stars[2]) 
  if (length(inds) > 0){
    stars_vec[inds] <- c("*")
  }
  
  pval_vec <- paste0(pval_round, stars_vec)
  
  return(pval_vec)  
}