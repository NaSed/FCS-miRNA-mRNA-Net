CoefAnalysis <- function(coef, n_top, valTarget)
{
  
  if (is.list(coef)) num <- length(coef) else num <- 1
  
  ident <- rep(0,num)
  
  for (i in 1:num)
  {
    if (num > 1) c <- coef[[i]] else c <- coef
    top <- n_top[i]
    coef_ind <- which(c!=0,arr.ind = T)
    coef_ind <- cbind(coef_ind,c[coef_ind])
    coef_ind <- coef_ind[order(-abs(coef_ind[,3])),] # Sort in descending order
    
    # Making 0-1 matrix which '1' show top-ranked interactions
    top_coef <- c * 0
    top_coef[coef_ind[1:top,1:2]] <- 1
    
    ident[i] <- sum(valTarget & top_coef)
  }
  return(ident)
}