EnrichmentAnalysis <- function(coef, n_top, target, valTarget, nboot, seed)
{
  #browser()
  set.seed(seed)
  if (is.list(coef)) 
  {
    coef_num <- length(coef) 
  }else{ coef_num <- 1}
  
  p_value <- rep(NA,coef_num)
  
  putative_tar_ind <- which(target==1,arr.ind=T)
  
  i=1
  for (i in 1:coef_num)
  {
    
    if (coef_num > 1)
    {
      c <- coef[[i]] 
    }else{ c <- coef}
    
    top <- n_top[i]
    
    # Find indices of non-zero coefs and sort them based on their absolute value
    coef_ind <- which(c!=0,arr.ind = T)
    coef_ind <- cbind(coef_ind,c[coef_ind])
    coef_ind <- coef_ind[order(-abs(coef_ind[,3])),] # Sort in descending order
    
    # Making 0-1 matrix which '1' show top-ranked interactions
    top_coef <- c *0
    top_coef[coef_ind[1:top,1:2]] <- 1
    
    # Real number for identified validate interactions using algorithm
    identified <- sum(valTarget & top_coef)
    
    num <- rep(0,nboot)
    for (k in 1:nboot)
    {
      # sampling n_top putative interactions among non-zero indices of target
      samp_ind <- sample(1:nrow(putative_tar_ind),top)
      
      rand_coef <- target * 0
      rand_coef[putative_tar_ind[samp_ind,1:2]] <- 1
      
      # Number of validated targets among n_top random selected putative interactions
      num[k] <- sum(valTarget & rand_coef)  
    }
    
    p_value[i] <- sum(num >= identified)/nboot
  }
  return(p_value)
}