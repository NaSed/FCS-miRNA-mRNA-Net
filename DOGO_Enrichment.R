DOGO_Enrichment <- function(real_ID, n_gene, n_rep, background_genes)
{
  res <- NULL
  TermType <- substr(real_ID,1,2)
  for (i in 1 : n_rep)
  {
    cat('rep ',i, '\n')
    g_l <- sample(background_genes, n_gene)
    if (TermType=="GO")
      l <- unlist(lapply(g_l, gene2GO))
    else l <- unlist(mclapply(g_l, gene2DO))
    res[i] <- is.element(real_ID, l)
  }
  pval <- mean(res)
  return(pval)
}