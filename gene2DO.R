# this is a hidden function from DOSE package
gene2DO <- function(gene) {
library(DOSE)
  if(!exists("DOSEEnv"))
  {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)
    
    tryCatch(utils::data(list="DOSEEnv", package="DOSE"))
  }
  EG2DO <- get("EG2DO", envir=DOSEEnv)
  DO <- EG2DO[[gene]]
  DO <- unlist(DO)
  if (is.null(DO)) {
    return(NA)
  }
  if (sum(!is.na(DO)) == 0) {
    return(NA)
  }
  DO <- DO[!is.na(DO)]
  if (length(DO) == 0) {
    return(NA)
  }
  return(DO)
}
