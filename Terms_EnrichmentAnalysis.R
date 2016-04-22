# This script performs DO/GO enrichment analysis
rm(list=ls())
setwd("ENTER PATH")
cat("\014")

library(parallel)
library(GOSemSim)
#library(data.table)

source("CoefAnalysis.R")
source("EnrichmentAnalysis.R")

Data_path <- "ENTER PATH"
Result_path <- "ENTER PATH"
Figure_path <- "ENTER PATH"


getALLEG <- function(organism) {
  annoDb <- getDb(organism)
  require(annoDb, character.only = TRUE)
  annoDb <- eval(parse(text=annoDb))
  eg=keys(annoDb, keytype="ENTREZID")
  return(eg)
}

background_genes <- getALLEG("human")



GO_common <- unlist(read.csv(paste(Result_path,"GO_CommonTerms.csv",sep="")))
length(GO_common)
DO_common <- unlist(read.csv(paste(Result_path,"DO_CommonTerms.csv",sep="")))
length(DO_common)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         DO and GO term enrichment analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("DOGO_Enrichment.R")

source("gene2DO.R")
source("gene2GO.R")

n_gene <- 13
n_rep <- 100

DO_pval <- mcapply(DO_common, function(x) DOGO_Enrichment(real_ID=x, n_gene, n_rep, background_genes))
GO_pval <- mclapply(GO_common, function(x) DOGO_Enrichment(real_ID=x, n_gene, n_rep, background_genes))

GO_info <- cbind(GOID=as.character(GO_common),GO_pval=unlist(GO_pval))
DO_info <- cbind(DOID=as.character(DO_common),DO_pval=unlist(DO_pval))

write.csv(GO_info, file=paste(Result_path,"GO_pval.csv",sep=""), row.names=F)
write.csv(DO_info, file=paste(Result_path,"DO_pval.csv",sep=""), row.names=F)
