# This script reads output of GenmiR and TaLasso and perform analysis
rm(list=ls())
setwd("ENTER PATH")
cat("\014")

library(R.matlab)
library(igraph)
library(ggplot2)
library(data.table)

source("CoefAnalysis.R")
source("EnrichmentAnalysis.R")

Data_path <- "ENTER PATH"
Result_path <- "ENTER PATH"
Figure_path <- "ENTER PATH"

Type <- "TGCT" #KIRC, TGCT, Prostate
#+++++++++++++++++++++++++++++++++++++++++++++++++
#           Loading gene names and their entrez IDs
#+++++++++++++++++++++++++++++++++++++++++++++++++
# converting entrez ids to  gene symbols
load(paste(Data_path,"geneSymb_Entrez.RData",sep=""))
colnames(names) <- c("Symbol","Entrez")
names_DT <- data.frame(Symbol=names[,1], Entrez=names[,2])

Symbol2Entrez <- function(symbols)
{
  return(as.character(sapply(symbols, function(x) names_DT[names_DT$Symbol %in% x,2])))
}

Entrez2Symbol <- function(entrez)
{
  return(as.character(sapply(entrez, function(x) names_DT[names_DT$Entrez %in% x,1])))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#           Load data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load(paste(Data_path,Type,"_Files.RData",sep=""))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Load returned results by MATLAB
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

genmir <- readMat(paste(Result_path, Type, "_GenMir.mat",sep=""))
if (Type=='TGCT') genmir <- genmir$GenMir
if (Type=='Prostate') genmir <- genmir$GenMiR

dim(genmir)

talasso <- readMat(paste(Result_path,"100_global_M_5_", Type, ".mat",sep=""))
pval <- talasso$pval
talasso <- talasso$solution
dim(talasso)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Check number of validated interactions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++
# which(genmir>0)

# sum((talasso > 0) & valTarget)
# sum((talasso > 0) & genmir)

t_min <- 100
t_max <- 1000
step <- 100
n <- length(seq(t_min, t_max, by=step))

ta_valid <- sapply(seq(t_min, t_max, by=step), function(x) CoefAnalysis(talasso, n_top=x, valTarget))
gen_valid <- sapply(seq(t_min, t_max, by=step), function(x) CoefAnalysis(genmir, n_top=x, valTarget))

ta_valid/seq(t_min, t_max, by=step)
gen_valid/seq(t_min, t_max, by=step)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Plot number of validated interactions against top predicted interactions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DT <- data.frame(N_top=rep(seq(t_min, t_max, by=step),2), 
                 method=c(rep("TaLasso",n), rep("GenMiR++",n)),N_Valid=c(ta_valid, gen_valid))

pdf(paste(Figure_path, "NumValid.pdf",sep=""))
ggplot(DT, aes(x=N_top, y=N_Valid, color=method), group=method) +
  geom_line() + geom_point(size=3)+
  theme(legend.position="top", legend.text=element_text(size=13), text = element_text(size=17), aspect.ratio=1)+ #
  ylab("Number of validated interactions")+
  xlab("Number of top predicted interactions")+
  ggtitle("Number of validated interactions among number of top predicted interactions")
dev.off()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#         Enrichment Analysis
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#step <- 200
nboot <- 100
ta_enrich <- sapply(seq(t_min, t_max, by=step), function(x) EnrichmentAnalysis (talasso, n_top=x, target, valTarget, nboot, seed=123))
genmir_enrich <- sapply(seq(t_min, t_max, by=step), function(x) EnrichmentAnalysis (genmir, n_top=x, target, valTarget, nboot, seed=123))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#               Find common interactions between TaLasso and GenMiR
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Retrieve_top_Genes_miRNAs <- function(coef, target, top)
{
  coef_ind <- which(coef!=0, arr.ind = T)
  coef_ind <- cbind(coef_ind, coef[coef_ind])
  coef_ind <- coef_ind[order(-abs(coef_ind[,3])),] # Sort in descending order
  head(coef_ind)
  coef_ind <- coef_ind[1:top, 1:2]
  interactions <- cbind(genes=rownames(target)[coef_ind[,1]], mirnas=colnames(target)[coef_ind[,2]])
  return(interactions)
}

top <- seq(t_min, t_max, by=step)

common_interactions <- list()
num_common <- NULL
for (i in 1: length(top))
{
  ta_interactions <- Retrieve_top_Genes_miRNAs(talasso, target, top[i])
  genmir_interactions <- Retrieve_top_Genes_miRNAs(genmir, target, top[i])
  
  # Common interactions between Talasso and genmir
  interactions <- rbind(ta_interactions, genmir_interactions)
  common_interactions[[i]] <- interactions[which(duplicated(interactions)),]
  num_common[i] <- nrow(common_interactions[[i]])  
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# find which of interactions are validated
#++++++++++++++++++++++++++++++++++++++++++++++++
valid <- NULL
for (j in 1:length(common_interactions))
{
  val <- NULL
  for(i in 1:nrow(common_interactions[[j]]))
  {
    r_ind <- which(rownames(valTarget)==common_interactions[[j]][i,1])
    c_ind <- which(colnames(valTarget)==common_interactions[[j]][i,2])
    val[i] <- valTarget[r_ind, c_ind]
  }
  common_interactions[[j]] <- cbind(common_interactions[[j]], val)
  colnames(common_interactions[[j]]) <- c("Entrez ID", "miRNA", "Validated")
  valid[j] <- length(which(common_interactions[[j]][,3]=="1"))
}

##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Saving top-ranked 500 interactions in TaLasso and GenMiR++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ta_interactions <- Retrieve_top_Genes_miRNAs(talasso, target, 500)
genmir_interactions <- Retrieve_top_Genes_miRNAs(genmir, target, 500)

# find which of interactions are validated

valid_ta <- NULL
valid_genmir <- NULL

for(i in 1:nrow(ta_interactions))
{
    r_ind <- which(rownames(valTarget)==ta_interactions[i,1])
    c_ind <- which(colnames(valTarget)==ta_interactions[i,2])
    valid_ta[i] <- valTarget[r_ind, c_ind]
    
    r_ind <- which(rownames(valTarget)==genmir_interactions[i,1])
    c_ind <- which(colnames(valTarget)==genmir_interactions[i,2])
    valid_genmir[i] <- valTarget[r_ind, c_ind]    
}
ta_interactions <- cbind(ta_interactions, valid_ta)
genmir_interactions <- cbind(genmir_interactions, valid_genmir)

ta_interactions[,1] <- Entrez2Symbol(ta_interactions[,1])
genmir_interactions[,1] <- Entrez2Symbol(genmir_interactions[,1])

colnames(ta_interactions) <- c("Gene Symbol", "miRNA", "Validated")
colnames(genmir_interactions) <- c("Gene Symbol", "miRNA", "Validated")

write.csv(ta_interactions, file=paste(Result_path,Type, "_TaLasso_TopRanked 500 interactions.csv", sep=""), row.names=F)
write.csv(genmir_interactions, file=paste(Result_path,Type, "_GenMiR++_TopRanked 500 interactions.csv", sep=""), row.names=F)






DT3 <- data.frame(N_top=rep(top,2), 
                  Type=c(rep("Common",length(top)), rep("Validated", length(top))), 
                  Num= c(num_common,valid))

pdf(paste(Figure_path, "CommonInteraction.pdf",sep=""))
ggplot(DT3, aes(x=N_top, y=Num, color=Type), group=Type) +
  geom_line() + geom_point(size=3)+
  theme(legend.title=element_blank(), legend.position="top", legend.text=element_text(size=13), text = element_text(size=17), aspect.ratio=1)+ 
  ylab("Number of common/validated interactions")+
  xlab("Number of top predicted interactions")+
  ggtitle("Common interactions between TaLasso and GenMiR")
dev.off()

# Select the first 500 top interactions
common_interactions <- common_interactions[[5]] 

# removing validated interactions
common_interactions <- common_interactions[-which(common_interactions[,3]=="1"),]

table(common_interactions[,2])
table(common_interactions[,1])

#Write edgelist for plotting network with CytoScape
EdgeList <- common_interactions[,1:2]
EdgeList[,1] <- Entrez2Symbol(EdgeList[,1])
write.csv(EdgeList, file=paste(Result_path, Type, "_Common_interactions.csv",sep=""), row.names=F)


mi_list <- unique(common_interactions[,2])
ge_list <- unique(common_interactions[,1])

ge_list_symb <- Entrez2Symbol(ge_list)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#related genes to TGCT from: http://www.malacards.org/card/testicular_germ_cell_tumor#related_genes
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
genecards_tgct <- c("BCL10", "KIT", "STK11", "FGFR3", "STK10")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++
#related genes to TGCT from: http://dga.nubic.northwestern.edu/pages/result.php
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
DGA_tgct <- c("ATF7IP", "DMRT1", "ESR1", "ESR2", "KITLG", "LHCGR", "POLG", "TERT")

tgct_genes_symb <- c(genecards_tgct, DGA_tgct)
tgct_genes_entz <- Symbol2Entrez(tgct_genes_symb)

# all_genes <- cbind(Query_genes=ge_list_symb, TGCT_Genes=tgct_genes_symb)
# write.csv(all_genes, file=paste(Data_path, "Genes_of_interest.csv",sep=""), row.names=F)

# All of genes
#Gene_table <- data.table(Symbol=c(ge_list_symb, tgct_genes_symb), Entrez=c(ge_list, tgct_genes_entz))

#+++++++++++++++++++++++++++++++++++++++++++++++++++
#         Plot DO (Disease Ontology) tems similarity matrix
#+++++++++++++++++++++++++++++++++++++++++++++++++++

library(DOSE)
library(corrplot)

DO <- DOSE::geneSim(tgct_genes_entz, ge_list, measure="Wang", combine="BMA")
dim(DO)
#View(DO)
DO <- DO[-5,-c(4,8)] # removing NA rows and columns

colnames(DO) <- Entrez2Symbol(colnames(DO))
rownames(DO) <- Entrez2Symbol(rownames(DO))

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", 
                           "#007FFF", "blue", "#00007F"))

pdf(paste(Figure_path, "DOTermSimilarity.pdf",sep=""))
corrplot(DO, method = "circle", col=col1(10), cl.lim = c(0, 1), 
         cl.ratio = 0.3, tl.col = "black", tl.cex=1.1, title="DO Term Similarity")
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++
#         Plot GO tems similarity matrix
#+++++++++++++++++++++++++++++++++++++++++++++++++++
library(GOSemSim)
library(corrplot)

GO <- GOSemSim::mgeneSim(genes=c(ge_list, tgct_genes_entz),
         ont="MF", organism="human", measure="Wang",
         verbose=TRUE)

GO <- GO[is.element(rownames(GO), tgct_genes_entz), is.element(colnames(GO), ge_list)]

colnames(GO) <- Entrez2Symbol(colnames(GO))
rownames(GO) <- Entrez2Symbol(rownames(GO))
# View(GO)
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", 
                           "#007FFF", "blue", "#00007F"))

pdf(paste(Figure_path, "GOTermSimilarity.pdf",sep=""))
corrplot(GO, method = "circle", col=col1(10), cl.lim = c(0, 1), 
         cl.ratio = 0.3, tl.col = "black", tl.cex=1.1, title="GO Term Similarity")
dev.off()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#        finding common GO Terms between genes and TGCT genes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("gene2GO.R")

GO_1 <- lapply(tgct_genes_entz, gene2GO)
GO_tgct <- NULL
i=1
for (i in 1:length(GO_1))
{
  GO <- GO_1[[i]]
  if (all(!is.na(GO)))
    GO_tgct <- rbind(GO_tgct, cbind(rep(tgct_genes_symb[i], length(GO)),GO))
}

## RPS17: GO:0044822" "GO:0003735

GO_2 <- lapply(ge_list, gene2GO)
GO_list <- NULL
for (i in 1:length(GO_2))
{
  GO <- GO_2[[i]]
  if (all(!is.na(GO)))
    GO_list <- rbind(GO_list, cbind(rep(ge_list_symb[i], length(GO)),GO))
}

GO_all <- rbind(GO_list, GO_tgct)
nrow(GO_all)

write.csv( intersect(GO_list[,2], GO_tgct[,2]), file=paste(Result_path,"GO_CommonTerms.csv",sep=""),row.names=F)
ind <- which(is.element(GO_all[,2], intersect(GO_list[,2], GO_tgct[,2])))
nrow(GO_all)
GO_all <- GO_all[ind,]
unique(GO_all[-ind,1])
length(unique(GO_list[,1]))

#Write edgelist for plotting network with CytoScape
write.csv(GO_all, file=paste(Result_path,"GO.csv",sep=""),row.names=F)

#++++++++++++++++++++++++++++++++++++++++++++++++
#           Finding common DO (Disease Ontology) Terms
#++++++++++++++++++++++++++++++++++++++++++++++++

source("gene2DO.R")

DO_1 <- lapply(tgct_genes_entz, gene2DO)
DO_tgct <- NULL
for (i in 1:length(DO_1))
{
  DO <- DO_1[[i]]
  if (all(!is.na(DO)))
    DO_tgct <- rbind(DO_tgct, cbind(rep(tgct_genes_symb[i], length(DO)),DO))
}


DO_2 <- lapply(ge_list, gene2DO)
DO_list <- NULL
i=1
for (i in 1:length(DO_2))
{
  DO <- DO_2[[i]]
  if (all(!is.na(DO)))
    DO_list <- rbind(DO_list, cbind(rep(ge_list_symb[i], length(DO)),DO))
}

DO_all <- rbind(DO_list, DO_tgct)
DO_all <- DO_all[!duplicated(DO_all),]

write.csv( intersect(DO_list[,2], DO_tgct[,2]), file=paste(Result_path,"DO_CommonTerms.csv",sep=""),row.names=F)

ind <- which(is.element(DO_all[,2], intersect(DO_list[,2], DO_tgct[,2])))
nrow(DO_all)
DO_all <- DO_all[ind,]
#Write edgelist for plotting network with CytoScape
write.csv(DO_all, file=paste(Result_path,"DO.csv",sep=""),row.names=F)

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#             Analysis of GeneMANIA networks
#++++++++++++++++++++++++++++++++++++++++++++++++++++
geneMANIA <- read.table(paste(Data_path,"GeneMANIA_interactions_list.txt", sep=""), header=F,
                        check.names=F, fill=T, stringsAsFactors = F)
dim(geneMANIA)
geneMANIA <- geneMANIA[-1,]
intr <- geneMANIA[,c(1,2,4)]
colnames(intr) <- c("Gene1", "Gene2", "Type")
geneMANIA_DT <- data.table(intr)

unique(geneMANIA_DT[,Type])
coexpression <- geneMANIA_DT[Type=="Co-expression", .(Gene1, Gene2)]
nrow(coexpression)

write.csv(coexpression, paste(Result_path, "GeneMANIA_Coexpression.csv", sep=""), row.names=F)

physical <- geneMANIA_DT[Type=="Physical&nbsp;interactions", .(Gene1, Gene2)]

#For cytoscape
write.csv(physical, paste(Result_path, "GeneMANIA_Physical.csv", sep=""), row.names=F)

