# Over-representation Analysis

# the hypergeometric distribution is traditionally explained in terms of
# drawing a sample of balls from an urn containing black and white balls.
# to keep the arguments straight (in my mind at least), I use these terms
# here also

rm(list=ls())
setwd("C:/Users/sedaghat/Dropbox/TGCT_miRNA-mRNA Network Analysis/Codes")
cat("\014")

set.seed(123)
library(KEGGREST)
library(org.Hs.eg.db)
library(ggplot2)

Data_path <- "C:/Users/sedaghat/Dropbox/TGCT_miRNA-mRNA Network Analysis/Data/"
Result_path <- "C:/Users/sedaghat/Dropbox/TGCT_miRNA-mRNA Network Analysis/Results/"
Figure_path <- "C:/Users/sedaghat/Dropbox/TGCT_miRNA-mRNA Network Analysis/Docs/Figures/"

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

genes <- read.csv(paste(Data_path, "Genes_of_interest.csv",sep=""))
dim(genes)
query_genes <- Symbol2Entrez(genes[,1])
tgct_genes <- Symbol2Entrez(genes[,2])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#             Read pathways
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# created named list, eg:  path:map00010: "Glycolysis / Gluconeogenesis" 

pathways.list <- keggList("pathway", "hsa")

# make them into KEGG-style human pathway identifiers
pathway.codes <- sub("path:", "", names(pathways.list))

# for demonstration, just use the first ten pathways
# not all pathways exist for human, so TODO: tryCatch the
# keggGet to be robust against those failures

# subsetting by c(TRUE, FALSE) -- which repeats
# as many times as needed, sorts through some
# unexpected packaging of geneIDs in the GENE element
# of each pw[[n]]

genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             pw[[1]]$GENE[c(TRUE, FALSE)]
                           })
#load(paste(Data_path, "KEGG_Gene_Pathways.RData", sep=""))

all.geneIDs <- keys(org.Hs.eg.db)

# chose one of these for demonstration.  the first (a whole genome random
# set of 100 genes)  has very little enrichment, the second, a random set
# from the pathways themsevles,  has very good enrichment

genes.of.interest <- tgct_genes

tgct_pvals <- sapply(names(genes.by.pathway),
                     function(pathway) {
                       pathway.genes <- genes.by.pathway[[pathway]]
                       white.balls.drawn <- length(intersect(genes.of.interest, pathway.genes))
                       white.balls.in.urn <- length(pathway.genes)
                       total.balls.in.urn <- length(all.geneIDs)
                       black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
                       total.balls.drawn.from.urn <- length(genes.of.interest)
                       if (white.balls.drawn > 0 ) # if genes of interest hit the pathway
                         1-sum(dhyper(0:white.balls.drawn,
                                      white.balls.in.urn,
                                      black.balls.in.urn,
                                      total.balls.drawn.from.urn)) else NA
                     })


genes.of.interest <- query_genes

query_pvals <- sapply(names(genes.by.pathway),
                      function(pathway) {
                        pathway.genes <- genes.by.pathway[[pathway]]
                        white.balls.drawn <- length(intersect(genes.of.interest, pathway.genes))
                        white.balls.in.urn <- length(pathway.genes)
                        total.balls.in.urn <- length(all.geneIDs)
                        black.balls.in.urn <- total.balls.in.urn - white.balls.in.urn
                        total.balls.drawn.from.urn <- length(genes.of.interest)
                        if (white.balls.drawn > 0 ) # if genes of interest hit the pathway
                        1-sum(dhyper(0:white.balls.drawn,
                                     white.balls.in.urn,
                                     black.balls.in.urn,
                                     total.balls.drawn.from.urn)) else NA
                      })

all_path <- cbind(Pathway_name=pathways.list, Pathway_code=pathway.codes, 
                  Querygenes_pValue= query_pvals, TGCTgenes_pValue= tgct_pvals)
#write.csv(all_path, file=paste(Result_path, "Pathway_PValues.csv", sep=""), row.names=F)

th <- 0.01

tgct_pval_shared <- tgct_pvals[intersect(which(tgct_pvals < th), which(query_pvals < th))]
query_pval_shared <- query_pvals[intersect(which(tgct_pvals < th), which(query_pvals < th))]

length(which(query_pvals < th))
length(which(tgct_pvals < th))
length(query_pval_shared)

Path.names <- pathways.list[paste("path:", names(query_pval_shared),sep="")]

query_hit <- sapply(names(query_pval_shared), function(pathway){ 
                                              pathway.genes <- genes.by.pathway[[pathway]]
                                              Entrez2Symbol(intersect(query_genes, pathway.genes))})

tgct_hit <- sapply(names(query_pval_shared), function(pathway){ 
                                              pathway.genes <- genes.by.pathway[[pathway]]
                                              Entrez2Symbol(intersect(tgct_genes, pathway.genes))})

query_nhit <- sapply(query_hit, length)
tgct_nhit <- sapply(tgct_hit, length)

shared_info <- cbind(Pathway.Name =Path.names,
                     Path.ID=names(tgct_pval_shared), 
                     TGCT_hit= tgct_nhit, TGCT_pvalue=tgct_pval_shared, 
                     Quey_hit = query_nhit, Query_pvalues=query_pval_shared )

tgct_info <- cbind(Pathway.Name =Path.names,
                     Path.ID=names(tgct_pval_shared), 
                     hit= tgct_nhit, pvalue=tgct_pval_shared, 
                     Type=rep("TGCT", length(Path.names)))

query_info <- cbind(Pathway.Name =Path.names,
                     Path.ID=names(tgct_pval_shared), 
                     hit = query_nhit, pvalues=query_pval_shared,
                    Type=rep("Query", length(Path.names)))
#++++++++++++++++++++++++++++++++++++++
#           Plotting
#++++++++++++++++++++++++++++++++++++++

cbPalette <- c("dodgerblue", "green3")

info <- rbind(tgct_info, query_info)
rownames(info) <- NULL
DT <- data.frame(info)

pdf(paste(Figure_path, "PathwayBarPlot.pdf",sep=""))

DT$pvalue <- as.numeric(as.character(DT$pvalue))
options(digits=3)

ggplot(DT, aes(x = factor(Path.ID), y = hit, fill=factor(Type), width=0.4)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=cbPalette , name="")+
  geom_text(aes(label=paste(format(DT$pvalue, digits=3, scientific=T))), position=position_dodge(width=1), size=5, angle=90) +
  scale_x_discrete("") +
  scale_y_discrete("Number of Genes")+
  theme(axis.text.x = element_text(size=14, angle = 90, hjust = 1), axis.text.y = element_text(size=14))


plot(x=0, y=0, axes=F, frame.plot=F, col="white", xaxt='n', ann=F, yaxt='n' )
leg_text <- paste(Path.ID=names(tgct_pval_shared),": ", sub(" - Homo sapiens (human)","",Path.names, fixed = T) , sep="")

leg_text <- sapply (1:length(query_hit), function(x) 
                                          paste(names(query_hit)[x],": ", sub(" - Homo sapiens (human)","", Path.names[x], fixed = T), 
                                          "\n", "                 ", 
                                          "Query genes:", paste(c(query_hit[[x]]), collapse = ', '),
                                          "\n", "                 ", 
                                          "TGCT genes:", paste(c(tgct_hit[[x]]), collapse = ', '), '\n'))

legend("topleft", leg_text, cex=0.9, box.col="white")

dev.off()
x=1

sapply (1:length(query_hit), function(x) paste(names(query_hit)[x], "=>", 
                                               "Query genes:", paste(c(query_hit[[x]]), collapse = ', '),
                                               "\n", "T"))


#########################################################################
#########################################################################

#           Calculating p-value for overlap

#########################################################################
#########################################################################

# length(B) < length(A)
A <- which(tgct_pvals < th)
B <- which(query_pvals < th)
N <- length(pathways.list)
# sum(sapply(length(intersect(A,B)):length(B), function(i) choose(length(A),i)*choose(N-length(A),length(B)-i)/choose(N,length(B))))

(overlap_pval <- sum(dhyper(length(intersect(A,B)):length(B), length(A), N-length(A), length(B))))

#N=10;M=5;K=3;l=2
#(choose(M,l)*choose(N-M,K-l))/choose(N,K)
#dhyper(l,M,N-M,K)
