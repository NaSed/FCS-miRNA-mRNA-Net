rm(list=ls())
setwd("ENTER PATH")
cat("\014")

Type <- "TGCT" #Prostate"
#+++++++++++++++++++++++++++++++++
#       Load data 
#+++++++++++++++++++++++++++++++++

Data_path <- "ENTER PATH"

load(paste(Data_path,Type, "_miRNA.RData",sep=""))
load(paste(Data_path,Type, "_mRNA.RData",sep=""))
load(paste(Data_path,"predicted.RData",sep=""))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Conversion of colnames of target from ensemble id to entrez id
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(org.Hs.eg.db)
Ens2Enz <- as.list(org.Hs.egENSEMBL2EG)
id <- Ens2Enz[match(colnames(target), names(Ens2Enz))]
id <- id[!duplicated(id)]
ind <- which(sapply(1:length(id), function(x) length(id[[x]]))==1)
id <- id[ind]

comm <- intersect(colnames(target), names(id))
target <- target[, match(comm, colnames(target))]
Entz <- as.numeric(id[match(comm,names(id))])
colnames(target) <- Entz

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Conversion of rownames of mRNA to Entrez ID from Symbol|Entrez
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

split <- sapply(rownames(mRNA), function(x)  which(strsplit(x,'')[[1]]=="|"))
x=1
# names <- lapply(1:length(split), function(x)  c(substr(rownames(mRNA)[x],1,split[x]-1),
#   substr(rownames(mRNA)[x],split[x]+1,nchar(rownames(mRNA)[x]))))
# names <- do.call(rbind, names)
# save(names, file=paste(Data_path,"geneSymb_Entrez.RData",sep=""))
rownames(mRNA) <- sapply(1:length(split), function(x)  substr(rownames(mRNA)[x],split[x]+1,nchar(rownames(mRNA)[x])))
#++++++++++++++++++++++++++++++++++++++++++
# sort patients in both mRNA and miRNA
#++++++++++++++++++++++++++++++++++++++++++

ind <- match(colnames(miRNA),colnames(mRNA))
mRNA <- mRNA[,ind]

#TEST: all(colnames(mRNA)==colnames(miRNA))

#+++++++++++++++++++++++++++++++++++++++++++++++++
# Check type of samples, i.e. normal or tumor
#+++++++++++++++++++++++++++++++++++++++++++++++++

#All of them are tumors
data.frame(table(substr(colnames(mRNA),14,15)))


#+++++++++++++++++++++++++++++++++++++++++++++++++
# Read Validated target (miRWalk 2.0)
#      Make validated target matrix     
#+++++++++++++++++++++++++++++++++++++++++++++++++
load(paste(Data_path,".RData",sep="hsa-vtm-entrez.rdata"))
names(id) <- gsub("mir","miR", names(id))

all_genes <- unique(unlist(id))

valTarget <- matrix(0, nrow=length(id), ncol=length(all_genes))
dim(valTarget)

i=1
for (i in 1: length(id))
{
  entz <- id[[i]]
  ind <- match(entz, all_genes)
  valTarget[i,ind] <- 1
}
rownames(valTarget) <- names(id)
colnames(valTarget) <- all_genes

#+++++++++++++++++++++++++++++++++++++++++++++
#       Intersect among all data
#+++++++++++++++++++++++++++++++++++++++++++++
rownames(miRNA) <- gsub("mir","miR", rownames(miRNA))

miR <- intersect(rownames(target),intersect(rownames(miRNA),rownames(valTarget)))

target <- target[match(miR,rownames(target)),]
valTarget <- valTarget[match(miR,rownames(valTarget)),]
miRNA <- miRNA[match(miR,rownames(miRNA)),]


mR <- intersect(rownames(mRNA),intersect(colnames(target),colnames(valTarget)))

target <- target[,match(mR,colnames(target))]
valTarget <- valTarget[,match(mR,colnames(valTarget))]
mRNA <- mRNA[match(mR,rownames(mRNA)),]

target <- t(target)
valTarget <- t(valTarget)

which(apply(target,2,sum)==0)
ind <- which(apply(target,1,sum)==0)
mRNA <- mRNA[-ind,]
target <- target[-ind,]
valTarget <- valTarget[-ind,]

target <- (target | valTarget)

dim(target)
dim(valTarget)
dim(mRNA)
dim(miRNA)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#       Scale read counts on each sample (columns of data)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
standardization <- function(x){(x-min(x))/(max(x)-min(x))}

mRNA <- apply(mRNA, 2, standardization)
miRNA <- apply(miRNA, 2, standardization)

#++++++++++++++++++++++++++++++++++++++
#       Log transform
#++++++++++++++++++++++++++++++++++++++
log_mRNA <- log2(mRNA + 1)
log_miRNA <- log2(miRNA +1)

dim(log_mRNA)
dim(log_miRNA)

save(miRNA,mRNA,target,valTarget,file=paste(Data_path,Type, "_Files.RData",sep=""))

library(R.matlab)
writeMat(paste(Data_path,Type,".mat",sep=""), mRNA=mRNA, miRNA=miRNA, target=((target+1)-1))
  