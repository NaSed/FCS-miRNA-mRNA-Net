rm(list=ls())
setwd("ENTER PATH")
cat("\014")

library(parallel)
save_path <- "ENTER PATH"
Type <- "Prostate"
#++++++++++++++++++++++++++++++++++++++++++++++++
# Read File_Sample_Map.text to find   	+
# matching between sample name and data files	+
#++++++++++++++++++++++++++++++++++++++++++++++++

if (Type=="TGCT") FileSampleMap_path <- "ENTER PATH"
if (Type=="Prostate") FileSampleMap_path <- "ENTER PATH"

file_name <- "FILE_SAMPLE_MAP.txt"
info <- read.csv(paste(FileSampleMap_path,file_name,sep=""),sep="	")
dim(info)
# info[,1] contains name of files
# info[,2] contains patient barcode

#+++++++++++++++++++++++++++++++++++++++++++++++++
# Find rows which correspond to RNASeqV2 files
#+++++++++++++++++++++++++++++++++++++++++++++++++
file_id <- 'rsem.genes.results'
#genes.results, genes.normalized
ind <- which(info[,1] %in% grep(file_id, info[,1], value=TRUE))
file_barcode <- info[ind,]
file_barcode[8,]

#++++++++++++++++++++++++++++++++++
#	Read RNASeqV2 Files
#++++++++++++++++++++++++++++++++++

if (Type=="TGCT") RNA_path <- "ENTER PATH/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
if (Type=="Prostate") RNA_path <- "ENTER PATH/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"

data <- mclapply(paste(RNA_path,file_barcode[,1],sep=""), function(x) read.table(x,sep="	", stringsAsFactors=F))

if (file_id == 'rsem.genes.results')
{
  # removing columns contain scaled_estimate and transcript id in case of raw data
  data <- mclapply(1:length(data), function(x) data[[x]][,-c(3,4)])
}

All_Expr <- do.call(cbind,data)
mRNA <- All_Expr
dim(mRNA)
mode(mRNA)
# remove iterated columns which contain gene names
mRNA <- mRNA[-1,]
gene_id <- mRNA[,1]
length(gene_id)

mRNA <- mRNA[,-seq(1,ncol(mRNA),by=2)] #Keep normalized counts

# To convert mRNA to numeric
mRNA <- lapply(1:ncol(mRNA), function(x) as.numeric(mRNA[,x]))
mRNA <- do.call(cbind,mRNA)

rownames(mRNA) <- gene_id #gene_id is "genesymbol|entrez ID"
colnames(mRNA) <- substr(file_barcode[,2],1,16)


dim(mRNA)
mode(mRNA)

save(mRNA,file=paste(save_path,Type, "_mRNA.RData",sep=""))
