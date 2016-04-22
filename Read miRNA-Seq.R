rm(list=ls())
cat("\014")
setwd("ENTER PATH")
library(parallel)

save_path <- "ENTER PATH"
Type <- "Prostate"
#++++++++++++++++++++++++++++++++++++++++++++++++
# Read File_Sample_Map.text to find     +
# matching between sample name and data files	+
#++++++++++++++++++++++++++++++++++++++++++++++++

if (Type=="TGCT") FileSampleMap_path <- "ENTER PATH"
if (Type=="Prostate") FileSampleMap_path <- "ENTER PATH"

file_name <- "FILE_SAMPLE_MAP.txt"

info <- read.csv(paste(FileSampleMap_path,file_name,sep=""),sep="	")
dim(info)
head(info)
# info[,1] contains name of files
# info[,2] contains patient barcode
#+++++++++++++++++++++++++++++++++++++++++++++++++
# Find rows which correspond to miRNA-Seq files
#+++++++++++++++++++++++++++++++++++++++++++++++++
file_id <- 'mirna.quantification'

ind <- which(info[,1] %in% grep(file_id, info[,1], value=TRUE))
file_barcode <- info[ind,]
file_barcode[8,]

#++++++++++++++++++++++++++++++++++
#  Read miRNASeq Files
#++++++++++++++++++++++++++++++++++

if (Type=="TGCT") RNA_path <- "ENTER PATH/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/"
if (Type=="Prostate") RNA_path <- "ENTER PATH/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/"

data <- mclapply(paste(RNA_path,file_barcode[,1],sep=""), function(x) read.table(x,sep="	"))


# removing columns contain reads_per_million_miRNA_mapped  and cross-mapped in case of raw data
data <- lapply(1:length(data), function(x) data[[x]][,-c(3,4)])

miRNA <- do.call(cbind,data)
dim(miRNA)

miRNA[1:4,1:10]
# remove iterated columns which contain gene names
miRNA <- miRNA[-1,]
mi_id <- miRNA[,1]
length(mi_id)

miRNA <- miRNA[,-seq(1,ncol(miRNA),by=2)] #Keep normalized counts

# To convert mRNA to numeric
miRNA <- lapply(1:ncol(miRNA), function(x) as.numeric(miRNA[,x]))
miRNA <- do.call(cbind,miRNA)

rownames(miRNA) <- mi_id #mi_id is miRNA name
colnames(miRNA) <- substr(file_barcode[,2],1,16)

dim(miRNA)

save(miRNA,file=paste(save_path, Type, "_miRNA.RData",sep=""))

