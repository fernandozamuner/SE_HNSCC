# This script processes RNA-seq data using the tximport package to import transcript-level quantification data from Salmon
# and aggregates it to the gene level. The resulting data is saved as an RDS file for further analysis.

# Load libraries
library(tximport)
library(readr)

# Specify the directory containing Salmon quantification files
salmonDir <- "data/salmon_output/"

sf.files<-dir(path=salmonDir,pattern="*quant.sf",recursive=TRUE)
sf.folders<-sapply(lapply(strsplit(sf.files,split="/",fixed=TRUE),rev),"[[",2)

tmpmat=read.table(paste0(salmonDir,sf.files[1]),sep="\t",header=TRUE,row.names=1)

gene.ensg.with.dot=sapply(strsplit(rownames(tmpmat),split="|",fixed=TRUE),"[[",2)
gene.ensg=sapply(strsplit(gene.ensg.with.dot,split=".",fixed=TRUE),"[[",1)

salmonfiles<-paste0(salmonDir,sf.files)
names(salmonfiles)<-sf.folders

tx2gene=cbind(rownames(tmpmat),gene.ensg)
colnames(tx2gene)=c('TXNAME',   'ENSEMBL')
tx2gene=as.data.frame(tx2gene)
txi <- tximport(salmonfiles, type = "salmon", tx2gene = tx2gene)
saveRDS(txi,"txi_gencode_19_ensembl_JQ1_treatment.rds")
