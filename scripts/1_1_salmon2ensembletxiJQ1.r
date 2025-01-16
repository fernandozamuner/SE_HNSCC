# This script processes RNA-seq data using the tximport package to import transcript-level quantification data from Salmon
# and aggregates it to the gene level. The resulting data is saved as an RDS file for further analysis.
#
# Steps:
# 1. Specify the directory containing the Salmon quantification files.
# 2. Identify all Salmon quantification files in the specified directory.
# 3. Extract folder names corresponding to the Salmon quantification files.
# 4. Read the first Salmon quantification file to extract gene and transcript information.
# 5. Extract gene Ensembl IDs without version numbers.
# 6. Create a data frame mapping transcript names to gene Ensembl IDs.
# 7. Import the Salmon quantification files and aggregate the data to the gene level using tximport.
# 8. Save the aggregated data as an RDS file for further analysis.

# Load libraries
library(tximport)
library(readr)

#fill in this directory, end with a
salmonDir = " "

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
