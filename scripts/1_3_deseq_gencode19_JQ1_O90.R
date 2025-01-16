
# This script performs differential expression analysis using DESeq2 on RNA-seq data.
# The data is filtered to include only the O90 cell line treated with JQ1.
# The results are saved as an RDS file.

# Steps:
# 1. Load RNA-seq data from a preprocessed RDS file.
# 2. Clean and standardize column names in the RNA-seq data.
# 3. Read sample annotations from a CSV file.
# 4. Extract sample IDs and filter for the O90 cell line.
# 5. Filter and order RNA-seq data to match the sample annotations.
# 6. Create a DESeq2 dataset and perform differential expression analysis.
# 7. Save the differential expression results to an RDS file.

#load libraries
library(DESeq2)
library(stringr)
library(dplyr)

#setdir
setwd("data/")

txi<-readRDS('data/txi_gencode_19_ensembl_JQ1_treatment.rds')
#here, we load the results from salmon -JQ1 treatment
####the names in txi are filename-based, we change them to more regular (and, the same with samples$ID)
matched_colnames<-colnames(txi$counts) %>% 
  str_replace("-trimmed$", "") %>% 
  str_replace("HGFCGBCX2-","") %>%
  str_replace("-","~")
colnames(txi$counts)<-matched_colnames
colnames(txi$abundance)<-matched_colnames
colnames(txi$length)<-matched_colnames
#txi colnames are ready

#read samples
samples<-read.csv('RNAseq-gencode19/Annotations_RNAseqJQ1_fzamune1jhmiedu_144134.csv',sep=",",header=TRUE,stringsAsFactors = FALSE)

samples$ID<-substr(samples$File_Name,13,29)

samples <- samples %>% filter(substr(File_Name,11,11) == "1")

#Filtering O90 cell line
samples <- samples %>% filter(Cell.line.name=="O90")

txi$counts=txi$counts[,samples$ID] #we filter them and order them
txi$abundance=txi$abundance[,samples$ID] #we filter them and order them
txi$length=txi$length[,samples$ID] #we filter them and order them

#deseq itself -- pay attention, columns in txi and samples (rows) in colData are in the same order

samples$condition<-factor(samples$Treatment,levels=c("DMSO (non-treatment control)","JQ1 (the treatment)"))
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~condition)
DE<-DESeq(ddsTxi)
JQ1_vs_DMSO_DE_results<-results(DE,contrast = c("condition","JQ1 (the treatment)","DMSO (non-treatment control)"))

#save results O90 DE
saveRDS(JQ1_vs_DMSO_DE_results,"data/JQ1_vs_DMSO_DE_results_O90_JQ1.rds")