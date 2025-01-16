
# This script performs differential expression analysis using DESeq2 on RNA-seq data.
# The data is related to the JQ1 treatment on the O47 cell line.

# Data Loading:
# - Loads transcript-level quantification data from Salmon (txi object).
# - Reads sample annotations from a CSV file.

# Data Preparation:
# - Cleans and standardizes column names in the txi object.
# - Extracts sample IDs from the file names in the annotations.
# - Filters samples to include only those with a specific condition (O47 cell line).

# Data Filtering:
# - Filters the txi object to include only the relevant samples based on the annotations.

# Differential Expression Analysis:
# - Creates a DESeq2 dataset (ddsTxi) using the txi object and sample annotations.
# - Performs differential expression analysis comparing JQ1 treatment to DMSO control.
# - Extracts the results of the differential expression analysis.

# Results Saving:
# - Saves the differential expression results to an RDS file.

# Load libraries
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
samples<-read.csv('data/Annotations_RNAseqJQ1_fzamune1jhmiedu_144134.csv',sep=",",header=TRUE,stringsAsFactors = FALSE)

samples$ID<-substr(samples$File_Name,13,29)

samples <- samples %>% filter(substr(File_Name,11,11) == "1")

#Filtering O47 cell line
samples <- samples %>% filter(Cell.line.name=="O47")

txi$counts=txi$counts[,samples$ID] #we filter them and order them
txi$abundance=txi$abundance[,samples$ID] #we filter them and order them
txi$length=txi$length[,samples$ID] #we filter them and order them

#deseq itself -- pay attention, columns in txi and samples (rows) in colData are in the same order

samples$condition<-factor(samples$Treatment,levels=c("DMSO (non-treatment control)","JQ1 (the treatment)"))
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~condition)
DE<-DESeq(ddsTxi)
JQ1_vs_DMSO_DE_results<-results(DE,contrast = c("condition","JQ1 (the treatment)","DMSO (non-treatment control)"))

#save results O47 DE
saveRDS(JQ1_vs_DMSO_DE_results,"data/JQ1_vs_DMSO_DE_results_O47_JQ1.rds")
