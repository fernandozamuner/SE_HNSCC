#' Analysis of Differential Gene Expression using DESeq2
#' 
#' This script performs differential expression analysis on RNA-seq data comparing
#' tumor (T) versus normal (N) samples using DESeq2 package.
#' 
#' The workflow includes:
#' 1. Reading pre-processed Salmon quantification results stored in txi format
#' 2. Filtering samples to include only tumor (T) and normal (N) samples
#' 3. Running DESeq2 analysis
#' 4. Saving differential expression results
#' 
#' @param txi Input transcriptome expression data from Salmon
#' @param samples Sample metadata containing condition information
#' 
#' @return DESeq2 results object comparing tumor vs normal samples
#' 
#' @note Results are saved as 'T_to_N_DE_results_hg19_ensembl_77_samples.rds'
#' @note Uses GENCODE 19/Ensembl annotation
#' @note Analysis includes 77 samples


#load libraries
library(DESeq2)
library(stringr)
library(dplyr)

#setdir
# changing path
setwd("data/")

txi<-readRDS('data/txi_gencode_19_ensembl_DGAY_77_samples.rds')
#here, we load the results from salmon
#the names in txi are filename-based, we change them to more regular (and, the same with samples$ID)
matched_colnames<-colnames(txi$counts) %>% str_replace("_trimed$", "") %>% str_replace("_2","")
colnames(txi$counts)<-matched_colnames
colnames(txi$abundance)<-matched_colnames
colnames(txi$length)<-matched_colnames
#txi colnames are ready

#read samples
samples<-read.csv('data/DGAY_77_sampleID_mapping.txt',sep='\t',header=TRUE,stringsAsFactors = FALSE)
#we nee only T and N for this run
samples <- samples %>% filter(class=="T" | class=="N")

txi$counts=txi$counts[,samples$ID] #we filter them and order them
txi$abundance=txi$abundance[,samples$ID] #we filter them and order them
txi$length=txi$length[,samples$ID] #we filter them and order them


#deseq itself -- pay attention, columns in txi and samples (rows) in colData are in the same order
samples$condition<-factor(samples$class,levels=c("N","T"))
ddsTxi <- DESeqDataSetFromTximport(txi,colData = samples,design = ~condition)
DE<-DESeq(ddsTxi)
T_to_N_DE_results<-results(DE,contrast = c("condition","T","N"))

#save results
saveRDS(T_to_N_DE_results,"RNAseq-gencode19/T_to_N_DE_results_hg19_ensembl_77_samples.rds")

