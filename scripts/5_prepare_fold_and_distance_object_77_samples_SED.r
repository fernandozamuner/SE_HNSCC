#' This script prepares a fold and distance object for 77 samples by processing 
#' super-enhancer (SE) domains and gene annotations, and integrating RNA-seq differential expression results.
#' 
#' @details 
#' The script performs the following steps:
#' 
#' 1. **Load Libraries**: Loads necessary libraries including `rtracklayer`, `dplyr`, `plyranges`, and `differential.coverage`.
#' 2. **Set Working Directory**: Sets the working directory to "data/".
#' 3. **Load SE Domains**: Loads SE domains from an RDS file and calculates the midpoint of each SE domain.
#' 4. **Get Gene Annotations**: Retrieves gene annotations from the `differential.coverage` package and filters out genes with non-unique names on chromosome Y.
#' 5. **Calculate TSS**: Calculates the Transcription Start Site (TSS) for each gene.
#' 6. **Calculate Distances**: Uses `distanceToNearest` to calculate the distance from each TSS to the nearest SE domain midpoint.
#' 7. **Prepare Data Frame**: Prepares a data frame containing distances, SE domain types, specificity, and IDs.
#' 8. **Integrate RNA-seq Data**: Loads RNA-seq differential expression results, removes rows with NA values, and integrates these results with the SE domain data.
#' 9. **Save Results**: Saves the final fold and distance object to an RDS file.

###Load libraries
library(rtracklayer)
library(dplyr)
library(plyranges)

#setdir
setwd("data/")

# load libraries
if (! require(differential.coverage)) {
  devtools::install_github("favorov/differential.coverage") 
  library(differential.coverage)
}

#load SEs
SED<-readRDS('data/SE_domains_ranges.rds')
#we need middle of SE to see distances
SEDmids<-SED
start(SEDmids)<-(start(SED)+end(SED))/2
width(SEDmids)<-0

#get the_genes
the_genes<-differential.coverage::get.Known.Gene.List(genome.annotation.id = 'gencode19')
#some genes exist in chrX and chrY
#we will leave only chrX coordinates
#count how many times a gene name appears
names_abundance_table<-table(the_genes$ensembl)
nonunique_names<-names(which(names_abundance_table>1))
#the unique gene names are OK, nonumique non-Y, too
the_genes <- the_genes %>% filter((!ensembl %in% nonunique_names) | seqnames!="chrY")

#and, TSS for each gene
TSS<-the_genes
start(TSS)<-ifelse(as.character(strand(the_genes))=='-',end(the_genes),start(the_genes))
width(TSS)<-0
###distanceToNearest: Returns the distance for each range in x to its nearest neighbor in the subject.

the_genes_to_closest_SED<-distanceToNearest(TSS,SEDmids,ignore.strand=TRUE)

distances <- mcols(the_genes_to_closest_SED)$distance

TSS_in_pairs<-TSS[queryHits(the_genes_to_closest_SED)]
SED_in_pairs<-SEDmids[subjectHits(the_genes_to_closest_SED)]

is_strand_down <- as.character(strand(TSS_in_pairs))=='-'
is_enhacer_upstream <- start(TSS_in_pairs) < start(SED_in_pairs)
is_distance_positive <- is_enhacer_upstream==is_strand_down

# Add a new column "specific" for "tumor" and "UPPP"
SED_in_pairs$specific <- ifelse(SED_in_pairs$source %in% c("tumor", "UPPP"), "specific", 
                                ifelse(SED_in_pairs$source == "mixed", "non_specific", NA))

SED_in_pairs$ID <- as.character(SED_in_pairs@ranges@NAMES)

closest_SED_per_gene_geometry <- data.frame(
  distance=ifelse(is_distance_positive,distances,-distances),
  SED_type=SED_in_pairs$source,SED_is_specific=SED_in_pairs$specific,
  SED_ID=SED_in_pairs$ID,stringsAsFactors=FALSE)

##Adding "ensembl" column to "closest_SED_per_gene_geometry"
rownames(closest_SED_per_gene_geometry)=TSS_in_pairs$ensembl

##Adding "gene_name" column to "closest_SED_per_gene_geometry"
closest_SED_per_gene_geometry$gene_name=TSS_in_pairs$gene_name

#===========================
# RNAseq data for 77 samples
#===========================	
T_to_N_DE_results<-readRDS('data/T_to_N_DE_results_hg19_ensembl_77_samples.rds')

# Remove rows with NA values
T_to_N_DE_results <- T_to_N_DE_results[complete.cases(T_to_N_DE_results), ]

#==============================
# we test names correspondence:
DE_names<-rownames(T_to_N_DE_results)
is_DE_name_bad <- is.na(closest_SED_per_gene_geometry[DE_names,1])
if (sum(is_DE_name_bad)>0) {warning(paste0(as.character(sum(is_DE_name_bad))," unannotated the_genes"))}

fold_and_distance <- as.data.frame(T_to_N_DE_results[!is_DE_name_bad,c("baseMean","log2FoldChange","lfcSE","padj")],stringsAsFactors=FALSE)
fold_and_distance <- cbind(fold_and_distance,closest_SED_per_gene_geometry[rownames(fold_and_distance),])

saveRDS(fold_and_distance,'data/fold_and_distance_77_samples_SED.rds')