# This script prepares domain regions for Super Enhancers (SE), Enhancers (E), and Promoters (P).
# It checks if the required RDS files already exist. If not, it processes the data to create them.

# Workflow:
# 1. Check if the RDS files for SE, E, and P domain ranges exist.
# 2. If the files do not exist:
#    a. Read the 'AllLiilyTissueIntervals.rds' file.
#    b. Filter and process the data for SE, E, and P separately:
#       - Filter intervals by name (SE, enhancer, promoter).
#       - Union the intervals based on origin samples (3, 4, 5, 6).
#       - Assign unique names to the domain ranges.
#       - Identify the source of the domain ranges (tumor, UPPP, mixed).
#    c. Save the processed domain ranges to RDS files.
# 3. If the files exist, read the RDS files for SE, E, and P domain ranges.

#load libraries
library(rtracklayer)
library(dplyr)
library(plyranges)
library(stringr)

todo<-file.exists('data/SE_domains_ranges.rds') & 
  file.exists('data/E_domains_ranges.rds') & 
  file.exists('data/P_domains_ranges.rds')

todo <- !todo

#if the files are not there, we prepre them
#we suppose AllLiilyTissueIntervals.rds to be ready at that moment
#we read SEs and uninon (stich) them
#so for enhancers and to promoters -- 
#all three separately

if (todo) {
  AllLillyIntervals<-readRDS('data/AllLiilyTissueIntervals.rds')
  
  SE_by_sample<-AllLillyIntervals %>% filter(name == "SE")
  origin_sample<-str_split_fixed(SE_by_sample$ID, "_",2)[,1]
  SE_domain_ranges<-GenomicRanges::union(SE_by_sample[origin_sample=="3"],SE_by_sample[origin_sample=="4"])
  SE_domain_ranges<-GenomicRanges::union(SE_domain_ranges,SE_by_sample[origin_sample=="5"])
  SE_domain_ranges<-GenomicRanges::union(SE_domain_ranges,SE_by_sample[origin_sample=="6"])
  names(SE_domain_ranges)=str_c("SE_domain_",as.character(seq(1,length(SE_domain_ranges))))
  tumor_id<-unique(queryHits(findOverlaps(SE_domain_ranges,SE_by_sample %>% filter(sample_type=='tumor'))))
  UPPP_id<-unique(queryHits(findOverlaps(SE_domain_ranges,SE_by_sample %>% filter(sample_type=='UPPP'))))
  SE_domain_ranges$source<-""
  SE_domain_ranges$source[tumor_id]<-"tumor"
  SE_domain_ranges$source[UPPP_id]<-"UPPP"
  SE_domain_ranges$source[intersect(tumor_id,UPPP_id)]<-"mixed"

  E_by_sample<-AllLillyIntervals %>% filter(name == "enhancer")
  origin_sample<-str_split_fixed(E_by_sample$ID, "_",2)[,1]
  E_domain_ranges<-GenomicRanges::union(E_by_sample[origin_sample=="3"],E_by_sample[origin_sample=="4"])
  E_domain_ranges<-GenomicRanges::union(E_domain_ranges,E_by_sample[origin_sample=="5"])
  E_domain_ranges<-GenomicRanges::union(E_domain_ranges,E_by_sample[origin_sample=="6"])
  names(E_domain_ranges)=str_c("E_domain_",as.character(seq(1,length(E_domain_ranges))))
  tumor_id<-unique(queryHits(findOverlaps(E_domain_ranges,E_by_sample %>% filter(sample_type=='tumor'))))
  UPPP_id<-unique(queryHits(findOverlaps(E_domain_ranges,E_by_sample %>% filter(sample_type=='UPPP'))))
  E_domain_ranges$source<-""
  E_domain_ranges$source[tumor_id]<-"tumor"
  E_domain_ranges$source[UPPP_id]<-"UPPP"
  E_domain_ranges$source[intersect(tumor_id,UPPP_id)]<-"mixed"
  
  P_by_sample<-AllLillyIntervals %>% filter(name == "promoter")
  origin_sample<-str_split_fixed(P_by_sample$ID, "_",2)[,1]
  P_domain_ranges<-GenomicRanges::union(P_by_sample[origin_sample=="3"],P_by_sample[origin_sample=="4"])
  P_domain_ranges<-GenomicRanges::union(P_domain_ranges,P_by_sample[origin_sample=="5"])
  P_domain_ranges<-GenomicRanges::union(P_domain_ranges,P_by_sample[origin_sample=="6"])
  names(P_domain_ranges)=str_c("P_domain_",as.character(seq(1,length(P_domain_ranges))))
  tumor_id<-unique(queryHits(findOverlaps(P_domain_ranges,P_by_sample %>% filter(sample_type=='tumor'))))
  UPPP_id<-unique(queryHits(findOverlaps(P_domain_ranges,P_by_sample %>% filter(sample_type=='UPPP'))))
  P_domain_ranges$source<-""
  P_domain_ranges$source[tumor_id]<-"tumor"
  P_domain_ranges$source[UPPP_id]<-"UPPP"
  P_domain_ranges$source[intersect(tumor_id,UPPP_id)]<-"mixed"
  
  saveRDS(SE_domain_ranges,'data/SE_domains_ranges.rds')
  saveRDS(E_domain_ranges,'data/E_domains_ranges.rds')
  saveRDS(P_domain_ranges,'data/P_domains_ranges.rds')
} else {
  SE_domain_ranges <- readRDS('data/SE_domains_ranges.rds')
  E_domain_ranges <- readRDS('data/E_domains_ranges.rds')
  P_domain_ranges <- readRDS('data/P_domains_ranges.rds')
}

export.bed(SE_domain_ranges,con = "data/SE_domains.bed")
export.bed(E_domain_ranges,con = "data/E_domains.bed")
export.bed(P_domain_ranges,con = "data/P_domains.bed")