# This script processes and organizes super-enhancer (SE) data.

# Workflow:
# 1. Import BED files for H3K27ac scores from four samples.
# 2. Assign strand information and add sample type and ID columns to each dataset.
# 3. Combine all imported datasets into a single object.
# 4. Filter the combined dataset to separate enhancers and promoters.
# 5. Save the filtered datasets and the combined dataset as RDS files.
#
# Output:
# - 'rds/SEnhancers.rds': RDS file containing super-enhancers.
# - 'rds/Enhancers.rds': RDS file containing enhancers.
# - 'rds/Promoters.rds': RDS file containing promoters.
# - 'rds/AllLiilyTissueIntervals.rds': RDS file containing all combined intervals.

# Load libraries
library(rtracklayer)
library(plyranges)

#setdir
setwd("data/")

# read in the SE, E and P lists
Enh3 <- import('h3k27ac_6samples/3_h3k27ac.scores.bed')
strand(Enh3) <- "*"
#add the "sample_type" column to know after merge where the enhancer from
Enh3$sample_type <- rep('UPPP', length(Enh3))
Enh3$ID <- paste0("3_", as.character(c(1:length(Enh3))))

Enh4 <- import('h3k27ac_6samples/4_h3k27ac.scores.bed')
strand(Enh4) <- "*"
Enh4$sample_type <- rep('UPPP', length(Enh4))
Enh4$ID <- paste0("4_", as.character(c(1:length(Enh4))))

Enh5 <- import('h3k27ac_6samples/5_h3k27ac.scores.bed')
strand(Enh5) <- "*"
Enh5$sample_type <- rep('tumor', length(Enh5))
Enh5$ID <- paste0("5_", as.character(c(1:length(Enh5))))

Enh6 <- import('h3k27ac_6samples/6_h3k27ac.scores.bed')
strand(Enh6) <- "*"
Enh6$sample_type <- rep('tumor', length(Enh6))
Enh6$ID <- paste0("6_", as.character(c(1:length(Enh6))))

AllLillyTissueIntervals <- c(Enh3, Enh4, Enh5, Enh6)

Enhancers <- AllLillyTissueIntervals %>% filter(name != 'promoter') # enhancers and SE
Promoters <- AllLillyTissueIntervals %>% filter(name == 'promoter') # only promoters

SEnhancers <- Enhancers %>% filter(name == 'SE') # only SE

saveRDS(SEnhancers, 'rds/SEnhancers.rds')
saveRDS(Enhancers, 'rds/Enhancers.rds')
saveRDS(Promoters, 'rds/Promoters.rds')

saveRDS(AllLillyTissueIntervals, 'rds/AllLiilyTissueIntervals.rds') # all - promoters,enhancer, SE