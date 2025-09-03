# This script processes RNA-seq data using Salmon and tximport.

# Steps:
# 1. Define the directory containing Salmon quantification files.
# 2. Identify all quant.sf files within the specified directory.
# 3. Extract folder names where quant.sf files are located.
# 4. Read the first quant.sf file to extract row names and create a tx2gene mapping.
# 5. Parse gene Ensembl IDs from the row names, omitting version numbers.
# 6. Construct full paths to all Salmon quant.sf files.
# 7. Assign sample names based on folder names.
# 8. Create a tx2gene data frame for aggregating transcripts to genes.
# 9. Import Salmon results using tximport.
# 10. Save the imported data as an RDS file.
#
# Note:
# - The reference transcriptome used is from Gencode release 19.
# - The script assumes a specific directory structure and naming convention for the quant.sf files.

# Load libraries
library(tximport)
library(readr)
#the reference is ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz

#fill in this directory, end with / -- it is the directory we start to look for quant.sf file in
salmonDir = "/data/salmon_files"

#look for them
sf.files<-dir(path=salmonDir,pattern="*quant.sf",recursive=TRUE)
#we want to know the names of the folders where the files are
#there are the second from the end in the strsplit by "/" 
sf.folders<-sapply(lapply(strsplit(sf.files,split="/",fixed=TRUE),rev),"[[",2)

#we read ano of quant.sf to extract row names and to from tx2gene for tximport
tmpmat=read.table(paste0(salmonDir,sf.files[1]),sep="\t",header=TRUE,row.names=1)

#the row names are full id's fron the reference, they look like
#trancript_eensembl_name|gene_embl_name|somthing|something_more|etc
#gene ensemblbl names are in field 2
#transcript.enst=sapply(strsplit(rownames(tmpmat),split="|",fixed=T),"[[",1)
gene.ensg.with.dot=sapply(strsplit(rownames(tmpmat),split="|",fixed=TRUE),"[[",2)
#and, they are like ENSG00121222.3 -- the finishing .N is version number, we omit it
gene.ensg=sapply(strsplit(gene.ensg.with.dot,split=".",fixed=TRUE),"[[",1)

#we need the full path to salmon file for txiimport 
salmonfiles<-paste0(salmonDir,sf.files)
#and their names are the sample names 
names(salmonfiles)<-sf.folders

#it is important -- how do we aggregate trancrits to genes?
tx2gene=cbind(rownames(tmpmat),gene.ensg)
colnames(tx2gene)=c('TXNAME',   'ENSEMBL')
tx2gene=as.data.frame(tx2gene)

#import salmon results
txi <- tximport(salmonfiles, type = "salmon", tx2gene = tx2gene)
#save!!!
saveRDS(txi,"txi_gencode_19_ensembl_DGAY_77_samples.rds")
