# Results_fastqtl_permute_Randolph2021_IAV_MEAN.R

## GENERATE A FILE WITH ALL SIGNIFICANT eQTL FROM JOINT PBMCS USING AVERAGE PSEUDOBULK PROFILE (FLU+MOCK)/2
## SIGNIFICANT WAS DEFINED BY QVALUE <=0.10 OF THE PERMUTATION P VALUE FROM FASTQTL 

## HOW TO RUN:
## Start interactive job

## Activate your conda enviroment: R version 4.1.1 for multiome single cell data
# conda activate newEnv
## RUN: 
# Rscript /path_to_directory/Results_fastqtl_permute_Randolph2021_IAV_MEAN.R > Results_fastqtl_permute_Randolph2021_IAV_MEAN_log.txt 2>&1

# NO: # bsub -Is -XF -R 'rusage[mem=5000]' -n 4 /bin/bash

## Load packages
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(R.utils) # to open gz files 
library(qvalue)
library(data.table)
data.table::setDTthreads(7)
options(stringsAsFactors=FALSE)


## Create directory for ANALYSIS 
main_dir <- '/path_to_directory/fastqtl_permute'
sub_dir <- 'eQTL'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')

## FILES:
all_df <- list.files(path="/path_to_directory/fastqtl_permute",pattern=".permute.txt.gz",full.names=TRUE)


## Empty df
df <- data.frame()

## add sign eQTLS
for(i in all_df){
	tmp <- fread(i)
	tmp <- as.data.frame(tmp)
	colnames(tmp) <- c("Gene", "nVar", "MLE1", "MLE2", "Dummy", "snp", "Dist", "NomPvalue", "Beta", "PermPvalue", "Pvalue") 
	tmp$permQ <- qvalue::qvalue(tmp$Pvalue)$qvalue
	tmp <- tmp[tmp$permQ<0.10,] 
	df <- rbind(df,tmp)
}


## SNP Info: chr, position, snp ID, reference, alternative
vcf <- fread('/path_to_directory/SNPinfo_VCF.txt.gz') %>% as.data.frame()
colnames(vcf) <- c('chr','pos','snp','ref','alt')

## Merge
df <- merge(df,vcf,by='snp')
df <- distinct(df, Gene,snp,Pvalue, .keep_all = TRUE)
df <- df[, c('Gene','chr','pos','snp','ref','alt', "nVar", "MLE1", "MLE2", "Dummy", "Dist", "NomPvalue", "Beta", "PermPvalue", "Pvalue" )]


## Write rds file
saveRDS(df,file=paste0(files.path,'Results_fastqtl_permute_Randolph2021_IAV_MEAN.rds') )

  
## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()