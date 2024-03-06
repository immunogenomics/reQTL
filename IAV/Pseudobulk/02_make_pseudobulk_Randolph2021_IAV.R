# 02_make_pseudobulk_Randolph2021_IAV.R

## HOW TO RUN:
## Start interactive job
# bsub -Is -XF -R 'rusage[mem=8000]' -n 4 /bin/bash 
## Activate your conda enviroment: R version 4.1.1 for multiome single cell data
# conda activate newEnv
## RUN: 
# Rscript ~/scratch/02_make_pseudobulk_Randolph2021_IAV.R > 02_make_pseudobulk_Randolph2021_IAV_log.txt 2>&1

## Load libraries
library(presto)
library(DESeq2)
library(edgeR)
library(MASS)
options(stringsAsFactors=FALSE)

# Script to make pseudobulk expression profiles
# I didn't run this as a bash script and just did it in R for TBRU, but this is similar code from AMP

## Create directory for ANALYSIS
main_dir <- '/path_to_directory/Randolph2021/PBMC_RNA'
sub_dir <- 'Pseudobulk'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')

files.data <- '/path_to_directory/Randolph2021/AllPBMC/' 


## Load data
load(paste0('/path_to_directory/Randolph2021/PBMC_RNA/',
	'Randolph2021_PBMC_Clean.Rda') )
# Contains: ge_counts, meta.data


## Make pseudobulk profiles by donor and stimulation condition (NP, P)
all_collapse <- collapse_counts(ge_counts, meta.data, c("SOC_indiv_ID","SOC_infection_status"))

# Save genes (rows) x samples (columns) matrix (from: $count_mat)
out <- as.data.frame(all_collapse$counts_mat)

## Add sample name to Pseudobulk (120 unique samples)
samples <- paste0(all_collapse$meta_data$SOC_indiv_ID,'_',all_collapse$meta_data$SOC_infection_status )
colnames(out) <- samples # 

## write file
#write.table(out, file=gzfile(paste0(files.path,'Randolph2021_Allclean_02_Pseudobulk.txt.gz')), sep = "\t", quote = F, row.names = FALSE)


## GENERATE AVERAGE BY STIM
## remove '_condition' from colnames
tmp <- out

sampleid <- gsub('_.*','',colnames(tmp)) %>% unique()


newout <- data.frame(id=rownames(out))

for(i in sampleid){
	ss <- i
	x <- tmp[,grepl(ss, colnames(tmp)) ]
	x$avg <- (x[,1]+x[,2])/2
	newout <- cbind(newout,x$avg)
	colnames(newout)[ncol(newout)] <- ss
}

rownames(newout) <- newout$id

## NORMALIZATION
sum_exp <- newout[complete.cases(newout$id),] # make sure no missing gene name
sum_exp <- sum_exp[,-1]

# Remove genes that with non-zero expression in <= half of the samples
sum_exp <- sum_exp[rowSums(sum_exp > 0) > .5*ncol(sum_exp),] #  


# Log2 CPM normalization
norm_exp <- log2(cpm(sum_exp)+1)

# Inverse normal transformation
rn<-apply(norm_exp,1,function(x){
    	qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) )
})
rn<-t(rn)

rn <- as.data.frame(rn)
temp <- cbind(gene=rownames(rn),rn)

saveRDS(temp,file=paste0(files.path,'rn_peernormalize_Randolph2021_All_Mean.rds'),version=2 ) # notice using version here


# Extract top genotype PCs (5 top PCs) 
pcfile <- readRDS(paste0(files.data,'Randolph2021_mis0.01_maf5_hwe_PCA.rds'))  
# ID=pcfile$sample.id
pcs <- data.frame(id=pcfile$sample.id ,PC1=pcfile$eigenvec[,1],PC2=pcfile$eigenvec[,2],PC3=pcfile$eigenvec[,3],PC4=pcfile$eigenvec[,4],PC5=pcfile$eigenvec[,5])
## save PC as txt:
write.table(pcs,file=paste0(files.path,'genoPC1_5_Randolph2021_mis0.01_maf5_hwe.txt'),sep="\t",quote=FALSE, row.names = FALSE)



## Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
devtools::session_info()