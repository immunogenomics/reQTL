# 03_peernormalize_Mean_Randolph2021_IAV.R


## HOW TO RUN:
## Start interactive job
## Activate your conda enviroment: Notice this env is specific for Peer
# conda activate peer

## RUN: 
# bsub < 03_peernormalize_Mean_Randolph2021_IAV.sh

# Script to PEER normalize pseudobulk samples
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865505/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398141/

## Load libraries
library(peer)


## Create directory for ANALYSIS
main_dir <- '/path_to_directory/Randolph2021/PBMC_RNA/Pseudobulk'
sub_dir <- 'AllPBMC'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')
files.data <- '/path_to_directory/Randolph2021/'

# Replace with reading in your donor-level covariate matrix: Randolph_Table_S1_Sample meta data.csv
pheno <- read.csv(paste0(files.data,'Randolph_Table_S1_Sample meta data.csv') )
# select any condition (no problem)
pheno <- pheno[pheno$infection_status=='NI',]

rownames(pheno) <- pheno$indiv_ID #pheno$infection_ID

# Replace with reading in your genotype PCs (5 top PCs)
pcs <- read.table(paste0(main_dir,'/genoPC1_5_Randolph2021_rawVCF.txt'),header=T)
row.names(pcs) <- pcs$id
pcs <- pcs[,-1]


# PEER Factor Analysis: regress out latent variables (Stegle 2010, PLoS Comp Bio)

# Based on GTEx 2020 Science/Battle 2017 Nature:
K = 15


# Replace with reading in your pseudobulk expression matrix
rn <- readRDS(paste0(files.path,'rn_peernormalize_Randolph2021_All_Mean.rds') )
rn <- rn[complete.cases(rn$gene),]
row.names(rn) <- rn$gene 
rn <- rn[,-1]

## MATCH IDS: pcs(91),pheno(90),rn(90) infection_ID indiv_ID 
samples <- Reduce(intersect,list(colnames(rn),row.names(pcs),unique(pheno$indiv_ID)) ) # 89 samples

## samples in pheno
pheno <- pheno[pheno$indiv_ID %in% samples, ]
pcs <- pcs[row.names(pcs) %in% samples, ]
rn <- rn[, colnames(rn) %in% samples]


# Select known covariates to regress out (here, age, sex, and 5 genotype PCs)
# note: all donors are male
covs <- cbind(pheno[colnames(rn),'age'],
			#pheno[colnames(rn),'Sex'],
            pcs[colnames(rn),1:5])

colnames(covs) <- c("age", "PC1", "PC2", "PC3", "PC4", "PC5") # 

model = PEER()
PEER_setPhenoMean(model, as.matrix(t(rn)))
dim(PEER_getPhenoMean(model))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,K)
PEER_getNk(model)
PEER_setCovariates(model, as.matrix(covs))
PEER_setNmax_iterations(model,10000) # 
    
#perform the inference
PEER_update(model)
# Converged (var(residuals)) 
    
# Save residuals
residuals = PEER_getResiduals(model)
colnames(residuals)  <- row.names(rn)
row.names(residuals) <- colnames(rn)
dump <- data.frame( ids = row.names(rn),t(residuals))

gz1 <- gzfile(paste0(files.path, "PeerNorm_Mean_Pseudobulk_residuals_Randolph2021_All.txt.gz"),"w")
write.table(dump, gz1, sep = "\t", quote = F, row.names = FALSE)
close(gz1)

## Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
devtools::session_info()