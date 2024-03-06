# sceQTL_Randolph_IAV_glmm_reQTL_2df_model.R

## SINGLE CELL PME MODEL FOR reQTL WHERE WE TESTED THE SIGNIFICANCE OF TWO INTERACTION TERMS:
# G:SCORE (SCORE=PERTURBATION SCORE)
# G:STIMULUS (NP/P)

## Start interactive job
#bsub -Is -XF -R 'rusage[mem=4000]' -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -n 4 /bin/bash
# conda activate newEnv
# bsub < sceQTL_Randolph_IAV_glmm_reQTL_2df_model.sh 


## Load libraries
library(dplyr)
library(data.table)
data.table::setDTthreads(7)
library(R.utils)
options(stringsAsFactors=FALSE)
library(lme4)
library(Matrix)
library(parallel)
library(future.apply)
library(foreach)
library(doParallel)

## FOREACH PARALLELE
registerDoParallel(cores=5) # match this in shell
options(future.globals.maxSize = 8000 * 1024^2) # to avoid error Global size exceeds maximum allowed size


## Create directory for ANALYSIS
main_dir <- '/path_to_directory'
sub_dir <- 'reQTL_2df_model'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
  print('Directory already exist')
}

files.path <- paste0(output.dir,'/')
files.data <- '/path_to_directory/Inputs/'

# OPEN DATA
load("/path_to_directory/Randolph2021_PBMC_Clean.Rda") 
rm(mrnaNorm);gc() #  ge_counts pca_mrnaScaled  meta.data

# expression PCA data
pca_mrnaScaled <- pca_mrnaScaled$x
row.names(pca_mrnaScaled) <- row.names(meta.data)

## CELLTYPES: we will consider B, CD4_T, CD8_T, monocytes, and NK
celltypes <- c('B', 'CD4_T', 'CD8_T', 'monocytes', 'NK')


## FILTER DATA
meta.data <- meta.data[meta.data$celltype %in% celltypes, ]
meta.data <- meta.data[meta.data$SOC_indiv_ID != "HMN52545",] # remove sample with no genotype data


## DEFINE DISCRETE PERTURBATION ('Status01') AND SCALE IT (Mean=0, Var=1)
meta.data$Status01 <- ifelse(meta.data$SOC_infection_status=='NI',0,1)
meta.data$Status01 <- scale(meta.data$Status01)


## OPEN PERTURBATION SCORE
lda <- readRDS(paste0('/path_to_directory/LDA/',
  'PerturbationScore_Randolph2021_IAV.rds')) 

### ADD PERTURBATION SCORE TO META.DATA
meta.data$lda <- lda[row.names(meta.data),'s0']
meta.data$SCORE <- scale(meta.data$lda) # SCALED SCORE

## REMOVE THE EFFECT OF Status01 SO WE CAN INCLUDE BOTH PERTURBATION VARIABLES INTO THE reQTL MODEL
#(*) rSCORE represents the independent effect of SCORE given its correlation with Status01
meta.data$rSCORE <- summary(lm(SCORE~Status01,data=meta.data))$residuals

rm(lda); gc()


## ALLELE DOSAGE GENOTYPE DATA 
geno <- fread(paste0(files.data,'Randolph2021_PBMC_Mean_Flu_Mock_sig_mis0.01_maf5_hwe_DS.gz'),stringsAsFactors=FALSE) %>% as.data.frame()


## Individual Meta data (Age)
meta <- read.csv(paste0('/data/srlab2/anathan/external_data/Randolph2021/','Randolph_Table_S1_Sample meta data.csv') )
row.names(meta) <- meta$infection_ID  #meta$indiv_ID
meta.data$age <- meta[meta.data$sample_condition,"age"] 


## SINGLE CELL EXPRESSION DATA
exprs_raw = ge_counts
class(exprs_raw)
dim(exprs_raw)

## Keep cells that pass QC
exprs_raw <- exprs_raw[,row.names(meta.data)] #


## Selected eGENES (FROM PSEUDOBULK ANALYSIS WITH FASTQTL)
genes <- readRDS('/data/srlab2/cristian/Randolph2021/AllPBMC/SCeQTLinteraction/inputs/GoodGenes_Results_fastqtl_permute_Randolph2021_MEAN_ALL.rds')


## MAKE SURE TO HAVE THE SAME eGENES IN EXPRESSION DATA AND genes
exprs_raw <- exprs_raw[row.names(exprs_raw) %in% intersect(genes$Gene, row.names(exprs_raw)),]
genes <- genes[genes$Gene %in% row.names(exprs_raw),]


## Genotype PCs
pcs <- read.table("/data/srlab2/anathan/external_data/Randolph2021/genoPC1_5_Randolph2021_rawVCF.txt", header = T)
row.names(pcs) <- pcs$id
pcs <- pcs[,2:ncol(pcs)]


## Expression PCs 
pca_mrnaScaled <- pca_mrnaScaled[row.names(meta.data) ,]


## Define number of test
end <-  nrow(genes)

## RUN TESTS
system.time(
res <- foreach(i=1:end, .combine=rbind) %dopar% {  
  E <- as.numeric(exprs_raw[genes$Gene[i],]) 
  G <- t(subset(geno, snp == genes$snp[i] ) [,intersect(colnames(geno), as.character(meta.data$sample_condition))])[,1]

  if(length(table(round(as.numeric(as.character(G))))) >=2 & sum(E) > 0){ #at least 3 different genotype is required. (sample seize is so small, this filter helps)
    IND <- factor(meta.data$sample_condition) #individuals as factor covariates  
    nUMI <- scale(log(meta.data$nCount_RNA)) # offset
    Status <- meta.data$Status01 # DISCRETE PERTURBATION
    PC <- pcs[gsub("_.*", "", as.character(IND)),1:5] #
    MT <- meta.data$percent.mt/100
    B <- meta.data$batchID
    expPC <- pca_mrnaScaled[,1:5]
    AGE <- scale(meta.data$age)
    SCORE <- meta.data$rSCORE # CONTINOUS PERTURBATION SCORE (*residual effect)

    data <- data.frame(E,IND,B,Status,AGE,nUMI,MT,PC1 = PC[,1],PC2 = PC[,2],PC3 = PC[,3],PC4 = PC[,4],PC5 = PC[,5],expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], SCORE)

    data$G <- t(subset(geno, snp == genes$snp[i] ) [,as.character(data$IND)] ) [,1]
    data$G <- as.numeric(as.character(data$G))

    data$IND <- gsub('_.*','',data$IND) %>% as.factor() # remove '_condition' from donor 

    
    x <- data.frame(gene=genes$Gene[i],snp=genes$snp[i],term=NA, Estimate=NA, Std.Error=NA, zvalue=NA, pval=NA, lrt_pval=NA)
    tryCatch({
    test <- lme4::glmer(formula = paste0("E ~ G + (1 | B) + (1 | IND) + Status + AGE + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 +  SCORE + G:SCORE + G:Status", collapse = " + "), family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap") ) 
    test_null <- lme4::glmer(formula = paste0("E ~ G + (1 | B) + (1 | IND) + Status + AGE + nUMI + MT + PC1 + PC2 + PC3 + PC4 + PC5 + expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + SCORE", collapse = " + "), family = "poisson", nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap") ) 
    model_lrt <- anova(test_null, test)

    x <- summary(test)$coefficients
    colnames(x) <- c("Estimate","Std.Error","zvalue","pval")
    x <- data.frame(gene=genes$Gene[i],snp=genes$snp[i],term=row.names(x),x, lrt_pval=model_lrt$`Pr(>Chisq)`[2])
    #x <- x[grep("^G",x$term),] # Select only main G effect and Interaction?
    }, error=function(e){} )
  }
} 
)
## Save results
write.table(res, file=gzfile(paste0(files.path,'sceQTL_Randolph_IAV_glmm_reQTL_2df_model.txt.gz')),row.names=FALSE,quote=FALSE )

rm(list=ls())

## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()