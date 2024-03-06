# Lasso_PerturbationScore_Randolph2021_IAV.R

## THIS SCRIPT GENERATES THE CONTINUOUS PERTURBATION SCORE USING A PENALIZED LOGISTIC REGRESSION 'LASSO REGRESSION'
## FOR INFECTION STATUS (MOCK/FLU) OF SINGLE CELLS FROM RANDOLPH et al
## LASSO REGRESSION IS BASED ON 20 HARMONYZED PCs

## Start interactive job
#bsub -Is -XF -R 'rusage[mem=4000]' -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -n 4 /bin/bash
# conda activate newEnv
## RUN: 
# Rscript /path_to_directory/Lasso_PerturbationScore_Randolph2021_IAV.R > Lasso_PerturbationScore_Randolph2021_IAV_log.txt 2>&1


# Load libraries
library(dplyr)
library(data.table)
data.table::setDTthreads(7)
library(R.utils)
library(MASS) # for LDA
library(glmnet) # for lasso
library(ggplot2)
library(ggrastr)
library(grid)
library(gridExtra)
theme_set(theme_classic()) # 
options(stringsAsFactors=FALSE)
source('~/10x_Multiome_GEaTac/10x_PBMC/Rfunctions.R') #

## Create directory for ANALYSIS
main_dir <- '/path_to_directory/'
sub_dir <- 'LDA'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')
files.data <- '/path_to_directory/Inputs/'

## FOR REPRODUCIBILITY
set.seed(100)

f <- function(x) factor(x, levels = unique(x))

# Load QC scMeta.data all PBMCs ; rows=cell that pass qc 208223 cells (we need age)
load("/path_to_directory/Randolph2021_PBMC_Clean.Rda") #  ; Randolph2021_MOCK_Clean.Rda
rm(ge_counts, mrnaNorm);gc() #  Keep "pca_mrnaScaled"  "meta.data"

## CELLTYPES: we will consider B, CD4_T, CD8_T, monocytes, and NK
celltypes <- c('B', 'CD4_T', 'CD8_T', 'monocytes', 'NK')

## FILTER DATA 202035 cells
## Select cell from ge_counts, meta.data
meta.data <- meta.data[meta.data$celltype %in% celltypes, ]
meta.data <- meta.data[meta.data$SOC_indiv_ID != "HMN52545",]

## Select columns or variables: Top 20 Harmony PCs and Perturbation condition (NP, P)
vars <- c(colnames(meta.data)[grep('harmonize',colnames(meta.data) )],'SOC_infection_status' )

## SELECT VARS AND ALL CELLS 
lda.set <- meta.data[,vars]


## LDA: Lasso regression

# find best lambda for cross-validation
cv.lasso <- cv.glmnet(as.matrix(lda.set[,-21]), lda.set$SOC_infection_status, alpha=1,family='binomial', type.measure='auc')

## PLOT Figure S
pdf(paste0(files.path,'Crossvalidation_LDALasso_Randolph2021_PBMC_Clean.pdf'),w=5,h=3.5)
plot(cv.lasso)
dev.off()


## WE USE THIS ONE BELOW: 0.002311521
# lambda 'lambda.1se' gives the simplest model but also lies within one standard error of the optimal value of lambda
cv.lasso$lambda.1se 

## COEFICIENTS
lasso_coef <- coef(cv.lasso, cv.lasso$lambda.1se)
lasso_coef <- as.data.frame(lasso_coef[-1,])
lasso_coef$PC <- gsub('.*_','',rownames(lasso_coef))
lasso_coef$PC <- f(lasso_coef$PC)

## PLOT: Figure 1
pdf(paste0(files.path,'Coef_LassoRegression_Randolph2021_IAV.pdf'), width=6,height=3)
ggplot(lasso_coef,aes(x=PC,y=lasso_coef[,1])) + geom_bar(stat='identity') +
labs(title='Randolph2021', x='Harmony PC' ,y='Coefficient from Lasso regression' ) +
theme(legend.position = "bottom",plot.title=element_text(size = 12,face="bold") ,
		axis.title=element_text(size=12,face="bold"), axis.text = element_text(size = 12, face="bold"),
		axis.text.x=element_text(angle=90) )
dev.off()

## Lasso model 'simplest model'
lasso.model <- glmnet(as.matrix(lda.set[,-21]), lda.set$SOC_infection_status, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.1se)

# Make prediction based on lasso 
lasso_pred <- lasso.model %>% predict(newx=as.matrix(lda.set[,-21]) )
lasso_pred <- as.data.frame(lasso_pred)

## SAVE THIS FILE
saveRDS(pred_lda, file=paste0(files.path,'PerturbationScore_Randolph2021_IAV.rds'))

## PLOT: Figure 1 c
lasso_pred$SOC_infection_status <- lda.set[row.names(lda.set),'SOC_infection_status']
lasso_pred$SOC_infection_status <- factor(lasso_pred$SOC_infection_status,levels=c('NI','flu'))

pdf(paste0(files.path,'lassoPerturbationScore_Hist_Randolph2021_IAV.pdf'),w=3,h=3)
ggplot(lasso_pred, aes(x=SCORE, color=SOC_infection_status)) + geom_histogram(fill="white", alpha=0.5, position="identity") + 
scale_color_manual(name='IAV perturbation',labels=c('NP','P'), values=c('#FFCC00','#AA4499')) +
labs(title ='', x='Perturbation score') +
theme(legend.position = "bottom",
    plot.title = element_text(color="black", size=14, face="italic"),
    axis.title=element_text(size=14,face="bold"), 
    axis.text = element_text(size = 12, face="bold"),
    legend.key.size = unit(0.5, "cm")) 
dev.off()



## Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
devtools::session_info()