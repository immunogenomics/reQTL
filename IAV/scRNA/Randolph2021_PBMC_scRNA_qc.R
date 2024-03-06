# Randolph2021_PBMC_scRNA_qc.R

## HOW TO RUN:
## Start interactive job
# bsub -Is -XF -R 'rusage[mem=5000]' -n 4 /bin/bash 

## Activate your conda enviroment: R version 4.1.1 
# conda activate newEnv
## RUN: 
# Rscript ~/scratch/Randolph2021_PBMC_rna_qc.R > Randolph2021_PBMC_rna_qc_log.txt 2>&1


## Load libraries
library(Seurat)
library(Signac)
library(singlecellmethods)
library(harmony)
library(parallel)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggrastr)
theme_set(theme_classic())
library(grid)
library(gridExtra)
source('~/10x_Multiome_GEaTac/10x_PBMC/Rfunctions.R')
options(stringsAsFactors=FALSE)

## Create directory for ANALYSIS /path_to_directory/Randolph2021/
main_dir <- '/path_to_directory/Randolph2021'
sub_dir <- 'PBMC_RNA'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')
files.data <- '/path_to_directory/Randolph2021/'

## Open files: expr_matrix  and metadata
# Randolph_Table_S1_Sample meta data.csv ; mergedAllCells_withCellTypeIdents_CLEAN.rds
# Raw counts
ge_counts <- readRDS(paste0(files.data,'randolph_raw.rds')) 
# Meta data
meta.data <- readRDS(paste0(files.data,'randolph_meta.rds')) 

## REMOVE IAV GENES from Gene count matrix
iav_genes <- c('PB2','PB1','PA','HA','NP','NA','M2','M1','NEP','NS1')

ge_counts <- ge_counts[!rownames(ge_counts) %in% iav_genes,]

## Estimate % of reads mapped to mito genes and add it to meta.data?
# Already in meta.data and NO NEED TO FILTER
mito_genes <- grep("^MT-", rownames(ge_counts), value = TRUE, ignore.case = TRUE)
#meta.data$percent_mitolab <- Matrix::colSums(ge_counts[mito_genes, ])/Matrix::colSums(ge_counts)
#plot(meta.data$percent_mitolab,meta.data$percent.mt)
#nrow(meta.data[meta.data$percent_mitolab <=0.1, ]) #
#meta.data <- meta.data[meta.data$percent_mitolab <=0.1, ]

## Ribosomal reads
genes_ribo <- grep("^RPS\\d+|^RPL\\d+", row.names(ge_counts), value = TRUE)
length(genes_ribo) # 99 meta.data$cell
meta.data$percent.ribo <- Matrix::colSums(ge_counts[genes_ribo, rownames(meta.data)]) / 
                                        Matrix::colSums(ge_counts[, rownames(meta.data)])
meta.data$percent.ribo <- meta.data$percent.ribo * 100


## NUMBER OF CELLS WITH < 500 GENES BY CONDITIONS
FU_cells <- nrow(meta.data[meta.data$SOC_infection_status=='flu' & meta.data$nFeature_RNA>500,])
NI_cells <- nrow(meta.data[meta.data$SOC_infection_status=='NI' & meta.data$nFeature_RNA>500,])

print('DISTRIBUTION OF UMIS PER CELL BY CONDITION')
tapply(meta.data$nCount_RNA,meta.data$SOC_infection_status,summary)


## Plot each feature by 'infection_status':
pdf(paste0(files.path,'Randolph2021_PBMC_qc.pdf'),w=10,h=10)
(ggplot(meta.data,aes(x=SOC_infection_status,y=nCount_RNA)) + geom_violin(fill = "lightpink2") + rasterize(geom_jitter(position=position_jitter(0.2),size = 0.5, stroke = 0.1, shape = 21 ),dpi=300) +
xlab('') + ylab('Number of UMIs per cell') +
    theme(legend.position = "none")
    )+
(ggplot(meta.data,aes(x=SOC_infection_status,y=nFeature_RNA)) + geom_violin(fill = "lightpink2") + rasterize(geom_jitter(position=position_jitter(0.2),size = 0.5, stroke = 0.1, shape = 21 ),dpi=300) +
geom_hline(aes(yintercept=500 ),color="red", linetype="dashed", size=0.5 ) +
#annotate("text", x ='Randolph2021', y =  236972, label =F_cells, color='red' ) +
xlab(paste0('Cells with >500 genes=',FU_cells,'/',NI_cells)) + ylab('Number of genes detected per cell') + 
    theme(legend.position = "none")
    )+
(ggplot(meta.data, aes(x=SOC_infection_status, y=percent.mt)) + geom_violin(fill = "lightpink2") + rasterize(geom_jitter(position=position_jitter(0.2),size = 0.5, stroke = 0.1, shape = 21 ),dpi=300) + 
  	xlab('') + ylab('Percentage of mito')+ geom_hline(aes(yintercept=20 ),color="red", linetype="dashed", size=0.5 ) +
  	theme(legend.position = "none" )
  	)+
(ggplot(meta.data, aes(x=SOC_infection_status, y=percent.ribo)) + geom_violin(fill = "lightpink2") + rasterize(geom_jitter(position=position_jitter(0.2),size = 0.5, stroke = 0.1, shape = 21 ),dpi=300) +
  	xlab('') + ylab('Percentage of ribo') +
  	theme(legend.position = "none" )
  	)+
plot_layout(nrow = 2, byrow = TRUE)
dev.off()

## PLOT X vs Y
pdf(paste0(files.path,'Randolph2021_PBMC_FeatVsFeat_qc.pdf'),w=12,h=8)
(ggplot(meta.data,aes(x=percent.ribo,y=percent.mt,color=SOC_infection_status)) + rasterize(geom_point(size = 0.4,stroke=0.25, shape = 21),dpi=300 ) + xlab('Percentage of ribo') + ylab('Percentage of mito') +
#geom_hline(aes(yintercept=20 ),color="red", linetype="dashed", size=0.5 ) +
    guides(color=guide_legend(title="Infectious status",override.aes = list(size = 5) ) ) +
theme(
    #legend.position = "none",
    axis.title.x = element_text(size = 16,face="bold"),
    axis.title.y = element_text(size = 16,face="bold"),
    axis.text.x = element_text(face="bold", size=14),
    axis.text.y = element_text(face="bold", size=14,) )
    ) +
(ggplot(meta.data,aes(x=nFeature_RNA,y=percent.mt,color=SOC_infection_status)) + rasterize(geom_point(size = 0.4,stroke=0.25, shape = 21),dpi=300 )  + xlab('Number of genes detected per cell') + ylab('Percentage of mito') +
geom_vline(aes(xintercept=500 ),color="red", linetype="dashed", size=0.5 ) +
guides(color=guide_legend(title="Infectious status",override.aes = list(size = 5) ) ) +
theme(
    #legend.position = "none",
    axis.title.x = element_text(size = 16,face="bold"),
    axis.title.y = element_text(size = 16,face="bold"),
    axis.text.x = element_text(face="bold", size=14),
    axis.text.y = element_text(face="bold", size=14,) ) 
)
(ggplot(meta.data,aes(x=nFeature_RNA,y=percent.ribo,color=SOC_infection_status)) + rasterize(geom_point(size = 0.4,stroke=0.25, shape = 21),dpi=300 )  + xlab('Number of genes detected per cell') + ylab('Percentage of ribo') +
geom_vline(aes(xintercept=500 ),color="red", linetype="dashed", size=0.5 ) +
guides(color=guide_legend(title="Infectious status",override.aes = list(size = 5) ) ) +
theme(
    #legend.position = "none",
    axis.title.x = element_text(size = 16,face="bold"),
    axis.title.y = element_text(size = 16,face="bold"),
    axis.text.x = element_text(face="bold", size=14),
    axis.text.y = element_text(face="bold", size=14,) )
    ) +
plot_layout(nrow = 1, byrow = TRUE)
dev.off()

## FIRST REMOVE CELL WITH MITO MAP > 20%? 
# No: DATA WAS ALREADY FILTERED BY MITO < 10%


## FILTER CELLS WITH DETECTED GENES < 500
meta.data <- meta.data[meta.data$nFeature_RNA>500,] # 28,395 cell removed and 208,223 final PBMC   

# ge_counts
ge_counts <- ge_counts[ ,colnames(ge_counts) %in% rownames(meta.data)]


## celltypes counts
temp <- as.data.frame(table(meta.data$celltype))
pdf(paste0(files.path,'barplotcelltypes_Randolph2021_PBMC_Clean.pdf'))
ggplot(temp,aes(x=Var1,y=Freq)) + geom_bar(stat = 'identity') +
labs(x='Cell types',y='Number of cells') +
theme(
    axis.title.x = element_text(size = 16,face="bold"),
    axis.title.y = element_text(size = 16,face="bold"),
    axis.text.x = element_text(face="bold", size=14, angle = 45, vjust = 0.5, hjust=0.5),
    axis.text.y = element_text(face="bold", size=14,) )
dev.off()

rm(temp)

## Normalization: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
# This is then natural-log transformed using log1p
mrnaNorm <- singlecellmethods::normalizeData(ge_counts, method = "log")


## Variable gene selection using VST to select most variable genes 
## on 'raw counts' instead of normalized count
var_genes <- singlecellmethods::vargenes_vst(object = ge_counts, topn = 3000) #

## SCALE DATA ONLY ON THE MOST EXPRESSED GENES. Adding Consine normalization
mrnaCosScaled <- mrnaNorm[rownames(mrnaNorm) %in% var_genes, ] %>% scaleData() %>% cosine_normalize(2)

## PCA: 20
pca_mrnaScaled <- irlba::prcomp_irlba(t(mrnaCosScaled ), 20)

## HARMONYZATION (Harmony): (batchID, SOC_indiv_ID)
harmony <- HarmonyMatrix(pca_mrnaScaled$x, meta.data, c("SOC_indiv_ID","batchID"), theta = c(1,1), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = 100, max.iter.harmony = 20,
                               max.iter.cluster = 20, do_pca = F, verbose = T)

colnames(harmony) <- paste0("harmonized_", colnames(harmony), sep="")
meta.data <- cbind(meta.data, harmony)

## UMAP: top 20 PCs using (PCA/Harmony_PCs) 
umap_mrnaScaled <- uwot::umap(harmony, n_neighbors = 30, metric = "euclidean", min_dist = .1)

meta.data$harmonized_UMAP1 <- umap_mrnaScaled[, 1]
meta.data$harmonized_UMAP2 <- umap_mrnaScaled[, 2]


## Save
save(meta.data,mrnaNorm,ge_counts,pca_mrnaScaled, file=paste0(files.path,'Randolph2021_PBMC_Clean.Rda'))



## Reproducibility info
print('Reproducibility information:')
Sys.time()
proc.time()
devtools::session_info()