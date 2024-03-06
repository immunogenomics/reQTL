# TopGenes_geneSetEnrichment_Randolph2021_IAV.R

## PLOT NORMALIZED SINGLE CELL GENE EXPRESSION VS SINGLE CELL PERTURBATION SCORE

# conda activate newENv

## TO RUN: 
# Rscript TopGenes_geneSetEnrichment_Randolph2021_IAV.R > TopGenes_geneSetEnrichment_Randolph2021_IAV_log.txt 2>&1

## Load libraries
library(plyr) # always first and then dplyr
library(dplyr)
library(tidyr)
library(singlecellmethods)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(ggrastr)
library(grid)
library(gridExtra)
library(ggpattern)
library(scales)
theme_set(theme_classic())
options(stringsAsFactors=FALSE)
source('~/10x_Multiome_GEaTac/10x_PBMC/Rfunctions.R')

## Create directory for ANALYSIS
main_dir <- '/path_to_directory'
sub_dir <- 'geneSetEnrichment'
output.dir <- file.path(main_dir,sub_dir)
if(!dir.exists(output.dir)){
dir.create(output.dir) 
} else{
	print('Directory already exist')
}

files.path <- paste0(output.dir,'/')
files.data <- '/path_to_directory/'


## SIGNIFICANT CORRELATION OF MOST VARIABLE GENES PBMCS
corrs.df <- readRDS(paste0(files.path,'corr_PerturbationScore_EXP_Randolph2021_IAV.rds'))

## ORDER BY CORR
corrs.df <- corrs.df[order(corrs.df$corrs,decreasing=TRUE), ]
select_genes <- corrs.df[1:20,'gene']


## OPEN DATA
load("/path_to_directory/Randolph2021_PBMC_Clean.Rda") #
rm(pca_mrnaScaled,ge_counts);gc() #  Keep ge_counts mrnaNorm  "meta.data"

## CELLTYPES:
celltypes <- c('B', 'CD4_T', 'CD8_T', 'monocytes', 'NK')

## FILTER DATA 202035 cells
## Select cell from ge_counts, meta.data
meta.data <- meta.data[meta.data$celltype %in% celltypes, ]
meta.data <- meta.data[meta.data$SOC_indiv_ID != "HMN52545",] # 199,879 cells 


## OPEN PERTURBATION SCORE
lda <- readRDS(paste0(files.data,'PerturbationScore_Randolph2021_IAV.rds')) #

meta.data$SCORE <- lda[row.names(meta.data),'s0']
rm(lda); gc()

## SCALE DATA ONLY ON THE MOST VARIABLE GENES. Adding Consine normalization
mrnaNorm <- mrnaNorm[, row.names(meta.data)]

norm_exp <- mrnaNorm[select_genes, ] %>% as.matrix() %>% t() %>% as.data.frame()
rm(mrnaNorm); gc()

## MERGE
tmp <- meta.data[,c('SCORE','SOC_indiv_ID','SOC_infection_status','celltype')]
tmp <- merge(tmp,norm_exp,by='row.names')

## PLOT TOP 10
list_genes <- as.list(select_genes[1:10])

## PLOTS
plots <- lapply(list_genes,function(i) {
	gene <- i 
	lowx <- floor(min(tmp$SCORE))
	highx <- ceiling(max(tmp$SCORE))
	corr <- cor(tmp$SCORE,tmp[,gene] ) %>% formatC(digits=2)
	plot_dens <- get_density(tmp$SCORE, tmp[,gene],n = 100 ) 
	scale.min <- 0
	scale.max <- formatC(max(plot_dens),digits=2) %>% as.numeric()
	#
	ggplot(data=tmp, aes(x = SCORE, y =tmp[,gene], fill = plot_dens )) + 
  	rasterize(geom_point(size = 0.75, shape = 21, stroke = 0.00001),dpi=300 ) + 
	scale_fill_viridis_c(rescaler = function(x, to = c(0, 1), from = NULL) {
    	ifelse(x<0.20, 
           scales::rescale(x,
                           to = to,
                           from = c(min(x, na.rm = TRUE), 0.20)),
           1)}, direction=1 ) +
  	annotate("text", x = Inf, y = Inf,hjust = 1, vjust = 1 + 0.15 * 3, label =paste('r=',corr ),size=6 ) +
  	scale_x_continuous(limits=c(lowx ,highx),breaks=scales::pretty_breaks(n=5) ) +
  	labs(title=gene, x='Perturbation score',y='Normalize expression', fill = "density" ) + 
  	theme(legend.position = "none",
    	plot.title = element_text(size=14,face="bold",hjust=0.5),
    	axis.title = element_text(size = 14,face="bold"),
    	axis.text = element_text(size = 12,face="bold"),
    	legend.key.size = unit(0.5, "cm") )
	})
## Produce this list of plot 10 genes each 3.5x4
# we want 5 plots per row
nplots <- length(plots)
nrowp <- ceiling(nplots/5) # number of rows 
wp=3.5*5
hp=3.5*nrowp
#
pdf(paste0(files.path,'TopGenes_geneSetEnrichment_Randolph2021_IAV.pdf'),w=wp,h=hp)
do.call('grid.arrange',c(plots,ncol=5,nrow=nrowp)) 
dev.off()

## clean
rm(list=ls())



## Reproducibility info
print('Reproducibility Information:')
Sys.time()
proc.time()
devtools::session_info()