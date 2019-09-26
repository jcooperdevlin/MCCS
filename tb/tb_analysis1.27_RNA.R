###### Full TB Analysis for Manuscript 1.26

##### 

#setwd("~/Google Drive/Png/cooper/Eblood_Submission/archive/tb_1.25/")

library(ggplot2)
library(reshape2)
library("illuminaHumanv4.db")
library("GEOquery")
source("CIBERSORT.R")
#source("ciber_plotter2.R")

args = commandArgs(trailingOnly=TRUE)

##GET full list of necessary genes
bases <- list.files("basis/", full.names = T)
test <- c("curated", "immunoStates", "LM22")
b_genes <- c()
for (i in 1:length(bases)){
  b <- read.table(bases[i], header=T, sep='\t')
  b_genes <- c(b_genes, as.character(b$GeneSymbol))
}

basis_genes <- unique(b_genes)


#
##
#
#

###### Step 1

##### Read in all Microarray datasets and convert to gene symbols
##### Also store metadata for later



#ma_sets <- list.files("raw", full.names = T)
ma_sets <- args[1]
print(ma_sets)
gse_names <- gsub(".*raw/", "", ma_sets)
gse_names <- gsub("_series_.*", "", gse_names)

for (i in 1:length(ma_sets)){
  curr <- getGEO(filename = ma_sets[i])
  expr <- curr@assayData$exprs
  col_num <- grep("symbol", colnames(curr@featureData@data), ignore.case=T)
  expr <- aggregate(expr, by = list(curr@featureData@data[,col_num]), "mean")
  colnames(expr)[1] <- "GeneSymbol"
  expr_keep <- expr[expr[,1] %in% basis_genes,]
  
  enamer <- paste0("genesymbols/",gse_names[i], "_genesymbols.txt")
  write.table(expr_keep, enamer, sep='\t', quote=F, row.names=F)
  
  meta <- curr@phenoData@data
  mnamer <- paste0("metadata/", gse_names[i], "_metadata.txt")
  write.table(meta, mnamer, sep='\t', quote=F, row.names=F)
  
  for (j in 1:length(bases)){
  tb_results <- CIBERSORT(bases[j], enamer, 
                          perm=100, QN=T, 
                          absolute=F, abs_method='sig_score')
  res <- data.frame(tb_results[[1]])
  res$sample <- rownames(res)
  
  namer <- paste0("ciber/", gse_names[i], "_", test[j], ".txt")
  write.table(res, file=namer, sep='\t', row.names=F, quote = F)
  }
}

