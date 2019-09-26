library(reshape)
library(ggpubr)
library(ggplot2)
library(gplots)
library(ComplexHeatmap)
library(matrixStats)
library(TCGA2STAT)
library(circlize)

library(ggfortify)
library(survival)
library(survminer)
library(reshape)
library(gplots)
library(ComplexHeatmap)

setwd("/Users/devlij03/Google Drive/Png/cooper/PNAS_manuscript/int/glioma")

source("/Users/devlij03/Google Drive/Png/cooper/ciber/ciber_source/CIBERSORT.R") ## Source code must be acquired from CIBERSORT Developers

source("ciber_plotter2.R")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cell_type_cols <- gg_color_hue(4)

better_col_vector <- c(Macs_IFN_B1a_Macs_IFN_y="lightcoral",Macs_IL_10="red2",
                       Macs_IL_4_Macs_IL_13 = "red4", 
                       DCs_IFN_B1a="lightskyblue",
                       DCs_IFN_B1a_DCs_IFN_y="dodgerblue2", DCs_IL_10="navy", 
                       Monos_IFN_B1a_Monos_IFN_y="lightgreen", Monos_IL_10="green3",
                       Monos_IL_4_Monos_IL_13="darkgreen", 
                       PMNs_IFN_B1a="plum2",
                       PMNs_IFN_B1a_PMNs_IFN_y="mediumorchid2",
                       PMNs_IL_4_PMNs_IL_13="purple4")




##### Step 1 load in all cancer data, but only save genes in basis matrix
##### compile rnaseq rsem values and metadata

#
sig_matrix <- "curated_basis12.txt"
sig_matrix2 <- "LM22.txt"
sig_matrix3 <- "immunoStates_basis.txt"
X <- read.table(sig_matrix,header=T,sep="\t")
X2 <- read.table(sig_matrix2,header=T,sep="\t")
X3 <- read.table(sig_matrix3,header=T,sep="\t")
basis <- unique(c(as.character(X$GeneSymbol), as.character(X2$GeneSymbol), as.character(X3$GeneSymbol)))

#pancancers
cancer_key <- read.table("pan_cancer_key.txt", header=T, sep='\t')


c_type <- "Glioma"
c_acronym <- "GBMLGG"

rsem.value <- getTCGA(disease = c_acronym, data.type = "RNASeq2",
                      clinical = T)
rsem_df <- rsem.value$dat
colnames(rsem_df) <- substring(colnames(rsem_df), 0,12)

meta <- data.frame(rsem.value$clinical)

keepers <- intersect(colnames(rsem_df), rownames(meta))
#

meta <- meta[keepers,]
meta$id <- rownames(meta)
meta$survival <- as.numeric(as.character(meta$daystodeath))

###

rsem_df2 <- cbind(GeneSymbol = rownames(rsem_df), rsem_df[,keepers])

rsem_df_keep <- rsem_df2[rsem_df2[,1] %in% basis,]
gnamer <- paste0("genesymbols/", c_type, "_genesymbols.txt")
write.table(rsem_df_keep, gnamer, row.names = F, quote = F, sep='\t')

#meta$age = as.numeric(as.character(meta$yearstobirth))


#print(colnames(meta))
###add to meta IDH mutation status

mut_df <- read.table("idh_status.txt", T, '\t')
rownames(mut_df) <- as.character(mut_df$Case)
mut_df <- mut_df[keepers,]

### add aditional metadata
meta <- cbind(meta, mut_df)

mnamer <- paste0("metadata/", c_type, "_metadata.txt")
write.table(meta, mnamer, row.names = F, quote = F, sep='\t')


##ciber

pan_results <- CIBERSORT(sig_matrix, gnamer, perm=100, QN=T, 
                         absolute=F, abs_method='sig_score')
rownames(pan_results[[1]]) <- gsub("\\.", "-", rownames(pan_results[[1]]))
results_write <- data.frame(id = rownames(pan_results[[1]]), pan_results[[1]])
cnamer <- paste0("ciber_ours/", c_type, "_results_ours.txt")
write.table(results_write, cnamer, row.names=F, sep='\t', quote = F)


#
#
#
#

#### survival curves
cibers <- list.files("ciber_ours", full.names = T)
metas <- list.files(path = 'metadata', full.names = T)

cts <- colnames(X)[-1]
full_list <- list()
tres <- data.frame(high_mean=NA, low_mean=NA, ct=NA, cancer=NA, pvalue=NA)
counter = 1
for (k in 1:length(cibers)){
  ciber_curr <- read.table(cibers[k], header=T, sep='\t')
  meta_curr <- read.table(metas[k], header=T, sep='\t')
  
  df <- merge(ciber_curr, meta_curr, by = 'id')
  for (i in 1:nrow(df)){
    if (is.na(df$survival[i])){
      df$survival[i] <- df$daystolastfollowup[i]
    }
  }
  
  cts2do <- c()
  for (j in 1:length(cts)){
    nums <- df[,cts[j]]
    
    if (any(nums > 0)){
      
      q3 <- summary(nums)[5]
      q1 <- summary(nums)[2]
      
      df$status <- rep(0, nrow(df))
      df$status[nums > q3] <- 1
      df$status[nums <= q1] <- -1
      #print(table(df$status))
      
      df_keep <- subset(df, status != 0)
      #print(nrow(df_keep))
      if (nrow(df_keep)>5){
        if(table(df$status)[2] > 2){
          cts2do <- c(cts2do, cts[j])
        }
      }
    }
  }
  glist = list()
  for (f in 1:length(cts2do)){
    nums <- df[,cts2do[f]]
    
    q3 <- summary(nums)[5]
    q1 <- summary(nums)[2]
    
    df$status <- rep(0, nrow(df))
    df$status[nums > q3] <- 1
    df$status[nums <= q1] <- -1
    #print(table(df$status))
    
    df_keep <- subset(df, status != 0)
    df_keep$survival <- df_keep$survival/365
    
    
    survplot <- try(do.call(survfit, 
                            list(formula = Surv(df_keep$survival, df_keep$vitalstatus==1) ~ df_keep$status,
                                 conf.type = "log-log", data = df_keep)))
    pp <- try(ggsurvplot(survplot, conf.int = TRUE, legend.labs=c("low", "high"),
                         ggtheme = theme_minimal(), pval = T, risk.table = F, risk.table.height =0.4,
                         title= paste0(cts2do[f])))
    
    survy <- survdiff(Surv(df_keep$survival, df_keep$vitalstatus==1) ~ df_keep$status)
    pval <- round(1 - pchisq(survy$chisq, 1), 3)
    
    df_keep$status <- factor(df_keep$status, levels = c(-1,1))
    huh <- t.test(df_keep$survival ~ df_keep$status)
    low_mean = huh$estimate[1]
    high_mean = huh$estimate[2]
    
    cts2 <- gsub("\\..*", "", basename(cibers[k]))
    adder <- data.frame(high_mean=high_mean, low_mean=low_mean, 
                        ct=cts2do[f], cancer=cts2, pvalue = pval)
    tres <- rbind(tres, adder)
    
    glist[[f]] <- pp
    full_list[[counter]] <- pp
    counter <- counter+1
  }
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves.pdf")
  #pdf(namer, height = 15, width = 20)
  what <- arrange_ggsurvplots(glist, print = F, nrow = 4, ncol = 3)
  ggsave(namer, what, height = 11, width = 11)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves.RData")
  save(glist, file=namer)
  
}
  

#
#
#
#
##
##

#
#### survival curves for diff cancer types

cibers <- list.files("ciber_ours", full.names = T)
metas <- list.files(path = 'metadata', full.names = T)

cts <- colnames(X)[-1]
full_list <- list()
tres <- data.frame(high_mean=NA, low_mean=NA, ct=NA, cancer=NA, pvalue=NA)
counter = 1
for (k in 1:length(cibers)){
  ciber_curr <- read.table(cibers[k], header=T, sep='\t')
  meta_curr <- read.table(metas[k], header=T, sep='\t')
  
  df <- merge(ciber_curr, meta_curr, by = 'id')
  for (i in 1:nrow(df)){
    if (is.na(df$survival[i])){
      df$survival[i] <- df$daystolastfollowup[i]
    }
  }
  
  cts2do <- c()
  for (j in 1:length(cts)){
    nums <- df[,cts[j]]
    
    if (any(nums > 0)){
      
      q3 <- summary(nums)[5]
      q1 <- summary(nums)[2]
      
      df$status <- rep(0, nrow(df))
      df$status[nums > q3] <- 1
      df$status[nums <= q1] <- -1
      #print(table(df$status))
      
      df_keep <- subset(df, status != 0)
      #print(nrow(df_keep))
      if (nrow(df_keep)>5){
        if(table(df$status)[2] > 2){
          cts2do <- c(cts2do, cts[j])
        }
      }
    }
  }
  wtlist = list()
  mutlist = list()
  glist = list()
  for (f in 1:length(cts2do)){
    nums <- df[,cts2do[f]]
    
    q3 <- summary(nums)[5]
    q1 <- summary(nums)[2]
    
    df$status <- rep(0, nrow(df))
    df$status[nums > q3] <- "high"
    df$status[nums <= q1] <- "low"
    #print(table(df$status))
    
    df_keep <- subset(df, status != 0)
    df_keep$survival <- df_keep$survival/365
    
    
    
    survplot <- try(do.call(survfit, 
                            list(formula = Surv(df_keep$survival, df_keep$vitalstatus==1) ~ df_keep$status+df_keep$IDH.status,
                                 conf.type = "log-log", data = df_keep)))
    pp <- try(ggsurvplot(survplot, conf.int = TRUE, xlim=c(0,5), break.x.by=1, legend='none',
                         ggtheme = theme_minimal(), pval = T, risk.table = T, risk.table.height =0.4,
                         title= paste0(cts2do[f])))
    
    glist[[f]] <- pp
    
    df_mut <- subset(df_keep, IDH.status == "Mutant")
    
    survplot <- try(do.call(survfit, 
                            list(formula = Surv(df_mut$survival, df_mut$vitalstatus==1) ~ df_mut$status+df_mut$IDH.status,
                                 conf.type = "log-log", data = df_mut)))
    pp <- try(ggsurvplot(survplot, conf.int = TRUE, xlim=c(0,5), break.x.by=1, legend='none',
                         ggtheme = theme_minimal(), pval = T, risk.table = T, risk.table.height =0.4,
                         title= paste0(cts2do[f])))
    
    mutlist[[f]] <- pp
    
    df_wt <- subset(df_keep, IDH.status == "WT")
    
    survplot <- try(do.call(survfit, 
                            list(formula = Surv(df_wt$survival, df_wt$vitalstatus==1) ~ df_wt$status+df_wt$IDH.status,
                                 conf.type = "log-log", data = df_wt)))
    pp <- try(ggsurvplot(survplot, conf.int = TRUE, xlim=c(0,5), break.x.by=1, legend='none',
                         ggtheme = theme_minimal(), pval = T, risk.table = T, risk.table.height =0.4,
                         title= paste0(cts2do[f])))
    
    wtlist[[f]] <- pp
  }
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idh.pdf")
  #pdf(namer, height = 15, width = 20)
  what <- arrange_ggsurvplots(glist, print = F, nrow = 4, ncol = 3)
  ggsave(namer, what, height = 20, width = 20)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idh.RData")
  save(glist, file=namer)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idhWT.pdf")
  #pdf(namer, height = 15, width = 20)
  what <- arrange_ggsurvplots(wtlist, print = F, nrow = 4, ncol = 3)
  ggsave(namer, what, height = 20, width = 20)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idhWT.RData")
  save(wtlist, file=namer)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idhMUT.pdf")
  #pdf(namer, height = 15, width = 20)
  what <- arrange_ggsurvplots(mutlist, print = F, nrow = 4, ncol = 3)
  ggsave(namer, what, height = 20, width = 20)
  
  namer <- paste0("survival/", gsub("\\..*", "", basename(cibers[k])), "_survival_curves_idhMUT.RData")
  save(mutlist, file=namer)
  
}
#
#
#
#
#
#

#
#
#

#
#
#
#
#### gene based survival curves



### load precog subset
precog_df_gse <- read.table("ind_test/precog_GSE16011_genesymbols.txt", header=T, sep='\t')
precog_meta_gse <- read.table("ind_test/GSE16011.info.tsv", header=T, sep='\t')
rownames(precog_meta_gse) <- precog_meta_gse$Array
precog_meta_gse <- subset(precog_meta_gse, !is.na(OS_Status))
alright <- intersect(precog_meta_gse$Array, colnames(precog_df_gse))
precog_filt_gse <- precog_df_gse[,alright]
rownames(precog_filt_gse) <- precog_df_gse$Description
precog_meta_gse <- precog_meta_gse[alright,]

table(precog_meta_gse$OS_Status)
## survival curve gene by gene cox regression
library(survminer)
library(survival)
library(cowplot)

pp <- list()
for (i in 1:nrow(precog_filt_gse)){
  gener <- rownames(precog_filt_gse)[i]
  precog_meta_gse$gExpr <- log2(unlist(precog_filt_gse[i,]))
  df <- precog_meta_gse
  nums <- precog_meta_gse$gExpr
  
  q3 <- summary(nums)[5]
  q1 <- summary(nums)[2]
  med <- summary(nums)[3]
  
  sd(nums)
  
  #meds <- median(lung_df$ct_ratio[is.finite(lung_df$ct_ratio)], na.rm = T)
  #high_n <- length(which(df$ct_ratio > meds))
  #low_n <- nrow(df) - high_n
  
  df$status <- rep(0, nrow(df))
  df$status[nums > q3] <- 1
  df$status[nums <= q1] <- -1
  #print(table(df$status))
  
  df_keep <- subset(df, status != 0)
  #surv_time <- Surv(df_keep$survival, df_keep$vitalstatus==1)
  #survplot <- survfit(surv_time~df_keep$status,conf.type = "log-log", data = df_keep)
  
  survplot <- try(do.call(survfit, 
                          list(formula = Surv(df_keep$OS_Time, df_keep$OS_Status==1) ~ df_keep$status,
                               conf.type = "log-log", data = df_keep)))
  pp[[i]] <- try(ggsurvplot(survplot, conf.int = TRUE, legend.labs=c("low", "high"),
                            ggtheme = theme_minimal(), pval = T, risk.table = T, risk.table.height =0.4,
                            title= gener))
}

pdf("ind_test/surv_test_GSE_quart.pdf", height = 25, width = 20)
arrange_ggsurvplots(pp, nrow = 6, ncol = 5)
dev.off()

save(pp, file="ind_test/surv_test_GSE_quart.RData")


#
#
#
# TCGA

### load TCGA and metadata
tcga_df <- read.table("genesymbols/Glioma_genesymbols.txt", header=T, sep='\t')
tcga_meta <- read.table("metadata/Glioma_metadata.txt", header=T, sep='\t')

rownames(tcga_df) <- tcga_df$GeneSymbol
tcga_df <- tcga_df[,-1]
colnames(tcga_df) <- gsub("\\.", "-", colnames(tcga_df))
cur_genes <- intersect(read.table("curated_basis12.txt", T, '\t')$GeneSymbol,rownames(tcga_df))
alright <- intersect(cur_genes, rownames(tcga_df))
tcga_df <- tcga_df[alright,]

for (i in 1:nrow(tcga_meta)){
  if (is.na(tcga_meta$survival[i])){
    tcga_meta$survival[i] <- tcga_meta$daystolastfollowup[i]
  }
}
tcga_meta <- subset(tcga_meta, !is.na(survival))
alright <- intersect(tcga_meta$id, colnames(tcga_df))
tcga_df_filt <- tcga_df[,alright]
rownames(tcga_meta) <- tcga_meta$id
tcga_meta <- tcga_meta[alright,]

## survival curve gene by gene cox regression
library(survminer)
library(survival)
library(cowplot)

pp <- list()
for (i in 1:nrow(tcga_df_filt)){
  gener <- rownames(tcga_df_filt)[i]
  tcga_meta$gExpr <- log2(unlist(tcga_df_filt[i,]))
  df <- tcga_meta
  nums <- tcga_meta$gExpr
  
  q3 <- summary(nums)[5]
  q1 <- summary(nums)[2]
  med <- summary(nums)[3]
  
  sd(nums)
  
  #meds <- median(lung_df$ct_ratio[is.finite(lung_df$ct_ratio)], na.rm = T)
  #high_n <- length(which(df$ct_ratio > meds))
  #low_n <- nrow(df) - high_n
  
  df$status <- rep(0, nrow(df))
  df$status[nums > q3] <- 1
  df$status[nums <= q1] <- -1
  #print(table(df$status))
  
  df_keep <- subset(df, status != 0)
  #surv_time <- Surv(df_keep$survival, df_keep$vitalstatus==1)
  #survplot <- survfit(surv_time~df_keep$status,conf.type = "log-log", data = df_keep)
  df_keep$survival <- df_keep$survival/365
  
  survplot <- try(do.call(survfit, 
                          list(formula = Surv(df_keep$survival, df_keep$vitalstatus) ~ df_keep$status,
                               conf.type = "log-log", data = df_keep)))
  pp[[i]] <- try(ggsurvplot(survplot, conf.int = TRUE, legend.labs=c("low", "high"),
                            ggtheme = theme_minimal(), pval = T, risk.table = T, risk.table.height =0.4,
                            title= gener))
}

pdf("ind_test/surv_test_TCGA_quart.pdf", height = 30, width = 20)
arrange_ggsurvplots(pp, nrow = 8, ncol = 4)
dev.off()

save(pp, file="ind_test/surv_test_TCGA_quart.RData")

