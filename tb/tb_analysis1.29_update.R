###### Full TB Analysis for Manuscript 1.26

##### 

setwd("~/Google Drive/Png/cooper/PNAS_manuscript/int/tb/")

library(ggplot2)
library(reshape2)
library(cowplot)
library("illuminaHumanv4.db")
library("GEOquery")
library(randomcoloR)
source("functions/microarray2gene.R")

source("/Users/devlij03/Google Drive/Png/cooper/ciber/ciber_source/CIBERSORT.R") ## Source code must be acquired from CIBERSORT Developers

source("ciber_plotter2.R")



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

###
##### Read in all Microarray datasets and convert to gene symbols
##### Also store metadata for later and run ciber

### ran on cluster

#ma_sets <- list.files("raw_ma", full.names = T)
#gse_names <- gsub("raw_ma/", "", ma_sets)
#gse_names <- gsub("_series_.*", "", gse_names)

# for (i in 1:length(ma_sets)){
#   curr <- getGEO(filename = ma_sets[i])
#   expr <- curr@assayData$exprs
#   col_num <- grep("symbol", colnames(curr@featureData@data), ignore.case=T)
#   expr <- aggregate(expr, by = list(curr@featureData@data[,col_num]), "mean")
#   colnames(expr)[1] <- "GeneSymbol"
#   expr_keep <- expr[expr[,1] %in% basis_genes,]
#   
#   enamer <- paste0("genesymbols/",gse_names[i], "_genesymbols.txt")
#   write.table(expr_keep, enamer, sep='\t', quote=F, row.names=F)
#   
#   meta <- curr@phenoData@data
#   mnamer <- paste0("metadata/", gse_names[i], "_metadata.txt")
#   write.table(meta, mnamer, sep='\t', quote=F, row.names=F)
#   
#   for (j in 1:length(bases)){
#     tb_results <- CIBERSORT(bases[j], enamer, 
#                             perm=100, QN=T, 
#                             absolute=F, abs_method='sig_score')
#     res <- data.frame(tb_results[[1]])
#     res$sample <- rownames(res)
#     
#     namer <- paste0("ciber/", gse_names[i], "_", test[j], ".txt")
#     write.table(res, file=namer, sep='\t', row.names=F, quote = F)
#   }
# }

###### Step 2

##### Read in all RNA datasets and convert to gene symbols
##### Also store metadata for later



# ma_sets <- list.files("raw_RNA", pattern = 'series', full.names = T)
# #ma_sets <- args[1]
# print(ma_sets)
# gse_names <- gsub(".*raw_RNA/", "", ma_sets)
# gse_names <- gsub("_series_.*", "", gse_names)
# 
# for (i in 1:length(ma_sets)){
#   curr <- getGEO(filename = ma_sets[i])
#   
#   if(i==1){
#     e_reader <- list.files("raw_RNA", pattern="xlsx", full.names=T)
#     e_use <- e_reader[grepl(gse_names[i], e_reader)]
#     expr <- data.frame(read_excel(e_use))
#     colnames(expr)[1] <- "GeneSymbol"
#     expr_keep <- expr[expr[,1] %in% basis_genes,]
#   }else if (i==2) {
#     expr_keep <- data.frame(GeneSymbol = basis_genes)
#     for(j in 2:5){
#       e_reader <- list.files("raw_RNA", pattern="xlsx", full.names=T)
#       e_use <- e_reader[j]
#       expr <- data.frame(read_excel(e_use))
#       expr <- expr[,c(2, 4:ncol(expr))]
#       colnames(expr)[1] <- "GeneSymbol"
#       expr_keeper <- expr[expr[,1] %in% basis_genes,]
#       rownames(expr_keeper) <- expr_keeper[,1]
#       expr_keeper <- expr_keeper[basis_genes,]
#       expr_keeper <- expr_keeper[,-1]
#       expr_keep <- cbind(expr_keep, expr_keeper)
#     }
#   } else if(i==3){
#     e_reader <- list.files("raw_RNA", pattern="xlsx", full.names=T)
#     e_use <- e_reader[grepl(gse_names[i], e_reader)]
#     expr <- data.frame(read_excel(e_use))
#     expr <- expr[,6:ncol(expr)]
#     colnames(expr)[1] <- "GeneSymbol"
#     expr_keep <- expr[expr[,1] %in% basis_genes,]
#   }  else if (i ==4){
#     e_reader <- list.files("raw_RNA", pattern="csv", full.names=T)
#     e_use <- e_reader[grepl(gse_names[i], e_reader)]
#     expr <- read.table(e_use, header=T, sep=',')
#     expr <- expr[,-1]
#     colnames(expr)[1] <- "GeneSymbol"
#     expr_keep <- expr[expr[,1] %in% basis_genes,]
#   } else if (i == 5){
#     e_reader <- list.files("raw_RNA", pattern="csv", full.names=T)
#     e_use <- e_reader[grepl(gse_names[i], e_reader)]
#     expr <- read.table(e_use, header=T, sep=',')
#     expr <- expr[,-1]
#     colnames(expr)[1] <- "GeneSymbol"
#     expr_keep <- expr[expr[,1] %in% basis_genes,]
#     
#   }
#   
#   expr_keep[is.na(expr_keep)] <- 0
#   
#   enamer <- paste0("genesymbols/",gse_names[i], "_genesymbols.txt")
#   write.table(expr_keep, enamer, sep='\t', quote=F, row.names=F)
#   
#   meta <- curr@phenoData@data
#   mnamer <- paste0("metadata/", gse_names[i], "_metadata.txt")
#   write.table(meta, mnamer, sep='\t', quote=F, row.names=F)
#   
#   
#   for (j in 1:length(bases)){
#     tb_results <- CIBERSORT(bases[j], enamer, 
#                             perm=100, QN=T, 
#                             absolute=F, abs_method='sig_score')
#     res <- data.frame(tb_results[[1]])
#     res$sample <- rownames(res)
#     
#     namer <- paste0("ciber/", gse_names[i], "_", test[j], ".txt")
#     write.table(res, file=namer, sep='\t', row.names=F, quote = F)
#   }
# }



##### read in results and think about plots
meta_specs <- list(
  GSE101705 = c("condition: latent TB infection", "condition: TB"),
  GSE107995 = c("group: Active_TB", "group: Control", "group: LTBI"),
  GSE19491 = c("Whole blood from healthy control", 
        "Whole Blood from healthy control", 
        "Whole Blood from patient with Latent TB", 
        "Whole Blood from patient with Active TB"),
  GSE28623 = c("Peripheral blood LTBI", "Peripheral blood NID","Peripheral blood TB"),
  GSE37250 = c("Whole Blood from patient with LTBI", "Whole Blood from patient with Active TB"),
  GSE39939 = c("Whole Blood from patient with LTBI", "Whole Blood from patient with Active TB"),
  GSE39940 = c("Whole Blood from patient with LTBI", "Whole Blood from patient with Active TB"),
  GSE40553 = c("disease: LTB", "disease: PTB"), 
  GSE41055 = c("whole blood, healthy control", "whole blood, latent TB infection", "whole blood, active TB infection"),
  GSE56153 = c("condition: Active","condition: Control","condition: Recover","condition: Treatment"),
  GSE79362 = c("group: case (TB progressor)", "group: control (non-progressor)"),
  GSE89403 = c("disease state: Healthy Controls", "disease state: TB Subjects"),
  GSE94438 = c("group: case (TB)", "group: Control")
)

id_cats <- c("characteristics_ch1", "characteristics_ch1.1", 
             "source_name_ch1", "source_name_ch1", "source_name_ch1", "source_name_ch1", "source_name_ch1",
             "characteristics_ch1", "source_name_ch1", "characteristics_ch1.2", "characteristics_ch1",
             "characteristics_ch1.3", "characteristics_ch1.6")

meta_plotter <- data.frame(sample=NA,condition=NA, region=NA, study = NA, data_type=NA)
cur_plotter <- read.table("ciber/GSE39940_curated.txt", header=T, sep='\t')[1,]
immuno_plotter <- read.table("ciber/GSE39940_immunoStates.txt", header=T, sep='\t')[1,]
lm22_plotter <- read.table("ciber/GSE39940_LM22.txt", header=T, sep='\t')[1,]

for (i in c(1:length(gse_names))){
  looker <- list.files("ciber", pattern = gse_names[i], full.names = T)
  glooker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  meta <- read.table(
    list.files("metadata", pattern = gse_names[i], full.names = T),
    header=T, sep='\t', quote="", fill=T)
  id_cat_num <- grep("geo_accession", colnames(meta))
  meta <- subset(meta, meta[,id_cats[i]] %in% meta_specs[[i]])
  
  if(i ==1){
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*healthy.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Latent.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*TB.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("geographical region: South India", nrow(meta))
    data_type_add <- rep("RNASeq", nrow(meta))
  }else if (i ==2){
    sub_check = meta$characteristics_ch1.12
    sub_check[sub_check == ""] <- "timepoint_months: Baseline"
    meta = subset(meta, sub_check == "timepoint_months: Baseline")
    meta <- meta[which(duplicated(meta$characteristics_ch1.3)==F),]
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*Control.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*ltbi.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- meta$title
    region_add <- gsub(".*London.*", "UK", region_add, ignore.case = T)
    region_add <- gsub(".*SouthAfrica.*", "South Africa", region_add, ignore.case = T)
    region_add <- gsub(".*Leicester.*", "UK", region_add, ignore.case = T)
    data_type_add <- rep("RNASeq", nrow(meta))
  }else if (i==3){
    meta[,id_cats[i]] <- gsub("Blood", "blood", meta[,id_cats[i]])
    meta <- meta[which(grepl("long", meta[,1])==F),]
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*healthy.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Latent.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- meta$characteristics_ch1.4
    region_add <- gsub(".*London.*", "UK", region_add, ignore.case = T)
    region_add <- gsub(".*South Africa.*", "South Africa", region_add, ignore.case = T)
    data_type_add <- rep("Microarray", nrow(meta))
    
  } else if (i==4) {
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*NID.*", "Healthy", thing)
    thing <- gsub(".*LTBI.*", "Latent", thing)
    thing <- gsub(".*TB.*", "Active", thing)
    meta$condition <- thing
    region_add <- rep("West Africa", nrow(meta))
    data_type_add <- rep("Microarray", nrow(meta))
  } else if (i == 5|i==6|i==7) {
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*LTBI.*", "Latent", thing)
    thing <- gsub(".*Active.*", "Active", thing)
    meta$condition <- thing
    region_add <- meta$characteristics_ch1.2
    region_add <- gsub("geographical region: ", "", region_add)
    data_type_add <- rep("Microarray", nrow(meta))
  } else if (i==8){
    sub_check = gsub(".*\\[", "", meta$title)
    sub_check = gsub("-.*", "", sub_check)
    sub_check = gsub("/.*", "", sub_check)
    meta$sub <- sub_check
    
    meta = subset(meta, 
                  characteristics_ch1.4 == "status: untreated latent TB" |
                  characteristics_ch1.4 == "status: active TB pre-treatment")
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*LTB.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*PTB.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- meta$characteristics_ch1.3
    region_add <- gsub(".*UK.*", "UK", region_add, ignore.case = T)
    region_add <- gsub(".*South Africa.*", "South Africa", region_add, ignore.case = T)
    data_type_add <- rep("Microarray", nrow(meta))
  } else if (i==9){
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*healthy.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Latent.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("Venezuela", nrow(meta))
    data_type_add <- rep("Microarray", nrow(meta))
  } else if (i==10){
    meta <- subset(meta,
                   characteristics_ch1.2 == "condition: Active" |
                   characteristics_ch1.2 == "condition: Control")
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*Control.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("Indonesia", nrow(meta))
    data_type_add <- rep("Microarray", nrow(meta))
  } else if (i==11){
    sub_check = gsub(".*\\(", "", meta$title)
    sub_check = gsub("PAXGENE.*", "", sub_check)
    sub_check = gsub("D0.*", "", sub_check)
    sub_check = gsub("_L.*", "", sub_check)
    sub_check = gsub("_", "", sub_check)
    meta$sub <- sub_check
    meta <- meta[which(duplicated(meta$sub)==F),]
    
    cur_add <- read.table(looker[1], header=T, sep='\t')
    cur_names <- gsub("[^[:alnum:] ]", "", cur_add$sample)
    cur_names <- sub('.', '', cur_names)
    
    m_names <- gsub(".*progressor ", "", as.character(meta$title))
    m_names <- gsub("[^[:alnum:] ]", "", m_names)
    rownames(meta) <- m_names
    meta$new_names <- m_names

    keep_keep <- intersect(m_names, cur_names)
    meta <- meta[keep_keep,]
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*Control.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*TB.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("South Africa", nrow(meta))
    data_type_add <- rep("RNASeq", nrow(meta))
  } else if (i==12){
    meta$sub <- meta$characteristics_ch1.2
    meta <- meta[grepl("DX", meta$characteristics_ch1.5),]
    meta <- meta[which(duplicated(meta$sub)==F),]
    meta$new_names <- gsub("sample_code: ", "", meta$characteristics_ch1.1)
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*Control.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*TB.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("South Africa", nrow(meta))
    data_type_add <- rep("RNASeq", nrow(meta))
  } else if (i==13){
    
    meta_c <- subset(meta, characteristics_ch1.6 == "group: Control")
    meta_c <- subset(meta_c, characteristics_ch1.7 == "time.from.exposure.months: 0")
    meta_c <- meta_c[which(duplicated(meta_c$characteristics_ch1.2)==F),]
    
    meta_a <- subset(meta, characteristics_ch1.6 == "group: case (TB)")
    meta_a <- subset(meta_a, characteristics_ch1.7 == "time.from.exposure.months: 0")
    meta_a <- meta_a[which(duplicated(meta_a$characteristics_ch1.2)==F),]
    
    meta <- rbind(meta_c, meta_a)
    meta$new_names <- paste0("X", gsub("code: ", "", meta$characteristics_ch1.1))
    
    thing <- meta[,id_cats[i]]
    thing <- gsub(".*Control.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*TB.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("Africa", nrow(meta))
    data_type_add <- rep("RNASeq", nrow(meta))
  }
  
  if(i %in% c(1,2)){
    #gene_add <- read.table(glooker, header=T, sep='\t')
    #gene_add2 <- aggregate(gene_add[,-1], by=list(gene_add$GeneSymbol), "mean")
    #rownames(gene_add2) <- gene_add2$Group.1
    #gene_add2 <- gene_add2[,as.character(meta$title)]
    #colnames(gene_add2) <- meta[,id_cat_num]
    #gene_add3 <- data.frame(GeneSymbol =rownames(gene_add2), gene_add2)
    #write.table(gene_add3, glooker, row.names=F, quote=F, sep='\t')
    
    cur_add <- read.table(looker[1], header=T, sep='\t')
    rownames(cur_add) <- cur_add$sample
    cur_add <- cur_add[as.character(meta$title),]
    cur_add$sample <- meta[,id_cat_num]
    
    
    immuno_add <- read.table(looker[2], header=T, sep='\t')
    rownames(immuno_add) <- immuno_add$sample
    immuno_add <- immuno_add[as.character(meta$title),]
    immuno_add$sample <- meta[,id_cat_num]
    
    lm22_add <- read.table(looker[3], header=T, sep='\t')
    rownames(lm22_add) <- lm22_add$sample
    lm22_add <- lm22_add[as.character(meta$title),]
    lm22_add$sample <- meta[,id_cat_num]
  } else if (i==11){
    #gene_add <- read.table(glooker, header=T, sep='\t')
    #gene_add <- aggregate(gene_add[,-1], by=list(gene_add$GeneSymbol), "mean")
    #gene_names <- gsub("[^[:alnum:] ]", "", colnames(gene_add))
    #gene_names <- sub('.', '', gene_names)
    #colnames(gene_add) <- gene_names
    #gene_add2 <- gene_add[,as.character(meta$new_names)]
    #colnames(gene_add2) <- meta[,id_cat_num]
    #gene_add3 <- data.frame(GeneSymbol =gene_add$roup1, gene_add2)
    #write.table(gene_add3, glooker, row.names=F, quote=F, sep='\t')
    
    cur_add <- read.table(looker[1], header=T, sep='\t')
    cur_names <- gsub("[^[:alnum:] ]", "", cur_add$sample)
    cur_names <- sub('.', '', cur_names)
    rownames(cur_add) <- cur_names
    cur_add <- cur_add[as.character(meta$new_names),]
    cur_add$sample <- meta[,id_cat_num]
    
    
    immuno_add <- read.table(looker[2], header=T, sep='\t')
    immuno_names <- gsub("[^[:alnum:] ]", "", immuno_add$sample)
    immuno_names <- sub('.', '', immuno_names)
    rownames(immuno_add) <- immuno_names
    immuno_add <- immuno_add[as.character(meta$new_names),]
    immuno_add$sample <- meta[,id_cat_num]
    
    lm22_add <- read.table(looker[3], header=T, sep='\t')
    lm22_names <- gsub("[^[:alnum:] ]", "", lm22_add$sample)
    lm22_names <- sub('.', '', lm22_names)
    rownames(lm22_add) <- lm22_names
    lm22_add <- lm22_add[as.character(meta$new_names),]
    lm22_add$sample <- meta[,id_cat_num]
    
  } else if (i==12|i==13){
    #gene_add <- read.table(glooker, header=T, sep='\t')
    #gene_add2 <- aggregate(gene_add[,-1], by=list(gene_add$GeneSymbol), "mean")
    #rownames(gene_add2) <- gene_add2$Group.1
    #gene_add2 <- gene_add2[,as.character(meta$new_names)]
    #colnames(gene_add2) <- meta[,id_cat_num]
    #gene_add3 <- data.frame(GeneSymbol =rownames(gene_add2), gene_add2)
    #write.table(gene_add3, glooker, row.names=F, quote=F, sep='\t')
    
    
    cur_add <- read.table(looker[1], header=T, sep='\t')
    rownames(cur_add) <- cur_add$sample
    cur_add <- cur_add[as.character(meta$new_names),]
    cur_add$sample <- meta[,id_cat_num]
    
    immuno_add <- read.table(looker[2], header=T, sep='\t')
    rownames(immuno_add) <- immuno_add$sample
    immuno_add <- immuno_add[as.character(meta$new_names),]
    immuno_add$sample <- meta[,id_cat_num]
    
    lm22_add <- read.table(looker[3], header=T, sep='\t')
    rownames(lm22_add) <- lm22_add$sample
    lm22_add <- lm22_add[as.character(meta$new_names),]
    lm22_add$sample <- meta[,id_cat_num]
    
  }else {
    cur_add <- read.table(looker[1], header=T, sep='\t')
    rownames(cur_add) <- cur_add$sample
    cur_add <- cur_add[as.character(meta[,id_cat_num]),]
    cur_add$sample <- meta[,id_cat_num]
    
    immuno_add <- read.table(looker[2], header=T, sep='\t')
    rownames(immuno_add) <- immuno_add$sample
    immuno_add <- immuno_add[as.character(meta[,id_cat_num]),]
    immuno_add$sample <- meta[,id_cat_num]
    
    lm22_add <- read.table(looker[3], header=T, sep='\t')
    rownames(lm22_add) <- lm22_add$sample
    lm22_add <- lm22_add[as.character(meta[,id_cat_num]),]
    lm22_add$sample <- meta[,id_cat_num]
  }

  meta_keep <- as.character(meta[,id_cat_num])
  study_add <- rep(gse_names[i], nrow(meta))
  
  adder <- data.frame(sample = meta[,id_cat_num], condition = meta$condition,
                      region = region_add, study = study_add, data_type = data_type_add)
  
  meta_plotter <- rbind(meta_plotter, adder)
  
  cur_plotter <- rbind(cur_plotter, cur_add)
  immuno_plotter <- rbind(immuno_plotter, immuno_add)
  lm22_plotter <- rbind(lm22_plotter, lm22_add)
}
meta_plotter <- meta_plotter[-1,]
cur_plotter <- cur_plotter[-1,]
immuno_plotter <- immuno_plotter[-1,]
lm22_plotter <- lm22_plotter[-1,]
  
  
meta_plotter$condition <- factor(meta_plotter$condition, levels = c("Healthy", "Latent", "Active"))
meta_plotter <- meta_plotter[order(meta_plotter$condition),]
meta_keep <- as.character(meta_plotter$sample)

#cur
rownames(cur_plotter) <- cur_plotter[,ncol(cur_plotter)]
cur_plotter <- cur_plotter[meta_keep,1:ncol(cur_plotter)-1]

#good cols
color_picker_ours <- c(Macs_IL_4_Macs_IL_13 = colerer[11], Macs_IL_10=colerer[2],
                       Macs_IFN_B1a_Macs_IFN_y=colerer[3], DCs_IFN_B1a_DCs_IFN_y=colerer[4],
                       DCs_IL_10=colerer[5], DCs_IFN_B1a=colerer[12],
                       Monos_IL_4_Monos_IL_13=colerer[7], Monos_IL_10=colerer[8],
                       Monos_IFN_B1a_Monos_IFN_y=colerer[9], PMNs_IFN_B1a_PMNs_IFN_y=colerer[10],
                       PMNs_IL_4_PMNs_IL_13=colerer[1],
                       PMNs_IFN_B1a=colerer[13])

cur_plots <- bar_plotter2(cur_plotter, meta_plotter, category = "condition",
                          id_cat_num = 1, color_picker = color_picker_ours)
curs <- plot_grid(plotlist = cur_plots, nrow = 1)

#Immuno
rownames(immuno_plotter) <- immuno_plotter[,ncol(immuno_plotter)]
immuno_plotter <- immuno_plotter[meta_keep,1:ncol(immuno_plotter)-1]

immuno_cols <- randomColor(20)
color_picker_immuno <- c(`CD14_positive_monocyte`=immuno_cols[1], `CD16_positive_monocyte`=immuno_cols[2],
                         `CD4_positive_alpha_beta_T_cell`=immuno_cols[3], `CD56bright_natural_killer_cell`=immuno_cols[4],
                         `CD56dim_natural_killer_cell`=immuno_cols[5], `CD8_positive_alpha_beta_T_cell`=immuno_cols[6],
                         `MAST_cell`=immuno_cols[7], `basophil`=immuno_cols[8],
                         `eosinophil`=immuno_cols[9], `gamma_delta_T_cell`=immuno_cols[10],
                         `hematopoietic_progenitor`=immuno_cols[11], `macrophage_m0`=immuno_cols[12],
                         `macrophage_m1`=immuno_cols[13], `macrophage_m2`=immuno_cols[14],
                         `memory_B_cell`=immuno_cols[15], `myeloid_dendritic_cell`=immuno_cols[16],
                         `naive_B_cell`=immuno_cols[17], `neutrophil`=immuno_cols[18],
                         `plasma_cell`=immuno_cols[19], `plasmacytoid_dendritic_cell`=immuno_cols[20])

cur_plots <- bar_plotter2(immuno_plotter, meta_plotter, category = "condition",
                          id_cat_num = 1, color_picker = color_picker_immuno)
immunos <- plot_grid(plotlist = cur_plots, nrow = 1)

#lm22
rownames(lm22_plotter) <- lm22_plotter[,ncol(lm22_plotter)]
lm22_plotter <- lm22_plotter[meta_keep,1:ncol(lm22_plotter)-1]

ciber_cols <- randomColor(22)
color_picker_theirs <- c(`B.cells.naive`=ciber_cols[1], `B.cells.memory`=ciber_cols[2],
                         `Plasma.cells`=ciber_cols[3], `T.cells.CD8`=ciber_cols[4],
                         `T.cells.CD4.naive`=ciber_cols[5], `T.cells.CD4.memory.resting`=ciber_cols[6],
                         `T.cells.CD4.memory.activated`=ciber_cols[7], `T.cells.follicular.helper`=ciber_cols[8],
                         `T.cells.regulatory`=ciber_cols[9], `T.cells.gamma.delta`=ciber_cols[10],
                         `NK.cells.resting`=ciber_cols[11], `NK.cells.activated`=ciber_cols[12],
                         `Monocytes`=ciber_cols[13], `Macrophages.M0`=ciber_cols[14],
                         `Macrophages.M1`=ciber_cols[15], `Macrophages.M2`=ciber_cols[16],
                         `Dendritic.cells.resting`=ciber_cols[17], `Dendritic.cells.activated`=ciber_cols[18],
                         `Mast.cells.resting`=ciber_cols[19], `Mast.cells.activated`=ciber_cols[20],
                         `Eosinophils`=ciber_cols[21], `Neutrophils`=ciber_cols[22])

cur_plots <- bar_plotter2(lm22_plotter, meta_plotter, category = "condition",
                          id_cat_num = 1, color_picker = color_picker_theirs)
lm22 <- plot_grid(plotlist = cur_plots, nrow = 1)

##
##plot it
pdf("all_plots_3_types.pdf", height = 20, width = 20)
pp = plot_grid(plotlist=list(curs, immunos, lm22), nrow = 3)
print(pp)
dev.off()

#
#
#
###
meta_plotter_rna <- subset(meta_plotter, data_type == "RNASeq")
cur_plotter_rna <- cur_plotter[as.character(meta_plotter_rna$sample),]

cur_boxers <- box_plotter(cur_plotter_rna, meta_plotter_rna, "condition", 1, 
                            color_picker_ours, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("rna_plots/cur_boxes.pdf",
    height = 7*ceiling(length(cur_boxers)/3), 
    width = 25)
plot_grid(plotlist = cur_boxers, ncol = 3)
dev.off()

immuno_plotter_rna <- immuno_plotter[as.character(meta_plotter_rna$sample),]
immuno_boxers <- box_plotter(immuno_plotter_rna, meta_plotter_rna, "condition", 1, 
                          color_picker = color_picker_immuno, 
                          cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("rna_plots/immuno_boxes.pdf",
    height = 7*ceiling(length(immuno_boxers)/3), 
    width = 25)
plot_grid(plotlist = immuno_boxers, ncol = 3)
dev.off()


lm22_plotter_rna <- lm22_plotter[as.character(meta_plotter_rna$sample),]
lm22_boxers <- box_plotter(lm22_plotter_rna, meta_plotter_rna, "condition", 1, 
                             color_picker_theirs, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("rna_plots/lm22_boxes.pdf",
    height = 7*ceiling(length(lm22_boxers)/3), 
    width = 25)
plot_grid(plotlist = lm22_boxers, ncol = 3)
dev.off()


RNA_boxes <- plot_grid(plotlist = list(
  cur_boxers[[10]]+theme(legend.position = 'none',
                         axis.line.x = element_line(size = 1, color = "black"),
                         axis.line.y = element_line(size = 1, color = "black")), 
  cur_boxers[[11]]+theme(legend.position = 'none',
                         axis.line.x = element_line(size = 1, color = "black"),
                         axis.line.y = element_line(size = 1, color = "black")),
  immuno_boxers[[18]]+theme(legend.position = 'none',
                            axis.line.x = element_line(size = 1, color = "black"),
                            axis.line.y = element_line(size = 1, color = "black")),
  immuno_boxers[[4]]+theme(legend.position = 'none',
                           axis.line.x = element_line(size = 1, color = "black"),
                           axis.line.y = element_line(size = 1, color = "black")),
  lm22_boxers[[22]]+theme(legend.position = 'none',
                          axis.line.x = element_line(size = 1, color = "black"),
                           axis.line.y = element_line(size = 1, color = "black")),
  lm22_boxers[[11]]+theme(legend.position = 'none',
                          axis.line.x = element_line(size = 1, color = "black"),
                          axis.line.y = element_line(size = 1, color = "black"))), 
  nrow = 3, ncol = 2)




### microarray

###
meta_plotter_ma <- subset(meta_plotter, data_type == "Microarray")
cur_plotter_ma <- cur_plotter[as.character(meta_plotter_ma$sample),]

cur_boxers <- box_plotter(cur_plotter_ma, meta_plotter_ma, "condition", 1, 
                          color_picker_ours, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("ma_plots/cur_boxes.pdf",
    height = 7*ceiling(length(cur_boxers)/3), 
    width = 25)
plot_grid(plotlist = cur_boxers, ncol = 3)
dev.off()

immuno_plotter_ma <- immuno_plotter[as.character(meta_plotter_ma$sample),]
immuno_boxers <- box_plotter(immuno_plotter_ma, meta_plotter_ma, "condition", 1, 
                             color_picker = color_picker_immuno, 
                             cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("ma_plots/immuno_boxes.pdf",
    height = 7*ceiling(length(immuno_boxers)/3), 
    width = 25)
plot_grid(plotlist = immuno_boxers, ncol = 3)
dev.off()


lm22_plotter_ma <- lm22_plotter[as.character(meta_plotter_ma$sample),]
lm22_boxers <- box_plotter(lm22_plotter_ma, meta_plotter_ma, "condition", 1, 
                           color_picker_theirs, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("ma_plots/lm22_boxes.pdf",
    height = 7*ceiling(length(lm22_boxers)/3), 
    width = 25)
plot_grid(plotlist = lm22_boxers, ncol = 3)
dev.off()


MA_boxes <- plot_grid(plotlist = list(
  cur_boxers[[10]]+theme(legend.position = 'none',
                         axis.line.x = element_line(size = 1, color = "black"),
                         axis.line.y = element_line(size = 1, color = "black")), 
  cur_boxers[[11]]+theme(legend.position = 'none',
                         axis.line.x = element_line(size = 1, color = "black"),
                         axis.line.y = element_line(size = 1, color = "black")),
  immuno_boxers[[18]]+theme(legend.position = 'none',
                            axis.line.x = element_line(size = 1, color = "black"),
                            axis.line.y = element_line(size = 1, color = "black")),
  immuno_boxers[[4]]+theme(legend.position = 'none',
                           axis.line.x = element_line(size = 1, color = "black"),
                           axis.line.y = element_line(size = 1, color = "black")),
  lm22_boxers[[22]]+theme(legend.position = 'none',
                          axis.line.x = element_line(size = 1, color = "black"),
                          axis.line.y = element_line(size = 1, color = "black")),
  lm22_boxers[[11]]+theme(legend.position = 'none',
                          axis.line.x = element_line(size = 1, color = "black"),
                          axis.line.y = element_line(size = 1, color = "black"))), 
  nrow = 3, ncol = 2)
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
### make heat map of neutrophil associated genes from our curated gene set
## microarray first
meta_keep <- as.character(meta_plotter_ma$sample)

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))

n_mat <- matrix(0,45,1)
for (i in c(3:10)){
  looker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  curr <- read.table(looker, header=T, sep='\t')
  rownames(curr) <- curr[,1]
  curr <- curr[as.character(n_genes$GeneSymbol),-1]
  curr <- t(scale(t(curr)))
  n_mat <- cbind(n_mat, curr)
}
n_mat <- n_mat[,meta_keep]
rownames(n_mat) <- as.character(n_genes$GeneSymbol)



library(ComplexHeatmap)

#n_scale <- t(scale(t(n_mat)))
n_scale[is.na(n_scale)] <- 0


healthies <- which(meta_plotter_ma$condition == "Healthy")
meta_healthies <- meta_plotter_ma[healthies,]
latents <- which(meta_plotter_ma$condition == "Latent")
meta_latents <- meta_plotter_ma[latents,]
actives <- which(meta_plotter_ma$condition == "Active")
meta_actives <- meta_plotter_ma[actives,]




ha <- columnAnnotation(df = data.frame(condition = as.character(meta_healthies$condition),
                                       region = as.character(meta_healthies$region)),
                       col = list(condition = c("Healthy" = "dodgerblue2")))
healthy_heat <- Heatmap(n_scale[,healthies],
        top_annotation = ha,
        heatmap_legend_param = list(title = "Scaled Expression",
                                    title_gp = gpar(fontsize = 15),
                                    legend_height = unit(3, "cm"),
                                    labels_gp = gpar(fontsize = 14)),
        split = n_genes$Source,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = T,
        column_title = "",
        column_title_gp = gpar(fontsize = 15),
        show_column_names = F, show_row_names = T,
        row_names_gp = gpar(fontsize = 8))

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_latents$condition),
                                       region = as.character(meta_latents$region)),
                       col = list(condition = c("Latent" = "red3")))
latent_heat <- Heatmap(n_scale[,latents],
                   top_annotation = ha,
                   heatmap_legend_param = list(title = "Scaled Expression",
                                               title_gp = gpar(fontsize = 15),
                                               legend_height = unit(3, "cm"),
                                               labels_gp = gpar(fontsize = 14)),
                   split = n_genes$Source,
                   col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                   cluster_rows = T, cluster_columns = T,
                   column_title = "",
                   column_title_gp = gpar(fontsize = 15),
                   show_column_names = F, show_row_names = T,
                   row_names_gp = gpar(fontsize = 8))

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_actives$condition),
                                       region = as.character(meta_actives$region)),
                       col = list(condition = c("Active" = "purple")))
active_heat <- Heatmap(n_scale[,actives],
                  top_annotation = ha,
                  heatmap_legend_param = list(title = "Scaled Expression",
                                              title_gp = gpar(fontsize = 15),
                                              legend_height = unit(3, "cm"),
                                              labels_gp = gpar(fontsize = 14)),
                  split = n_genes$Source,
                  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                  cluster_rows = T, cluster_columns = T,
                  column_title = "",
                  column_title_gp = gpar(fontsize = 15),
                  show_column_names = F, show_row_names = T,
                  row_names_gp = gpar(fontsize = 8))

pdf("ma_plots/neutrophil_gene_matrix.pdf", height = 10, width = 12)
healthy_heat + latent_heat + active_heat
dev.off()

pdf("ma_plots/neutrophil_gene_matrix_healthy.pdf", height = 10, width = 7)
healthy_heat
dev.off()

pdf("ma_plots/neutrophil_gene_matrix_latent.pdf", height = 10, width = 10)
latent_heat
dev.off()

pdf("ma_plots/neutrophil_gene_matrix_active.pdf", height = 10, width = 10)
active_heat
dev.off()

### 
#### heatmap in order of PMNS_ifn_beta_gamma

### make heat map of neutrophil associated genes from our curated gene set

#
#
#

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))
n_genes <- subset(n_genes, GeneSymbol != "INAFM2" & GeneSymbol != "FLJ23867")

n_mat <- matrix(0,43,1)
for (i in c(3:10)){
  looker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  curr <- read.table(looker, header=T, sep='\t')
  rownames(curr) <- curr[,1]
  curr <- curr[as.character(n_genes$GeneSymbol),-1]
  curr <- t(scale(t(curr)))
  n_mat <- cbind(n_mat, curr)
}
n_mat <- n_mat[,meta_keep]
rownames(n_mat) <- as.character(n_genes$GeneSymbol)
#
#

rownames(meta_plotter_ma) <- meta_plotter_ma$sample
meta_plotter_ma <- meta_plotter_ma[meta_keep,]
meta_plotter_ma$IFN_level <- cur_plotter_ma$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter_ma$condition2 <- factor(meta_plotter_ma$condition, levels = rev(levels(meta_plotter_ma$condition)))
meta_plotter_ma <- meta_plotter_ma[order(meta_plotter_ma$condition2),]
meta_plotter_ma <- meta_plotter_ma[order(meta_plotter_ma$IFN_level, decreasing = T),]

n_mat <- n_mat[,as.character(meta_plotter_ma$sample)]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter_ma$condition),
                                       #region = as.character(meta_plotter_ma$region),
                                       Neutrophil_IFN = meta_plotter_ma$IFN_level),
                       col = list(condition = c("Healthy" = "dodgerblue2",
                                                "Latent" = "red3",
                                                "Active" = "purple"),
                                  Neutrophil_IFN = colorRamp2(c(0, 0.0001, 0.5), c("grey","white","mediumvioletred"))))
#n_mat <- t(scale(t(n_mat)))
n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                       top_annotation = ha,
                       heatmap_legend_param = list(title = "Scaled Expression",
                                                   #title_gp = gpar(fontsize = 15),
                                                   legend_height = unit(3, "cm"),
                                                   labels_gp = gpar(fontsize = 14)),
                       split = n_genes$Source,
                       col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                       cluster_rows = T, cluster_columns = F,
                       column_title = "",
                       column_title_gp = gpar(fontsize = 15),
                       show_column_names = F, show_row_names = T,
                       row_names_gp = gpar(fontsize = 8))

pdf("ma_plots/neut_heat.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
##### heatmap with each group clustered


rownames(meta_plotter_ma) <- meta_plotter_ma$sample
meta_plotter_ma <- meta_plotter_ma[meta_keep,]
cur_plotter_ma <- cur_plotter_ma[meta_keep,]
meta_plotter_ma$IFN_level <- cur_plotter_ma$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter_ma$condition2 <- factor(meta_plotter_ma$condition, levels = rev(levels(meta_plotter_ma$condition)))
meta_plotter_ma <- meta_plotter_ma[order(meta_plotter_ma$condition2),]
meta_plotter_ma <- meta_plotter_ma[order(meta_plotter_ma$IFN_level, decreasing = T),]

# reorder each group
gs <- levels(meta_plotter_ma$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter_ma, condition2 == gs[i])
  h_gs <- hclust(dist(t(n_mat[,as.character(g1$sample)])))
  adder <- rev(h_gs$labels[h_gs$order])
  new_order <- c(new_order, adder)
}

meta_plotter_ma <- meta_plotter_ma[new_order,]
n_mat <- n_mat[,new_order]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter_ma$condition),
                                       study = as.character(meta_plotter_ma$study),
                                       Neutrophil_IFN = meta_plotter_ma$IFN_level),
                       col = list(condition = c("Healthy" = "dodgerblue2",
                                                "Latent" = "red3",
                                                "Active" = "purple"),
                                  Neutrophil_IFN = colorRamp2(c(0, 0.0001, 0.5), c("grey","white","mediumvioletred"))))
#n_mat <- t(scale(t(n_mat)))
n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                     top_annotation = ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("ma_plots/neut_heat_clust.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
#
#### heatmap with barplot instead of other stuff

bars <- cur_plotter_ma[new_order,11]
bars <- cbind(bars, 0.6-bars)
colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter_ma$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

n_mat[is.na(n_mat)] <- 0
neut_heat_bars <- Heatmap(n_mat,
                     top_annotation = column_ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     column_title_gp = gpar(fontsize = 15),
                     top_annotation_height = unit(30, "mm"),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 18))

pdf("ma_plots/neut_heat_clust_bars.pdf", height = 5, width = 10)
neut_heat_bars
dev.off()



neut_heat <- grid.grabExpr(draw(neut_heat_bars))

pdf("ma_plots/MA_figure.pdf", height = 20, width = 22)
MA_fig <- grid.arrange(MA_boxes, neut_heat, nrow = 1, widths = c(1,1.5))
print(MA_fig)
dev.off()


#### clustered bars in decreaing order

# reorder each group
gs <- levels(meta_plotter_ma$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter_ma, condition2 == gs[i])
  g1 <- g1[order(g1$IFN_level, decreasing = T),]
  adder <- as.character(g1$sample)
  new_order <- c(new_order, adder)
}

meta_plotter_ma <- meta_plotter_ma[new_order,]
n_mat <- n_mat[,new_order]



bars <- cur_plotter_ma[new_order,1:12]
#bars <- cbind(bars, 0.6-bars)
#colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter_ma$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(Macs_IFN_B1a_Macs_IFN_y="purple", 
                                                                                                 Macs_IL_10=colerer[2], Macs_IL_4_Macs_IL_13 = colerer[11],
                                                                                                 DCs_IFN_B1a=colerer[12], DCs_IFN_B1a_DCs_IFN_y=colerer[4],
                                                                                                 DCs_IL_10=colerer[5],
                                                                                                 Monos_IFN_B1a_Monos_IFN_y="hotpink",
                                                                                                 Monos_IL_10=colerer[8], Monos_IL_4_Monos_IL_13=colerer[7],
                                                                                                 PMNs_IFN_B1a=colerer[13], PMNs_IFN_B1a_PMNs_IFN_y=colerer[10],
                                                                                                 PMNs_IL_4_PMNs_IL_13=colerer[1]))),
                              df = data.frame(condition = as.character(meta_plotter_ma$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                     top_annotation = column_ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     top_annotation_height = unit(40, "mm"),
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("ma_plots/neut_heat_clust_bars_decreasing.pdf", height = 5, width = 10)
neut_heat
dev.off()



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
#
#
#
#
#
######
### make heat map of neutrophil associated genes from our curated gene set
## Now RNA samples!!!


meta_keep <- as.character(meta_plotter_rna$sample)

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))

n_mat <- matrix(0,45,1)
for (i in c(1,2,11,12,13)){
  looker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  curr <- read.table(looker, header=T, sep='\t')
  rownames(curr) <- curr[,1]
  curr <- curr[as.character(n_genes$GeneSymbol),-1]
  curr <- t(scale(t(curr)))
  n_mat <- cbind(n_mat, curr)
}
n_mat <- n_mat[,meta_keep]
rownames(n_mat) <- as.character(n_genes$GeneSymbol)



library(ComplexHeatmap)

#n_scale <- t(scale(t(n_mat)))
n_scale[is.na(n_scale)] <- 0


healthies <- which(meta_plotter_rna$condition == "Healthy")
meta_healthies <- meta_plotter_rna[healthies,]
latents <- which(meta_plotter_rna$condition == "Latent")
meta_latents <- meta_plotter_rna[latents,]
actives <- which(meta_plotter_rna$condition == "Active")
meta_actives <- meta_plotter_rna[actives,]




ha <- columnAnnotation(df = data.frame(condition = as.character(meta_healthies$condition),
                                       region = as.character(meta_healthies$region)),
                       col = list(condition = c("Healthy" = "dodgerblue2")))
healthy_heat <- Heatmap(n_scale[,healthies],
                        top_annotation = ha,
                        heatmap_legend_param = list(title = "Scaled Expression",
                                                    title_gp = gpar(fontsize = 15),
                                                    legend_height = unit(3, "cm"),
                                                    labels_gp = gpar(fontsize = 14)),
                        split = n_genes$Source,
                        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                        cluster_rows = T, cluster_columns = T,
                        column_title = "",
                        column_title_gp = gpar(fontsize = 15),
                        show_column_names = F, show_row_names = T,
                        row_names_gp = gpar(fontsize = 8))

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_latents$condition),
                                       region = as.character(meta_latents$region)),
                       col = list(condition = c("Latent" = "red3")))
latent_heat <- Heatmap(n_scale[,latents],
                       top_annotation = ha,
                       heatmap_legend_param = list(title = "Scaled Expression",
                                                   title_gp = gpar(fontsize = 15),
                                                   legend_height = unit(3, "cm"),
                                                   labels_gp = gpar(fontsize = 14)),
                       split = n_genes$Source,
                       col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                       cluster_rows = T, cluster_columns = T,
                       column_title = "",
                       column_title_gp = gpar(fontsize = 15),
                       show_column_names = F, show_row_names = T,
                       row_names_gp = gpar(fontsize = 8))

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_actives$condition),
                                       region = as.character(meta_actives$region)),
                       col = list(condition = c("Active" = "purple")))
active_heat <- Heatmap(n_scale[,actives],
                       top_annotation = ha,
                       heatmap_legend_param = list(title = "Scaled Expression",
                                                   title_gp = gpar(fontsize = 15),
                                                   legend_height = unit(3, "cm"),
                                                   labels_gp = gpar(fontsize = 14)),
                       split = n_genes$Source,
                       col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                       cluster_rows = T, cluster_columns = T,
                       column_title = "",
                       column_title_gp = gpar(fontsize = 15),
                       show_column_names = F, show_row_names = T,
                       row_names_gp = gpar(fontsize = 8))

pdf("rna_plots/neutrophil_gene_matrix.pdf", height = 10, width = 12)
healthy_heat + latent_heat + active_heat
dev.off()

pdf("rna_plots/neutrophil_gene_matrix_healthy.pdf", height = 10, width = 10)
healthy_heat
dev.off()

pdf("rna_plots/neutrophil_gene_matrix_latent.pdf", height = 10, width = 10)
latent_heat
dev.off()

pdf("rna_plots/neutrophil_gene_matrix_active.pdf", height = 10, width = 10)
active_heat
dev.off()

### 
#### heatmap in order of PMNS_ifn_beta_gamma

### make heat map of neutrophil associated genes from our curated gene set

#
#
#

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))
n_genes <- subset(n_genes, GeneSymbol != "INAFM2" & GeneSymbol != "FLJ23867")

n_mat <- matrix(0,43,1)
for (i in c(1,2,11,12,13)){
  looker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  curr <- read.table(looker, header=T, sep='\t')
  rownames(curr) <- curr[,1]
  curr <- curr[as.character(n_genes$GeneSymbol),-1]
  curr <- t(scale(t(curr)))
  n_mat <- cbind(n_mat, curr)
}
n_mat <- n_mat[,meta_keep]
rownames(n_mat) <- as.character(n_genes$GeneSymbol)
#
#

rownames(meta_plotter_rna) <- meta_plotter_rna$sample
meta_plotter_rna <- meta_plotter_rna[meta_keep,]
cur_plotter_rna <- cur_plotter_rna[meta_keep,]
meta_plotter_rna$IFN_level <- cur_plotter_rna$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter_rna$condition2 <- factor(meta_plotter_rna$condition, levels = rev(levels(meta_plotter_rna$condition)))
meta_plotter_rna <- meta_plotter_rna[order(meta_plotter_rna$condition2),]
meta_plotter_rna <- meta_plotter_rna[order(meta_plotter_rna$IFN_level, decreasing = T),]

n_mat <- n_mat[,as.character(meta_plotter_rna$sample)]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter_rna$condition),
                                       #study = as.character(meta_plotter_rna$study),
                                       Neutrophil_IFN = meta_plotter_rna$IFN_level),
                       col = list(condition = c("Healthy" = "dodgerblue2",
                                                "Latent" = "red3",
                                                "Active" = "purple"),
                                  Neutrophil_IFN = colorRamp2(c(0, 0.0001, 0.5), c("grey","white","mediumvioletred"))))
#n_mat <- t(scale(t(n_mat)))
n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                     top_annotation = ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("rna_plots/neut_heat.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
##### heatmap with each group clustered


rownames(meta_plotter_rna) <- meta_plotter_rna$sample
meta_plotter_rna <- meta_plotter_rna[meta_keep,]
cur_plotter_rna <- cur_plotter_rna[meta_keep,]
meta_plotter_rna$IFN_level <- cur_plotter_rna$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter_rna$condition2 <- factor(meta_plotter_rna$condition, levels = rev(levels(meta_plotter_rna$condition)))
meta_plotter_rna <- meta_plotter_rna[order(meta_plotter_rna$condition2),]
meta_plotter_rna <- meta_plotter_rna[order(meta_plotter_rna$IFN_level, decreasing = T),]

# reorder each group
gs <- levels(meta_plotter_rna$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter_rna, condition2 == gs[i])
  h_gs <- hclust(dist(t(n_mat[,as.character(g1$sample)])))
  adder <- rev(h_gs$labels[h_gs$order])
  new_order <- c(new_order, adder)
}

meta_plotter_rna <- meta_plotter_rna[new_order,]
n_mat <- n_mat[,new_order]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter_rna$condition),
                                       study = as.character(meta_plotter_rna$study),
                                       Neutrophil_IFN = meta_plotter_rna$IFN_level),
                       col = list(condition = c("Healthy" = "dodgerblue2",
                                                "Latent" = "red3",
                                                "Active" = "purple"),
                                  Neutrophil_IFN = colorRamp2(c(0, 0.0001, 0.5), c("grey","white","mediumvioletred"))))
#n_mat <- t(scale(t(n_mat)))
n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                     top_annotation = ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("rna_plots/neut_heat_clust.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
#
#### heatmap with barplot instead of other stuff

bars <- cur_plotter_rna[new_order,11]
bars <- cbind(bars, 1-bars)
colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter_rna$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

n_mat[is.na(n_mat)] <- 0
neut_heat_bars <- Heatmap(n_mat,
                     top_annotation = column_ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     top_annotation_height = unit(30, "mm"),
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 18))

pdf("rna_plots/neut_heat_clust_bars.pdf", height = 5, width = 10)
neut_heat_bars
dev.off()


neut_heat <- grid.grabExpr(draw(neut_heat_bars))

pdf("rna_plots/RNA_figure.pdf", height = 20, width = 22)
RNA_fig <- grid.arrange(RNA_boxes, neut_heat, nrow = 1, widths = c(1,1.5))
print(RNA_fig)
dev.off()




#### clustered bars in decreaing order

# reorder each group
gs <- levels(meta_plotter_rna$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter_rna, condition2 == gs[i])
  g1 <- g1[order(g1$IFN_level, decreasing = T),]
  adder <- as.character(g1$sample)
  new_order <- c(new_order, adder)
}

meta_plotter_rna <- meta_plotter_rna[new_order,]
n_mat <- n_mat[,new_order]



bars <- cur_plotter_rna[new_order,11]
bars <- cbind(bars, 1-bars)
colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter_rna$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

bars <- cur_plotter_rna[new_order,1:12]
#bars <- cbind(bars, 1-bars)
#colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(Macs_IFN_B1a_Macs_IFN_y="purple", 
                                                                                                   Macs_IL_10=colerer[2], Macs_IL_4_Macs_IL_13 = colerer[11],
                                                                                                   DCs_IFN_B1a=colerer[12], DCs_IFN_B1a_DCs_IFN_y=colerer[4],
                                                                                                   DCs_IL_10=colerer[5],
                                                                                                   Monos_IFN_B1a_Monos_IFN_y="hotpink",
                                                                                                   Monos_IL_10=colerer[8], Monos_IL_4_Monos_IL_13=colerer[7],
                                                                                                   PMNs_IFN_B1a=colerer[13], PMNs_IFN_B1a_PMNs_IFN_y=colerer[10],
                                                                                                   PMNs_IL_4_PMNs_IL_13=colerer[1]))),
                              df = data.frame(condition = as.character(meta_plotter_rna$condition)),
                              col = list(condition = c("Healthy" = "dodgerblue2",
                                                       "Latent" = "red3",
                                                       "Active" = "purple")))

n_mat[is.na(n_mat)] <- 0
neut_heat <- Heatmap(n_mat,
                     top_annotation = column_ha,
                     heatmap_legend_param = list(title = "Scaled Expression",
                                                 #title_gp = gpar(fontsize = 15),
                                                 legend_height = unit(3, "cm"),
                                                 labels_gp = gpar(fontsize = 14)),
                     split = n_genes$Source,
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     cluster_rows = T, cluster_columns = F,
                     column_title = "",
                     top_annotation_height = unit(40, "mm"),
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("rna_plots/neut_heat_clust_bars_decreasing.pdf", height = 5, width = 10)
neut_heat
dev.off()




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
#### final TB Fig

pdf("ALL_TB_FIG_1.30.pdf", height = 20, width = 50)
grid.arrange(MA_fig, RNA_fig, nrow = 1)
dev.off()
