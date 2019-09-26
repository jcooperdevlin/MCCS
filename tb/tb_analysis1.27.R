###### Full TB Analysis for Manuscript 1.26

##### 

setwd("~/Google Drive/Png/cooper/Eblood_Submission/archive/tb_1.25/")

library(ggplot2)
library(reshape2)
library("illuminaHumanv4.db")
library("GEOquery")
library(randomcoloR)
source("functions/microarray2gene.R")
source("/Users/devlij03/Google Drive/Png/cooper/ciber/ciber_source/CIBERSORT.R")
source("/Users/devlij03/Google Drive/Png/cooper/Eblood_Submission/ciber_plotter2.R")



##GET full list of necessary genes
bases <- list.files("~/Google Drive/Png/cooper/cibersort1.23/basis/", full.names = T)
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

ma_sets <- list.files("raw", full.names = T)
gse_names <- gsub("raw/", "", ma_sets)
gse_names <- gsub("_series_.*", "", gse_names)

# for (i in 1:length(ma_sets)){
#   curr <- getGEO(filename = ma_sets[i])
#   expr <- curr@assayData$exprs
#   gener <- strsplit(curr@featureData@data$gene_assignment, "//")
#   gener2 <- unlist(lapply(gener, "[", 2))
#   gener2 <- gsub(" ", "", gener2)
#   expr <- aggregate(expr, by = list(gener2), "mean")
#   colnames(expr)[1] <- "GeneSymbol"
#   expr_keep <- expr[expr[,1] %in% basis_genes,]
# 
#   enamer <- paste0("genesymbols/", gse_names[i], "_genesymbols.txt")
#   write.table(expr_keep, enamer, sep='\t', quote=F, row.names=F)
# 
#   meta <- curr@phenoData@data
#   mnamer <- paste0("metadata/", gse_names[i], "_metadata.txt")
#   write.table(meta, mnamer, sep='\t', quote=F, row.names=F)
# 
#   for (j in 1:length(bases)){
#   tb_results <- CIBERSORT(bases[j], enamer,
#                           perm=100, QN=T,
#                           absolute=F, abs_method='sig_score')
#   res <- data.frame(tb_results[[1]])
#   res$sample <- rownames(res)
# 
#   namer <- paste0("ciber/", gse_names[i], "_", test[j], ".txt")
#   write.table(res, file=namer, sep='\t', row.names=F, quote = F)
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

for (i in c(1:length(gse_names))){
  looker <- list.files("ciber", pattern = gse_names[i], full.names = T)
  meta <- read.table(
    list.files("metadata", pattern = gse_names[i], full.names = T),
    header=T, sep='\t', quote="", fill=T)
  id_cat_num <- grep("geo_accession", colnames(meta))
  meta <- subset(meta, meta[,id_cat] %in% meta_specs[[i]])
  if (i==1){
    meta[,id_cat] <- gsub("Blood", "blood", meta[,id_cat])
  }
  
  meta_keep <- as.character(meta[,id_cat_num])
  
  #cur
  
  cur <- read.table(looker[1], header=T, sep='\t')
  rownames(cur) <- cur[,ncol(cur)]
  cur <- cur[meta_keep,1:ncol(cur)-1]

  #good cols
  color_picker_ours <- c(Macs_IL_4_Macs_IL_13 = colerer[11], Macs_IL_10=colerer[2],
                         Macs_IFN_B1a_Macs_IFN_y=colerer[3], DCs_IFN_B1a_DCs_IFN_y=colerer[4],
                         DCs_IL_10=colerer[5], DCs_IFN_B1a=colerer[12],
                         Monos_IL_4_Monos_IL_13=colerer[7], Monos_IL_10=colerer[8],
                         Monos_IFN_B1a_Monos_IFN_y=colerer[9], PMNs_IFN_B1a_PMNs_IFN_y=colerer[10],
                         PMNs_IL_4_PMNs_IL_13=colerer[1],
                         PMNs_IFN_B1a=colerer[13])
  
  cur_plots <- bar_plotter2(cur, meta, category = id_cat,
                           id_cat_num = id_cat_num, color_picker = color_picker_ours)
  curs <- plot_grid(plotlist = cur_plots, nrow = 1)
  
  #Immuno
  cur <- read.table(looker[2], header=T, sep='\t')
  rownames(cur) <- cur[,ncol(cur)]
  cur <- cur[meta_keep,1:ncol(cur)-1]
  
  cur_plots <- bar_plotter2(cur, meta, category = id_cat,
                            id_cat_num = id_cat_num, color_picker = c(good_cols, randomColor(ncol(cur))))
  immunos <- plot_grid(plotlist = cur_plots, nrow = 1)
  
  #lm22
  
  cur <- read.table(looker[3], header=T, sep='\t')
  rownames(cur) <- cur[,ncol(cur)]
  cur <- cur[meta_keep,1:ncol(cur)-1]
  
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
  
  cur_plots <- bar_plotter2(cur, meta, category = id_cat,
                            id_cat_num = id_cat_num, color_picker = color_picker_theirs)
  lm22 <- plot_grid(plotlist = cur_plots, nrow = 1)
  
  ##
  ##plot it
  namer <- paste0("plots/", gse_names[i], ".pdf")
  pdf(namer, height = 20, width = 20)
  pp = plot_grid(plotlist=list(curs, immunos, lm22), nrow = 3)
  print(pp)
  dev.off()
}



#
#
#
#
### try to combine studies active, latent and control only!

meta_plotter <- data.frame(sample=NA,condition=NA, region=NA, study = NA)
cur_plotter <- read.table("ciber/GSE39940_curated.txt", header=T, sep='\t')[1,]
immuno_plotter <- read.table("ciber/GSE39940_immunoStates.txt", header=T, sep='\t')[1,]
lm22_plotter <- read.table("ciber/GSE39940_LM22.txt", header=T, sep='\t')[1,]
for (i in c(1:7)){
  looker <- list.files("ciber", pattern = gse_names[i], full.names = T)
  meta <- read.table(
    list.files("metadata", pattern = gse_names[i], full.names = T),
    header=T, sep='\t', quote="", fill=T)
  id_cat_num <- grep("geo_accession", colnames(meta))
  id_cat <- "source_name_ch1"
  
  meta <- subset(meta, meta[,id_cat] %in% meta_specs[[i]])
  if (i==1){
    meta[,id_cat] <- gsub("Blood", "blood", meta[,id_cat])
    meta <- meta[which(grepl("long", meta[,1])==F),]
    
    thing <- meta[,id_cat]
    thing <- gsub(".*healthy.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Latent.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- meta$characteristics_ch1.4
  } else if (i==2) {
    thing <- meta[,id_cat]
    thing <- gsub(".*NID.*", "Healthy", thing)
    thing <- gsub(".*LTBI.*", "Latent", thing)
    thing <- gsub(".*TB.*", "Active", thing)
    meta$condition <- thing
    region_add <- rep("geographical region: West Africa", nrow(meta))
  } else if (i==7){
    thing <- meta[,id_cat]
    thing <- gsub(".*healthy.*", "Healthy", thing, ignore.case = T)
    thing <- gsub(".*Latent.*", "Latent", thing, ignore.case = T)
    thing <- gsub(".*Active.*", "Active", thing, ignore.case = T)
    meta$condition <- thing
    region_add <- rep("geographical region: Venezuela", nrow(meta))
    
  } else {
    thing <- meta[,id_cat]
    thing <- gsub(".*LTBI.*", "Latent", thing)
    thing <- gsub(".*Active.*", "Active", thing)
    meta$condition <- thing
    region_add <- meta$characteristics_ch1.2
  }
  meta_keep <- as.character(meta[,id_cat_num])
  study_add <- rep(gse_names[i], nrow(meta))
  
  adder <- data.frame(sample = meta[,id_cat_num], condition = meta$condition,
                      region = region_add, study = study_add)
 
  meta_plotter <- rbind(meta_plotter, adder)
  
  cur_plotter <- rbind(cur_plotter, read.table(looker[1], header=T, sep='\t'))
  immuno_plotter <- rbind(immuno_plotter, read.table(looker[2], header=T, sep='\t'))
  lm22_plotter <- rbind(lm22_plotter, read.table(looker[3], header=T, sep='\t'))
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
cur_boxers <- box_plotter(cur_plotter, meta_plotter, "condition", 1, 
                            color_picker_ours, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("cur_boxes.pdf",
    height = 7*ceiling(length(cur_boxers)/3), 
    width = 25)
plot_grid(plotlist = cur_boxers, ncol = 3)
dev.off()

immuno_boxers <- box_plotter(immuno_plotter, meta_plotter, "condition", 1, 
                          color_picker = color_picker_immuno, 
                          cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("immuno_boxes.pdf",
    height = 7*ceiling(length(immuno_boxers)/3), 
    width = 25)
plot_grid(plotlist = immuno_boxers, ncol = 3)
dev.off()


lm22_boxers <- box_plotter(lm22_plotter, meta_plotter, "condition", 1, 
                             color_picker_theirs, cat_cols = c("dodgerblue2", "red3", "purple"))

pdf("lm22_boxes.pdf",
    height = 7*ceiling(length(lm22_boxers)/3), 
    width = 25)
plot_grid(plotlist = lm22_boxers, ncol = 3)
dev.off()



#
#
#
#

#
### make heat map of neutrophil associated genes from our curated gene set

#
#
#

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))

n_mat <- matrix(0,45,1)
for (i in c(1:7)){
  looker <- list.files("genesymbols", pattern = gse_names[i], full.names = T)
  curr <- read.table(looker, header=T, sep='\t')
  rownames(curr) <- curr[,1]
  curr <- curr[as.character(n_genes$GeneSymbol),-1]
  n_mat <- cbind(n_mat, curr)
}
n_mat <- n_mat[,meta_keep]
rownames(n_mat) <- as.character(n_genes$GeneSymbol)



library(ComplexHeatmap)

n_scale <- t(scale(t(n_mat)))
n_scale[is.na(n_scale)] <- 0


healthies <- which(meta_plotter$condition == "Healthy")
meta_healthies <- meta_plotter[healthies,]
latents <- which(meta_plotter$condition == "Latent")
meta_latents <- meta_plotter[latents,]
actives <- which(meta_plotter$condition == "Active")
meta_actives <- meta_plotter[actives,]




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

pdf("neutrophil_gene_matrix.pdf", height = 10, width = 12)
healthy_heat + latent_heat + active_heat
dev.off()

pdf("neutrophil_gene_matrix_healthy.pdf", height = 10, width = 7)
healthy_heat
dev.off()

pdf("neutrophil_gene_matrix_latent.pdf", height = 10, width = 10)
latent_heat
dev.off()

pdf("neutrophil_gene_matrix_active.pdf", height = 10, width = 10)
active_heat
dev.off()


#### heatmap in order of PMNS_ifn_beta_gamma

### make heat map of neutrophil associated genes from our curated gene set

#
#
#

curated_basis <- read.table("archive/curated_basis12_source.txt", header=T, sep='\t')
n_genes <- subset(curated_basis, grepl("PMNs", curated_basis$Source))
n_genes <- subset(n_genes, GeneSymbol != "INAFM2" & GeneSymbol != "FLJ23867")

n_mat <- matrix(0,43,1)
for (i in c(1:7)){
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

rownames(meta_plotter) <- meta_plotter$sample
meta_plotter <- meta_plotter[meta_keep,]
meta_plotter$IFN_level <- cur_plotter$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter$condition2 <- factor(meta_plotter$condition, levels = rev(levels(meta_plotter$condition)))
meta_plotter <- meta_plotter[order(meta_plotter$condition2),]
meta_plotter <- meta_plotter[order(meta_plotter$IFN_level, decreasing = T),]

n_mat <- n_mat[,as.character(meta_plotter$sample)]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter$condition),
                                       #region = as.character(meta_plotter$region),
                                       Neutrophil_IFN = meta_plotter$IFN_level),
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

pdf("neut_heat.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
##### heatmap with each group clustered


rownames(meta_plotter) <- meta_plotter$sample
meta_plotter <- meta_plotter[meta_keep,]
meta_plotter$IFN_level <- cur_plotter$PMNs_IFN_B1a_PMNs_IFN_y
meta_plotter$condition2 <- factor(meta_plotter$condition, levels = rev(levels(meta_plotter$condition)))
meta_plotter <- meta_plotter[order(meta_plotter$condition2),]
meta_plotter <- meta_plotter[order(meta_plotter$IFN_level, decreasing = T),]

# reorder each group
gs <- levels(meta_plotter$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter, condition2 == gs[i])
  h_gs <- hclust(dist(t(n_mat[,as.character(g1$sample)])))
  adder <- rev(h_gs$labels[h_gs$order])
  new_order <- c(new_order, adder)
}

meta_plotter <- meta_plotter[new_order,]
n_mat <- n_mat[,new_order]

ha <- columnAnnotation(df = data.frame(condition = as.character(meta_plotter$condition),
                                       #region = as.character(meta_plotter$region),
                                       Neutrophil_IFN = meta_plotter$IFN_level),
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

pdf("neut_heat_clust.pdf", height = 5, width = 10)
neut_heat
dev.off()


#
#
#
#### heatmap with barplot instead of other stuff

bars <- cur_plotter[new_order,11]
bars <- cbind(bars, 0.6-bars)
colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter$condition)),
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
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("neut_heat_clust_bars.pdf", height = 5, width = 10)
neut_heat
dev.off()



#### clustered bars in decreaing order

# reorder each group
gs <- levels(meta_plotter$condition2)
new_order <- c()
for ( i in 1:length(gs)){
  g1 <- subset(meta_plotter, condition2 == gs[i])
  g1 <- g1[order(g1$IFN_level, decreasing = T),]
  adder <- as.character(g1$sample)
  new_order <- c(new_order, adder)
}

meta_plotter <- meta_plotter[new_order,]
n_mat <- n_mat[,new_order]



bars <- cur_plotter[new_order,11]
bars <- cbind(bars, 0.6-bars)
colnames(bars) <- c("IFN_level", "Other")
column_ha = HeatmapAnnotation(foo1 = anno_barplot(bars, bar_width = 0.85,
                                                  gp = gpar(border='blank',lty="blank", fill = c(IFN_level = "mediumvioletred", Other="grey"))),
                              df = data.frame(condition = as.character(meta_plotter$condition)),
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
                     column_title_gp = gpar(fontsize = 15),
                     show_column_names = F, show_row_names = T,
                     row_names_gp = gpar(fontsize = 8))

pdf("neut_heat_clust_bars_decreasing.pdf", height = 5, width = 10)
neut_heat
dev.off()





### with boxplot
# new_neut <- Heatmap(n_mat,
#                      #top_annotation = ha,
#                      show_heatmap_legend = F,
#                      split = n_genes$Source,
#                      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#                      cluster_rows = T, cluster_columns = F,
#                      column_title = "",
#                      column_title_gp = gpar(fontsize = 15),
#                      show_column_names = F, show_row_names = T,
#                      row_names_gp = gpar(fontsize = 8))
# neut_grob = grid.grabExpr(draw(new_neut))
# 
