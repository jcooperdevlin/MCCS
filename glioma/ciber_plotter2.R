#### cibersort figure generator

# load required libraries
library(DESeq2)
library(gplots)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(randomcoloR)
library(ComplexHeatmap)
library(viridis)
library(kohonen)
library(Rtsne)
library(psych)
library(eulerr)
library(data.table)
library(randomcoloR)
library(ggpubr)
library(cowplot)
library(circlize)
library(dendextend)
source("/Users/devlij03/Google Drive/Png/geom_signif_CD.R")


good_cols <- c("#f70c1c", "#6d89e8", "#91dd68", "#482e91", "#fc9207",
               "#fcdb23", "#e87ac3", "#5beabd", "#01871e", "#a0080f", "grey80", "black", "grey50")
colerer <- c(good_cols)

ciber_tb_plotter <- function(results_df, 
                          sig_mat, bulk_mat,
                          cts, total_cts,
                          color_cat, color_groups, color_ids,
                          metadata) {
  
  ## set colors
  good_cols <- c("#f70c1c", "#6d89e8", "#91dd68", "#482e91", "#fc9207",
                 "#fcdb23", "#e87ac3", "#5beabd", "#01871e", "#a0080f", "grey80", "black", "grey50")
  colerer <- c(good_cols)
  # cus_cell_cols <- c(Macs_IL_4_Macs_IL_13 = colerer[11], Macs_IL_10=colerer[2],
  #                    Macs_IFN_B1a=colerer[3], DCs_IFN_B1a_DCs_IFN_y=colerer[4],
  #                    DCs_IL_10=colerer[5], DCs_IFN_B1a=colerer[12],
  #                    Monos_IL_4_Monos_IL_13=colerer[7], Monos_IL_10=colerer[8],
  #                    Monos_IFN_B1a=colerer[9], PMNs_IFN_B1a_PMNs_IFN_y=colerer[10],
  #                    PMNs_IL_4_PMNs_IL_13=colerer[1], PMNs_IFN_lambda2=colerer[6],
  #                    PMNs_IFN_B1a=colerer[13])
  
  cus_cell_cols <- c(Macs_IFN_B1a_Macs_IFN_y="purple", Macs_IFN_B1a=colerer[3], 
              Macs_IL_10=colerer[2], Macs_IL_4_Macs_IL_13 = colerer[11],
              DCs_IFN_B1a_DCs_IFN_y=colerer[4], DCs_IFN_B1a=colerer[12],
              DCs_IL_10=colerer[5],
              Monos_IFN_B1a_Monos_IFN_y="hotpink", Monos_IFN_B1a=colerer[9],
              Monos_IL_10=colerer[8], Monos_IL_4_Monos_IL_13=colerer[7],
              PMNs_IFN_B1a_PMNs_IFN_y=colerer[10], PMNs_IFN_B1a=colerer[13],
              PMNs_IFN_lambda2=colerer[6], PMNs_IL_4_PMNs_IL_13=colerer[1])


#### bulk heatmap!

X <- read.table(sig_mat,header=T,sep="\t",row.names=1,check.names=F)
X <- X[grep(cts, X$group),]
basis <- rownames(X)
Y <- read.table(bulk_mat, header=T, sep="\t",check.names=F)

Y_keep <- Y[which(Y$GeneSymbol %in% basis),]

print(dim(Y_keep))


kk <- data.matrix(Y_keep[,-1])
rownames(kk) = Y_keep$GeneSymbol

print(rownames(kk))
## heat for each group
groups = color_groups

plist = list()
for (i in 1:length(groups)){
  g_meta <- metadata[grep(groups[i], metadata[,color_cat]),]
  g_ids <- rownames(g_meta)

  g_bulk <- kk[,g_ids]
  g_results <- results_df[g_ids,1:total_cts]
  
  ### reorder by decreasing means?
  ct_mean <- colMeans(g_results)
  ct_mean <- ct_mean[order(ct_mean, decreasing = T)]
  
  g_results <- g_results[,names(ct_mean)]
  
  for (j in ncol(g_results):1){
    g.sorted <- g_results[order(g_results[,j], decreasing = T),]
    g_results <- g.sorted
  }
  
  #heatmap.2(g_results, trace='none', density.info = 'none',
  #          Colv = F, Rowv = F)
  
  # reorder errything
  
  new_order <- rownames(g_results)
  
  g_meta_new <- g_meta[new_order,]
  g_bulk_new <- g_bulk[,new_order]
  
  filler <- c(Macs_IFN_B1a_Macs_IFN_y="purple", 
                               Macs_IL_10=colerer[2], Macs_IL_4_Macs_IL_13 = colerer[11],
                               DCs_IFN_B1a_DCs_IFN_y=colerer[4], DCs_IFN_B1a=colerer[12],
                               DCs_IL_10=colerer[5],
                               Monos_IFN_B1a_Monos_IFN_y="hotpink",
                               Monos_IL_10=colerer[8], Monos_IL_4_Monos_IL_13=colerer[7],
                               PMNs_IFN_B1a_PMNs_IFN_y=colerer[10], PMNs_IFN_B1a=colerer[13],
                               PMNs_IL_4_PMNs_IL_13=colerer[1])
  new_filler <- filler[names(ct_mean)]
  
  ## annoatation
  df3 <- data.frame(Category = g_meta_new[,color_cat])

  add <- color_ids
  names(add) <- color_groups

  ###alt top
  column_ha = HeatmapAnnotation(foo1 = anno_barplot(g_results, bar_width = 0.85,
                                                    gp = gpar(border='blank',lty="blank", 
                                                              fill = new_filler)),
                                df = df3,
                                col = list(Category=add),
                                annotation_height = c(unit(1,"cm"), unit(20,"cm")))
  
  ## side
  df <- data.frame(id = rownames(X), 
                   group = X$group)
  rownames(df) <- df$id
  heat_ids <- rownames(g_bulk_new)
  df2 <- df[heat_ids,]
  df3 <- data.frame(Condition=df2[,-1])
  
  va <- rowAnnotation(df = df3,
                      col = list(Condition=new_filler))
  
  g_bulk_scale <- t(scale(t(g_bulk_new)))
  #g_bulk_scale <- g_bulk_new
  
  if(i == length(color_groups)){
    plist[[i]] <- Heatmap(g_bulk_scale, top_annotation = column_ha, 
                          heatmap_legend_param = list(title = "Scaled Expression"),
                          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                          cluster_rows = T, cluster_columns = F,
                          show_column_names = F, show_row_names = T,
                          top_annotation_height = unit(10, "cm"),
                          row_names_gp = gpar(fontsize = 10)) + va
  } else{
    plist[[i]] <- Heatmap(g_bulk_scale, top_annotation = column_ha, 
                          heatmap_legend_param = list(title = "Scaled Expression"),
                          #col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
                          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                          cluster_rows = T, cluster_columns = F,
                          show_column_names = F, show_row_names = F,
                          top_annotation_height = unit(10, "cm"),
                          row_names_gp = gpar(fontsize = 10))
  }
  
}


#
#
#

#### boxes
# res <- results_df[,1:total_cts]
# res <- data.frame(res, metadata[,color_cat])
# plotter <- melt(res)
# colnames(plotter) <- c("Cluster", "X2", "value")
# 
# cell_types <- unique(plotter$X2)
# groups <- combn(unique(plotter$Cluster),2)
# 
# 
# kruskal_results <- data.frame(ct = NA, g1 = NA, g2 = NA, p.val=NA)
# for (i in 1:length(cell_types)){
#   curr <- subset(plotter, X2 == cell_types[i])
#   for (j in 1:ncol(groups)){
#     curr2 <- subset(curr, Cluster == groups[1,j] | Cluster == groups[2,j])
#     p = kruskal.test(value ~ Cluster, data=curr2)$p.value
#     add <- data.frame(ct = cell_types[i],
#                       g1 = groups[1,j],
#                       g2 = groups[2,j],
#                       p.val = p)
#     kruskal_results <- rbind(kruskal_results, add)
#   }
# }
# kruskal_results <- kruskal_results[-1,]
# kruskal_results <- subset(kruskal_results, p.val < 0.05)
# 
# #
# 
# ######## indiv plot looped
# #source("/Users/jcooperdevlin/Google Drive/Png/geom_signif_CD.R")
# 
# title_cols <- c("purple", colerer[3], 
#                colerer[2], colerer[11],
#                 colerer[4], colerer[12],
#                 colerer[5],
#                 "hotpink", colerer[9],
#                 colerer[8], colerer[7],
#                 colerer[10], colerer[13],
#                 colerer[6], colerer[1])
# 
# plotter_nz <- subset(plotter, value > 0)
# plist_boxes <- list()
# for (i in 1:length(cell_types)){
#   ct_z <- subset(plotter, X2 == cell_types[i])
#   ct_nz <- subset(plotter_nz, X2 == cell_types[i])
#   
#   test_2_keep <- subset(kruskal_results, ct == cell_types[i])
#   
#   pp <- ggplot(ct_nz, aes(x=Cluster, y=value, color=Cluster)) + 
#     ylab("Proportion") + xlab("Cluster") +
#     geom_boxplot(alpha=0.75, outlier.size=0, show.legend = F) +
#     geom_jitter(alpha=0.75, width = 0.2, size=2, show.legend = F) +
#     ggtitle(cell_types[i]) +
#     #annotate("rect", xmin=0, xmax=c(length(unique(ct_nz$Disease_Condition))+1), 
#     #        ymin=c(max(ct_nz$value)+0.01), ymax = c(max(ct_nz$value)+0.02),
#     #        fill = title_cols[i]) +
#     scale_color_manual(values = c(color_ids)) +
#     theme(axis.line.x = element_line(size = 3, color = title_cols[i]),
#           axis.line.y = element_line(size = 3, color = title_cols[i]),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.text.y = element_text(size = 20, color = "black"),
#           axis.title = element_text(size = 22),
#           plot.title = element_text(size = 22, color="black", hjust = 0.5),
#           legend.position = 'none')
#   if(nrow(test_2_keep)>0){
#     test_2_do <- list()
#     for (j in 1:nrow(test_2_keep)){
#       new <- list(unlist(test_2_keep[j,2:3]))
#       test_2_do[[j]] <- unlist(new)
#     }
#     pp <- pp + geom_signif(aes(x=Cluster, y=value), data = ct_z,
#                            comparisons = test_2_do,
#                            map_signif_level = TRUE, test = kruskal.test, textsize=6,
#                            step_increase = 0.1, color = "black")
#   }
#   
#   plist_boxes[[i]] <- pp
# }
# 
# plist[[length(color_ids)+1]] <- ggarrange(plotlist = plist_boxes, nrow = 1, ncol = 1)

return(plist)
}

gene_plotter <- function(curated, Y, color_picker){
  Yer <- data.frame(gene = rownames(Y))
  combo <- merge(curated, Yer, by.x = "GeneSymbol", by.y = "gene")
  combo <- combo[order(combo$group),]
  
  orig <- table(curated_basis$group)
  new <- table(combo$group)
  compare_plot <- data.frame(condition = names(orig),
                             value = (as.numeric(new)/as.numeric(orig))*100,
                             numerator = as.numeric(new), denominator = as.numeric(orig))
  compare_plot$labeler <- paste(compare_plot$numerator, compare_plot$denominator, sep="/")
  gg = ggplot(compare_plot, aes(condition, value, fill = condition, label = labeler))+
    geom_col() + scale_fill_manual(values = color_picker) +
    geom_text(nudge_y = 1) +
    theme(axis.text.x = element_blank())
  return(gg)
}


bar_plotter <- function(results_df, metadata, category, id_cat_num, color_picker){
  
  res <- results_df[,1:c(ncol(results_df)-3)]
  
  groups <- as.character(unique(metadata[,category]))
  glist <- list()
  for (j in 1:length(groups)){
    keep_ids <- as.character(subset(metadata, metadata[,category] == groups[j])[,id_cat_num])
    
    res_sub <- res[keep_ids,]
    
    ct_mean <- colMeans(res_sub)
    ct_mean <- ct_mean[order(ct_mean, decreasing = T)]
    
    res_sub <- res_sub[,names(ct_mean)]
    
    for (k in ncol(res_sub):1){
      g.sorted <- res_sub[order(res_sub[,k], decreasing = T),]
      res_sub <- g.sorted
    }
    
    res_plotter <- melt(res_sub)
    colnames(res_plotter) <- c("id", "condition", "value")
    res_plotter$condition <- factor(res_plotter$condition, levels = colnames(res_sub))
    
    glist[[j]] <- ggplot(res_plotter, aes(x=id, y = value, fill = condition)) +
      geom_bar(position='stack', stat = 'identity') + ylab("Proportion") +
      scale_fill_manual(values = color_picker) + ggtitle(groups[j]) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
    
  }
  return(glist)
  
}

bar_plotter2 <- function(results_df, metadata, category, id_cat_num, color_picker){
  
  res <- results_df[,1:c(ncol(results_df)-3)]
  
  groups <- as.character(unique(metadata[,category]))
  glist <- list()
  for (j in 1:length(groups)){
    keep_ids <- as.character(subset(metadata, metadata[,category] == groups[j])[,id_cat_num])
    
    res_sub <- res[keep_ids,]
    
    ct_mean <- colMeans(res_sub)
    ct_mean <- ct_mean[order(ct_mean, decreasing = T)]
    
    res_sub <- res_sub[,names(ct_mean)]
    
    for (k in ncol(res_sub):1){
      g.sorted <- res_sub[order(res_sub[,k], decreasing = T),]
      res_sub <- g.sorted
    }
    res_sub <- cbind(res_sub, rownames(res_sub))
    
    res_plotter <- melt(res_sub)
    colnames(res_plotter) <- c("id", "condition", "value")
    res_plotter$condition <- factor(res_plotter$condition, levels = colnames(res_sub))
    
    glist[[j]] <- ggplot(res_plotter, aes(x=id, y = value, fill = condition)) +
      geom_bar(position='stack', stat = 'identity') + ylab("Proportion") +
      scale_fill_manual(values = color_picker) + ggtitle(groups[j]) +
      guides(fill=guide_legend(ncol=1))+
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
    
  }
  return(glist)
  
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, c(lapply(plots, function(x)
      x + theme(legend.position="none")), nrow = 1)),
    legend,
    nrow = 2, heights = c(2,1))#,
  #heights = unit.c(unit(1, "npc") - lheight, lheight))
  #heights = unit(1, "npc") - lheight)
}


box_plotter <- function(results_df, metadata, category, id_cat_num, color_picker, cat_cols){
  res <- results_df[,1:c(ncol(results_df)-3)]
  res <- data.frame(res, metadata[,category])
  plotter <- melt(res)
  colnames(plotter) <- c("Group", "X2", "value")
  
  #plotter$X2 <- gsub("\\.", " ", plotter$X2)

  cell_types <- as.character(unique(plotter$X2))
  groups <- combn(unique(plotter$Group),2)


  kruskal_results <- data.frame(ct = NA, g1 = NA, g2 = NA, p.val=NA)
  for (i in 1:length(cell_types)){
    curr <- subset(plotter, X2 == cell_types[i])
    for (j in 1:ncol(groups)){
      curr2 <- subset(curr, Group == groups[1,j] | Group == groups[2,j])
      p = kruskal.test(value ~ Group, data=curr2)$p.value
      add <- data.frame(ct = cell_types[i],
                        g1 = groups[1,j],
                        g2 = groups[2,j],
                        p.val = p)
      kruskal_results <- rbind(kruskal_results, add)
    }
  }
  kruskal_results <- kruskal_results[-1,]
  kruskal_results <- subset(kruskal_results, p.val < 0.05)

  #

  ######## indiv plot looped
  #source("/Users/jcooperdevlin/Google Drive/Png/geom_signif_CD.R")

  plotter_nz <- subset(plotter, value > 0)
  plist_boxes <- list()
  for (i in 1:length(cell_types)){
    ct_z <- subset(plotter, X2 == cell_types[i])
    ct_nz <- subset(plotter_nz, X2 == cell_types[i])

    test_2_keep <- subset(kruskal_results, ct == cell_types[i])

    title_cols <- color_picker[which(names(color_picker) == cell_types[i])]
    
    pp <- ggplot(ct_z, aes(x=Group, y=value, color=Group)) +
      ylab("Proportion") + xlab("Group") +
      geom_boxplot(alpha=0, outlier.size=NA, show.legend = F) +
      geom_jitter(alpha=0.75, width = 0.2, size=2, show.legend = T) +
      ggtitle(cell_types[i]) +
      guides(color = guide_legend(override.aes = list(size=3))) +
      #annotate("rect", xmin=0, xmax=c(length(unique(ct_nz$Disease_Condition))+1),
      #        ymin=c(max(ct_nz$value)+0.01), ymax = c(max(ct_nz$value)+0.02),
      #        fill = title_cols[i]) +
      scale_color_manual(values = c(cat_cols)) +
      theme(axis.line.x = element_line(size = 3, color = title_cols),
            axis.line.y = element_line(size = 3, color = title_cols),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 20, color = "black"),
            axis.title = element_text(size = 22),
            plot.title = element_text(size = 22, color="black", hjust = 0.5),
            legend.position = 'right',
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 22))
    if(nrow(test_2_keep)>0){
      test_2_do <- list()
      for (j in 1:nrow(test_2_keep)){
        new <- list(unlist(test_2_keep[j,2:3]))
        test_2_do[[j]] <- unlist(new)
      }
      pp <- pp + geom_signif(aes(x=Group, y=value), data = ct_z,
                             comparisons = test_2_do,
                             map_signif_level = TRUE, test = kruskal.test, textsize=6,
                             step_increase = 0.1, color = "black")
    }

    plist_boxes[[i]] <- pp
  }
  return(plist_boxes)
}
