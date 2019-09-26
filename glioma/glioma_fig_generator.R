###### 
####### Glioma Survival Manuscript figure 5

### 2.17.19

library(ggplot2)
library(gridExtra)
library(cowplot)
library(caret)
library(precrec)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(ggplot2)
library(plotROC)
library(reshape)
library(ggsignif)

normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}

better_col_vector <- c(Macs_IFN_B1a_Macs_IFN_y="#FFDDCD",Macs_IL_10="#FE9A66",
                       Macs_IL_4_Macs_IL_13 = "#FF5600", 
                       DCs_IFN_B1a="#9AD4F2",
                       DCs_IFN_B1a_DCs_IFN_y="#0593E0", DCs_IL_10="#03496F", 
                       Monos_IFN_B1a_Monos_IFN_y="#A0D9C6", Monos_IL_10="#41B48B",
                       Monos_IL_4_Monos_IL_13="#0C5037", 
                       PMNs_IFN_B1a="#EC9FC1",
                       PMNs_IFN_B1a_PMNs_IFN_y="#D94183",
                       PMNs_IL_4_PMNs_IL_13="#9B0A4B")

setwd("/Users/devlij03/Google Drive/Png/cooper/PNAS_manuscript/int/glioma")


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
load("/Users/devlij03/Google Drive/Png/cooper/PNAS_manuscript/int/glioma/survival/Glioma_results_ours_survival_curves.RData")

MonosIFN_curve <- glist[[7]]
MonosIFN_curve$plot$labels$title <- "Monos_IFN_B1a_Monos_IFN_y"
MonosIFN_curve$plot$labels$colour <- "Proportion"
MonosIFN_curve$plot$labels$fill <- "Proportion"
MonosIFN_curve$table$labels$y <- "Proportion"

PMNsIL4.13_curve <- glist[[12]]
PMNsIL4.13_curve$plot$labels$title <- "PMNs_IL_4_PMNs_IL_13"
PMNsIL4.13_curve$plot$labels$fill <- "Proportion"
PMNsIL4.13_curve$plot$labels$colour <- "Proportion"
PMNsIL4.13_curve$table$labels$y <- "Proportion"
PMNsIL4.13_curve$plot$labels

figA <- arrange_ggsurvplots(list(MonosIFN_curve, PMNsIL4.13_curve), print = F, nrow = 1, ncol = 2)

pdf("fig/fig_5A.pdf", height = 4, width = 8)
figA
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

### Figure 5B
###
### use machine learning results to select features of most importance for survival prediction
#
#
####load in necessary files and process them exactly as when models were generated

trains_list <- c("glioma_train_list2yr.txt","glioma_train_list5yr.txt","glioma_train_list10yr.txt")
tests_list <- c("glioma_test_list2yr.txt","glioma_test_list5yr.txt","glioma_test_list10yr.txt")

cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_metadata <- "metadata/Glioma_metadata.txt"
sig_matrix <-list(
  c("curated_basis12.txt"),
  c("immunoStates_basis.txt"),
  c("LM22.txt"))

model_list <- list.files("curated_results/", pattern = ".rds", full.names=T)

all_auc <- data.frame(Model = NA, Basis = NA,
                      truePos = NA, falsePos = NA, falseNeg = NA, trueNeg = NA,
                      TestAccuracy = NA, Sensitivity = NA, Specificity = NA, Precision = NA,
                      Recall = NA, BalancedAccuracy = NA, TestROC_AUC = NA, TestPR_AUC = NA)
ssfort_total <- data.frame(x=NA, y=NA, modname=NA, dsid = NA,
                           dsid_modname = NA, curvetype = NA, AUC = NA, model = NA, Basis = NA)
for (i in 1:length(model_list)){
  c_mod <- model_list[[i]]
  i_mod <- paste0("immunoStates_results/", basename(model_list[[i]]))
  l_mod <- paste0("lm22_results/", basename(model_list[[i]]))
  
  av_models <- list(c_mod, i_mod, l_mod)
  
  if(grepl("none", model_list[[i]])){
    f=1
  }else {f=2}
  
  surv_num <-  gsub("yrs.*", "", basename(model_list[[i]]))
  surv_num <- as.numeric(sub(".*_", "", surv_num))
  
  for (j in 1:3){
    
    cancer_df <- read.table(cancer_files, header=T)
    
    X <- read.table(sig_matrix[[j]],header=T,sep="\t")
    basis <- unique(as.character(X$GeneSymbol))
    
    cancer_df <- subset(cancer_df, GeneSymbol %in% basis)
    colnames(cancer_df) <- gsub("\\.", "-", colnames(cancer_df))
    colnames(cancer_df) <- make.unique(substring(colnames(cancer_df), 0, 12))
    rownames(cancer_df) <- cancer_df$GeneSymbol
    cancer_df <- cancer_df[,-1]
    
    cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
    rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))
    
    cancer_meta <- cancer_meta
    cancer_df <- cancer_df
    
    cancer_meta$orig_vital <- cancer_meta$vitalstatus
    
    new_vital <- cancer_meta$vitalstatus
    new_vital[cancer_meta$daystodeath > 365*surv_num] <- 0
    cancer_meta$newvital <- new_vital
    
    test_read <- tests_list[grep(surv_num, tests_list)]
    test_group <- as.character(read.table(test_read, T, '\t')$Sample)
    train_read <- trains_list[grep(surv_num, trains_list)]
    train_group <- as.character(read.table(train_read, T, '\t')$Sample)
    
    test_meta <- cancer_meta[test_group,]
    train_meta <- cancer_meta[train_group,]
    
    cancer_test <- cancer_df[,test_group]
    test_class <- test_meta$newvital
    cancer_train <- cancer_df[,train_group]
    train_class <- train_meta$newvital
    
    ### rename
    data.train <- cancer_train
    data.test <- cancer_test
    #
    ### normalization methods
    
    data_train_none <- data.train
    data_test_none <- data.test
    
    data_train_norm <- data.frame(t(apply(data.train, 1, normalize))) #min max between 0-1
    data_train_norm[is.na(data_train_norm)] <- 0
    data_test_norm <- data.frame(t(apply(data.test, 1, normalize))) #min max between 0-1
    data_test_norm[is.na(data_test_norm)] <- 0
    
    # make list for ML
    data_train_list <- list(data_train_none, data_train_norm)
    data_test_list <- list(data_test_none, data_test_norm)
    
    
    norm_names <- c("none", "min.max_norm")
    
    data.traink <- data_train_list[[f]]
    train_class <- make.names(train_class)
    data_run1 <- data.frame(t(data.traink), condition=train_class)
    data_test1 <- data_test_list[[f]]
    test_class <- make.names(test_class)
    data_test1 <- data.frame(t(data_test1), condition=test_class)
    
    mod_test <- readRDS(av_models[[j]])
    
    if(grepl("lasso", model_list[[i]])){
      
      preds1 <- data.frame(X1=predict(mod_test, newdata = data_test1))
      
      
    } else {
    preds1 <- predict(mod_test, newdata = data_test1, type = 'prob')
    }
    #preds1$obs <- predict(mod_test, newdata = data_test1)
    preds1$true <- data_test1$condition
    
    
    sscurves <- evalmod(scores = preds1$X1, labels = preds1$true)
    aucs <- round(auc(sscurves)$aucs,4)

    obs <- preds1$X1
    obs[obs>0.5]<- "X1"
    obs[obs<0.5]<- "X0"
    cf <- confusionMatrix(data = factor(obs, levels=c("X0","X1")), 
                          reference = as.factor(data_test1$condition), 
                          mode = "prec_recall")
    c_mat <- c(truePos=cf$table[1,1], falsePos=cf$table[2,1], 
               falseNeg=cf$table[1,2], trueNeg=cf$table[2,2])
    accuracies <- c(cf$overall[1], cf$byClass[1], cf$byClass[2], 
                    cf$byClass[5], cf$byClass[6], cf$byClass[11])
    names(accuracies)[1] <- "TestAccuracy"
    names(accuracies)[6] <- "BalancedAccuracy"
    
    ssfort <- fortify(sscurves)
    ssfort$AUC <- c(rep(aucs[1], length(which(ssfort$curvetype=="ROC")==T)), rep(aucs[2], length(which(ssfort$curvetype=="PRC")==T)))
    
    ssfort$model <- rep(basename(model_list[[i]]), nrow(ssfort))
    basis <- gsub(".txt", "",sig_matrix[[j]])
    ssfort$Basis <- rep(gsub("_basis", "", basis), nrow(ssfort))
    
    ssfort_total <- rbind(ssfort_total, ssfort)
    
    curr_model <- gsub(".rds", "",basename(model_list[[i]]))
    
    
    auc_adder <- cbind(data.frame(Model = curr_model, Basis = basis), t(c_mat), t(accuracies),
                            data.frame(TestROC_AUC = aucs[1], TestPR_AUC = aucs[2]))
    
    all_auc <- rbind(all_auc, auc_adder)
  }
}
ssfort_total <- ssfort_total[-1,]
all_auc <- all_auc[-1,]
write.table(all_auc, "AUC_stats_all.txt", sep='\t', quote=F, row.names = F)

year <- gsub("yrs.*", "", ssfort_total$model)
ssfort_total$Year <- sub(".*_", "", year)

modtype1  <- gsub("10_(.+)yrs.*", "10_", ssfort_total$model)
modtype1  <- gsub("7_(.+)yrs.*", "7_", modtype1)
modtype1  <- gsub("5_(.+)yrs.*", "5_", modtype1)
modtype1  <- gsub("3_(.+)yrs.*", "3_", modtype1)

modtype2  <- gsub(".*yrs_", "", ssfort_total$model)

ssfort_total$modType <- paste0(modtype1, modtype2)

all_list <- list()
counter = 1
uniqs <- unique(ssfort_total$modType)
for( i in 1:length(uniqs)){
  curr <- subset(ssfort_total, modType == uniqs[i])
  
  glist=list()
  for ( j in c(2,5,10)){
    curr2 <- subset(curr, Year == j)
    auc_c <- subset(curr2, Basis == "curated12" & curvetype == "ROC")$AUC[1]
    auc_i <- subset(curr2, Basis == "immunoStates" & curvetype == "ROC")$AUC[1]
    auc_l <- subset(curr2, Basis == "LM22" & curvetype == "ROC")$AUC[1]
    
    #test <- subset(curr2, Basis == "curated12" & curvetype == "ROC")
    #roc <- roc.curve(scores.class0 = test$y, scores.class1 = test$x, curve = T)
    #plot(roc)
    
    # PR Curve
    #pr <- pr.curve(scores.class0 = test$y, scores.class1 = test$x, curve = T)
    #plot(pr)
    
    gROC <- ggplot(subset(curr2, curvetype=="ROC"), aes(x=x, y=y, color = Basis)) + 
      geom_line() +
      annotate("text", 0.75, 0.55, label = auc_c) +
      annotate("text", 0.75, 0.5, label = auc_i) +
      annotate("text", 0.75, 0.45, label = auc_l) +
      ggtitle(paste0(gsub(".rds", "", uniqs[i]), "_", j, "yrs"))
    
    auc_c <- subset(curr2, Basis == "curated12" & curvetype == "PRC")$AUC[1]
    auc_i <- subset(curr2, Basis == "immunoStates" & curvetype == "PRC")$AUC[1]
    auc_l <- subset(curr2, Basis == "LM22" & curvetype == "PRC")$AUC[1]
    gPRC <- ggplot(subset(curr2, curvetype=="PRC"), aes(x=x, y=y, color = Basis)) + 
      geom_line() +
      annotate("text", 0.75, 0.35, label = auc_c) +
      annotate("text", 0.75, 0.3, label = auc_i) +
      annotate("text", 0.75, 0.25, label = auc_l) +
      ggtitle(paste0(gsub(".rds", "", uniqs[i]), "_", j, "yrs"))
    
    glist[[j]] <- grid.arrange(gROC, gPRC, nrow = 2)
    all_list[[counter]] <- grid.arrange(gROC, gPRC, nrow = 2)
    counter=+1
  }
  
  newlist <- list(glist[[2]], glist[[5]], glist[[10]])
  namer <- paste0("plots/", gsub(".rds", "", uniqs[i]), "_AUC.pdf")
  pdf(namer, height = 10, width = 15)
  pp = plot_grid(plotlist = newlist, nrow = 1)
  print(pp)
  dev.off()
}

#

### clean figure just curated 3 diff years ###3 4.18 update only 2 and 5 year models
curr3 <- subset(ssfort_total, Basis == "curated12")
#curr3$Year <- factor(curr3$Year, levels = c(2,5,10))
curr3$Year <- factor(curr3$Year, levels = c(2,5))
auc_2 <- subset(curr3, Year == 2 & curvetype == "ROC")$AUC[1]
auc_5 <- subset(curr3, Year == 5 & curvetype == "ROC")$AUC[1]
#auc_10 <- subset(curr3, Year == 10 & curvetype == "ROC")$AUC[1]

gROC <- ggplot(subset(curr3, curvetype=="ROC"), aes(x=x, y=y, linetype = Year)) + 
  geom_line() + scale_linetype_manual(values = c("solid", "longdash", "dotted"))+
  annotate("text", 0.75, 0.55, label = auc_2) +
  annotate("text", 0.75, 0.5, label = auc_5) +
  annotate("text", 0.75, 0.45, label = auc_10) +
  ggtitle("Lasso Model\n7-fold Cross Validation\nMin-Max Normalization")

auc_2 <- subset(curr3, Year == 2 & curvetype == "PRC")$AUC[1]
auc_5 <- subset(curr3, Year == 5 & curvetype == "PRC")$AUC[1]
auc_10 <- subset(curr3, Year == 10 & curvetype == "PRC")$AUC[1]
gPRC <- ggplot(subset(curr3, curvetype=="PRC"), aes(x=x, y=y, linetype = Year)) + 
  geom_line() + scale_linetype_manual(values = c("solid", "longdash", "dotted"))+
  annotate("text", 0.75, 0.5, label = auc_2) +
  annotate("text", 0.75, 0.45, label = auc_5) +
  annotate("text", 0.75, 0.4, label = auc_10) +
  ggtitle("Lasso Model\n7-fold Cross Validation\nMin-Max Normalization")


pdf("clean_AUC_curated.pdf", height = 5, width = 10)
grid.arrange(gROC, gPRC, nrow = 1)
dev.off()

pdf("fig/fig_5B.pdf", height = 5, width = 5)
gROC
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
### Figure 5C
###
### use machine learning results to select features of most importance for survival prediction

### load best model

bestMod <- readRDS("curated_results/lasso_model_7_5yrs_min.max_norm.rds")
bestModParam <- unlist(bestMod$bestTune)
selectedIndices <- bestMod$pred$mtry == bestModParam
huh <- bestMod$pred[selectedIndices, ]


# Best features bars
bestImp <- varImp(bestMod)$importance
bestImp$Genes <- rownames(bestImp)
bestImp <- bestImp[order(bestImp$Overall, decreasing = T),]
bestImpHigh <- bestImp[1:20,]
#bestImpHigh$Genes <- factor(bestImpHigh$Genes, levels = rev(bestImpHigh$Genes))

#merge with cell type + stim colors
basisSource <- read.table("curated_basis12_source.txt",
                          header=T, sep='\t')
rownames(basisSource) <- basisSource$GeneSymbol
basisSource2 <- subset(basisSource, GeneSymbol %in% bestImpHigh$Genes)
#basisSource2 <- basisSource2[as.character(bestImpHigh$Genes),]
bestImpHigh <- bestImpHigh[as.character(basisSource2$GeneSymbol),] #for alt fig

bestImpHigh2 <- cbind(bestImpHigh, basisSource2)
bestImpHigh2$Source <- factor(bestImpHigh2$Source, levels = unique(bestImpHigh2$Source))
bestImpHigh2$Genes <- factor(bestImpHigh2$Genes, levels = rev(bestImpHigh2$Genes))

featBars <- ggplot(bestImpHigh2, aes(x=Genes, y = Overall, fill=Source)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=better_col_vector)+
  guides(fill=guide_legend(ncol = 1)) +
  ylab("Importance") +
  theme(legend.position = 'right',
        legend.text = element_text(size=12),
        axis.text.x = element_blank())


###
pdf("fig/fig_5C.pdf", height = 6, width = 8)
featBars
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
#### figure 5D
#### back to gene expression in surviving and deceased populations
cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_metadata <- "metadata/Glioma_metadata.txt"
cancer_df <- read.table(cancer_files, header=T)

feats <- as.character(bestImpHigh2$Genes)

cancer_df <- subset(cancer_df, GeneSymbol %in% feats)
colnames(cancer_df) <- gsub("\\.", "-", colnames(cancer_df))
colnames(cancer_df) <- make.unique(substring(colnames(cancer_df), 0, 12))
rownames(cancer_df) <- cancer_df$GeneSymbol
cancer_df <- cancer_df[feats,]
cancer_df <- cancer_df[,-1]

cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))

cancer_melt <- data.frame(t(cancer_df), vitalStatus=cancer_meta$vitalstatus, id = colnames(cancer_df))
gorder <- colMeans(cancer_melt[,1:20])
gorder <- names(gorder[order(gorder, decreasing = T)])

plotter <- melt(cancer_melt, id.vars = c("id", "vitalStatus"), measure.vars = colnames(cancer_melt)[1:20])
plotter$logExpr <- log2(plotter$value+1)
plotter$variable <- factor(plotter$variable, levels = gorder)
plotter$Source <- rep(bestImpHigh2$Source, each = 674)
plotter$vitalStatus[plotter$vitalStatus==0] <- "Alive"
plotter$vitalStatus[plotter$vitalStatus==1] <- "Deceased"

figD <- ggplot(plotter, aes(x=variable, y=logExpr,
                            fill=as.factor(vitalStatus))) +
  geom_boxplot(alpha=0.8) +
  scale_fill_manual(values = c("grey50", "grey95")) +
  xlab("Genes") + ylab("Log2 Expression") +
  guides(fill=guide_legend(title="Vital Status")) +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = 'black', hjust=1,angle = 45),
        axis.text.y = element_text(size = 12, color = 'black'))


pdf("fig/fig_5D.pdf", height = 6, width = 8)
figD
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
###### gene gene correlation in top features and diff in idh by clusters?

##### play with IDH mutation status and other variables
cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_metadata <- "metadata/Glioma_metadata.txt"
cancer_df <- read.table(cancer_files, header=T)
rownames(cancer_df) <- cancer_df$GeneSymbol

feats <- as.character(bestImp$Genes)

cancer_full <- subset(cancer_df, GeneSymbol %in% feats)
colnames(cancer_full) <- gsub("\\.", "-", colnames(cancer_full))
colnames(cancer_full) <- make.unique(substring(colnames(cancer_full), 0, 12))
rownames(cancer_full) <- cancer_full$GeneSymbol
cancer_full <- cancer_full[feats,]
cancer_full <- cancer_full[,-1]

## prep metadata for heat
bestImp <- varImp(bestMod)$importance
bestImp$Genes <- rownames(bestImp)
bestImp <- bestImp[order(bestImp$Overall, decreasing = T),]

#bestImpHigh$Genes <- factor(bestImpHigh$Genes, levels = rev(bestImpHigh$Genes))

#merge with cell type + stim colors
basisSource <- read.table("curated_basis12_source.txt",
                          header=T, sep='\t')
rownames(basisSource) <- basisSource$GeneSymbol
basisSource2 <- subset(basisSource, GeneSymbol %in% bestImp$Genes)
#basisSource2 <- basisSource2[as.character(bestImpHigh$Genes),]
bestImp <- bestImp[as.character(basisSource2$GeneSymbol),] #for alt fig

bestImp2 <- cbind(bestImp, basisSource2)
bestImp2$Source <- factor(bestImp2$Source, levels = unique(bestImp2$Source))
bestImp2$Genes <- factor(bestImp2$Genes, levels = rev(bestImp2$Genes))
bestImp2 <- bestImp2[order(bestImp2$Overall, decreasing = T),]

### filter most import by plot(bestImp2$Overall)
### sharp decline after gene #44 so thats the cutoff

bestImp_top <- bestImp2[1:40,]

top_genes <- as.character(bestImp_top$Genes)
### wait are these genes intercorrelated
library(psych)
library(circlize)
cor_mat <- corr.test(t(cancer_full[top_genes,]), ci = F)
rr <- cor_mat$r

ha <- columnAnnotation(df=data.frame(Source=bestImp_top$Source),
                       col=list(Source=c(Macs_IFN_B1a_Macs_IFN_y="#FFDDCD",Macs_IL_10="#FE9A66",
                                         Macs_IL_4_Macs_IL_13 = "#FF5600", 
                                         DCs_IFN_B1a="#9AD4F2",
                                         DCs_IFN_B1a_DCs_IFN_y="#0593E0", DCs_IL_10="#03496F", 
                                         Monos_IFN_B1a_Monos_IFN_y="#A0D9C6", Monos_IL_10="#41B48B",
                                         Monos_IL_4_Monos_IL_13="#0C5037", 
                                         PMNs_IFN_B1a="#EC9FC1",
                                         PMNs_IFN_B1a_PMNs_IFN_y="#D94183",
                                         PMNs_IL_4_PMNs_IL_13="#9B0A4B")))

pdf("meta_plots/top_gene_by_gene_corr_in_TCGA.pdf", height = 5, width = 9)
Heatmap(rr, top_annotation = ha, name="Correlation",
        show_column_names = F, show_row_names = T,
        row_names_gp = gpar(fontsize=9),
        row_names_side = 'left',
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
) #+ va
dev.off()


bestImp_top$cluster <- kmeans(rr, 2)$cluster

top_genes_1 <- as.character(subset(bestImp_top, cluster == "1")$Genes)
top_genes_2 <- as.character(subset(bestImp_top, cluster == "2")$Genes)

#
#### extra meta
### 
cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))

cancer_full <- cbind(t(cancer_full[as.character(bestImp_top$Genes),]), cancer_meta)

cancer_fit <- subset(cancer_full, !is.na(IDH.status))
cancer_fit$test <- as.factor(paste0(cancer_fit$IDH.status,"_", cancer_fit$vitalstatus))
keep_ids <- gsub("-", "\\.", as.character(cancer_fit$id))

### plot gene expression for each cluster
idh_cluster_plot <- melt(cancer_fit, measure.vars = colnames(cancer_fit)[c(1:40)],
                         id.vars = c("id", "IDH.status"))
idh_cluster_plot$cluster <- 0
idh_cluster_plot$cluster[idh_cluster_plot$variable %in% top_genes_1] <- "cluster_1"
idh_cluster_plot$cluster[idh_cluster_plot$variable %in% top_genes_2] <- "cluster_2"

c1 <- subset(idh_cluster_plot, cluster == "cluster_1")
c1_aov = aov(value ~ IDH.status, data=c1)
c1_pval <- summary(c1_aov)[[1]][["Pr(>F)"]][1]

c2 <- subset(idh_cluster_plot, cluster == "cluster_2")
c2_aov = aov(value ~ IDH.status, data=c2)
c2_pval <- summary(c2_aov)[[1]][["Pr(>F)"]][1]

annot_df <- data.frame(pvals = c(c1_pval, c2_pval),
                       y_pos = c(15.4,15.4), x_min = c(0.75,1.75),
                       x_max = c(1.25, 2.25), sym = c("****", "****"))
annot_df$sym <- as.character(annot_df$sym)

g=ggplot(idh_cluster_plot, aes(cluster, log2(value+1), color=IDH.status)) +
  geom_boxplot(alpha=0.6) +
  geom_signif(annotations = c(annot_df$sym), y_position = c(annot_df$y_pos),
              xmin=c(annot_df$x_min), xmax=c(annot_df$x_max),textsize=5, color='black') +
  scale_color_manual(values=c("dodgerblue2", "navy")) +
  #theme_bw() + 
  ylab("Log2 Expression") +
  theme(axis.title.x = element_blank())

pdf("fig/fig_5D_small.pdf", height = 2.6, width = 3.5)
g
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
##### play with IDH mutation status and other variables
cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_metadata <- "metadata/Glioma_metadata.txt"
cancer_df <- read.table(cancer_files, header=T)
rownames(cancer_df) <- cancer_df$GeneSymbol

feats <- as.character(bestImp$Genes)

cancer_full <- subset(cancer_df, GeneSymbol %in% feats)
colnames(cancer_full) <- gsub("\\.", "-", colnames(cancer_full))
colnames(cancer_full) <- make.unique(substring(colnames(cancer_full), 0, 12))
rownames(cancer_full) <- cancer_full$GeneSymbol
cancer_full <- cancer_full[feats,]
cancer_full <- cancer_full[,-1]

## prep metadata for heat
bestImp <- varImp(bestMod)$importance
bestImp$Genes <- rownames(bestImp)
bestImp <- bestImp[order(bestImp$Overall, decreasing = T),]

#bestImpHigh$Genes <- factor(bestImpHigh$Genes, levels = rev(bestImpHigh$Genes))

#merge with cell type + stim colors
basisSource <- read.table("curated_basis12_source.txt",
                          header=T, sep='\t')
rownames(basisSource) <- basisSource$GeneSymbol
basisSource2 <- subset(basisSource, GeneSymbol %in% bestImp$Genes)
#basisSource2 <- basisSource2[as.character(bestImpHigh$Genes),]
bestImp <- bestImp[as.character(basisSource2$GeneSymbol),] #for alt fig

bestImp2 <- cbind(bestImp, basisSource2)
bestImp2$Source <- factor(bestImp2$Source, levels = unique(bestImp2$Source))
bestImp2$Genes <- factor(bestImp2$Genes, levels = rev(bestImp2$Genes))
bestImp2 <- bestImp2[order(bestImp2$Overall, decreasing = T),]

### wait are these genes intercorrelated
library(psych)
library(circlize)
cor_mat <- corr.test(t(cancer_full), ci = F)
rr <- cor_mat$r

ha <- columnAnnotation(df=data.frame(Source=bestImp2$Source),
                       col=list(Source=c(Macs_IFN_B1a_Macs_IFN_y="#FFDDCD",Macs_IL_10="#FE9A66",
                                  Macs_IL_4_Macs_IL_13 = "#FF5600", 
                                  DCs_IFN_B1a="#9AD4F2",
                                  DCs_IFN_B1a_DCs_IFN_y="#0593E0", DCs_IL_10="#03496F", 
                                  Monos_IFN_B1a_Monos_IFN_y="#A0D9C6", Monos_IL_10="#41B48B",
                                  Monos_IL_4_Monos_IL_13="#0C5037", 
                                  PMNs_IFN_B1a="#EC9FC1",
                                  PMNs_IFN_B1a_PMNs_IFN_y="#D94183",
                                  PMNs_IL_4_PMNs_IL_13="#9B0A4B")))
va <- rowAnnotation(df=data.frame(Source=bestImp2$Source),
                       col=list(Source=c(Macs_IFN_B1a_Macs_IFN_y="#FFDDCD",Macs_IL_10="#FE9A66",
                                         Macs_IL_4_Macs_IL_13 = "#FF5600", 
                                         DCs_IFN_B1a="#9AD4F2",
                                         DCs_IFN_B1a_DCs_IFN_y="#0593E0", DCs_IL_10="#03496F", 
                                         Monos_IFN_B1a_Monos_IFN_y="#A0D9C6", Monos_IL_10="#41B48B",
                                         Monos_IL_4_Monos_IL_13="#0C5037", 
                                         PMNs_IFN_B1a="#EC9FC1",
                                         PMNs_IFN_B1a_PMNs_IFN_y="#D94183",
                                         PMNs_IL_4_PMNs_IL_13="#9B0A4B")))

va <- rowAnnotation(df=data.frame(Importance=bestImp2$Overall),
                    col=list(Importance=colorRamp2(c(0, 100), c("white", "navyblue"))))

pdf("meta_plots/gene_by_gene_corr_in_TCGA.pdf", height = 5, width = 10)
Heatmap(rr, top_annotation = ha, name="Correlation",
        show_column_names = F, show_row_names = F,
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        ) + va
dev.off()

##

#

#
#

### 
cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))

cancer_full <- cbind(t(cancer_full), cancer_meta)

cancer_fit <- subset(cancer_full, !is.na(IDH.status))
cancer_fit$test <- as.factor(paste0(cancer_fit$IDH.status,"_", cancer_fit$vitalstatus))
keep_ids <- gsub("-", "\\.", as.character(cancer_fit$id))

##option 1

## make a table of high and low expression and category according to subtypes

impFeats <- subset(bestImp2, Overall > 58.5)

subtypes_all <- data.frame(Gene=NA, Expression=NA, Classical=NA, G.CIMP=NA,
                           IDHmut.codel=NA, IDHmut.non.codel=NA, IDHwt=NA, 
                           Mesenchymal=NA, Neural=NA, Proneural=NA)
mol_marks_all <- data.frame(Gene=NA, Expression=NA, IDH_wt=NA, IDH_mutant=NA, 
                            codel=NA, non_codel=NA, Methylated=NA, Unmethylated=NA)
for(i in 1:length(unique(impFeats$Genes))){
  gene <- as.character(unique(impFeats$Genes[i]))
  
  q3 <- summary(cancer_fit[,gene])[5]
  q1 <- summary(cancer_fit[,gene])[2]
  med <- summary(cancer_fit[,gene])[3]
  
  new_meta <- data.frame(id=cancer_fit$id, Expr=cancer_fit[,gene])
  new_meta$Expr_class <- "mid"
  new_meta$Expr_class[new_meta$Expr > q3] <- 'high'
  new_meta$Expr_class[new_meta$Expr < q1] <- 'low'
  
  #idh
  idh_tab <- table(new_meta$Expr_class, cancer_fit$IDH.status)
  idh_df <- data.frame(IDH_wt=idh_tab[1:2,2],
                       IDH_mutant = idh_tab[1:2,1])
  idh_df <- round((idh_df/rowSums(idh_df))*100,2)
  
  #codel
  codel_tab <- table(new_meta$Expr_class, cancer_fit$X1p.19q.codeletion)
  codel_df <- data.frame(codel=codel_tab[1:2,1],
                       non_codel = codel_tab[1:2,2])
  codel_df <- round((codel_df/rowSums(codel_df))*100,2)
  
  #MGMT promoter
  mgmt_tab <- table(new_meta$Expr_class, cancer_fit$MGMT.promoter.status)
  mgmt_df <- data.frame(Methylated=mgmt_tab[1:2,1],
                         Unmethylated = mgmt_tab[1:2,2])
  mgmt_df <- round((mgmt_df/rowSums(mgmt_df))*100,2)
  
  mol_marks <- data.frame(Gene=rep(gene,2), Expression=c("High", "Low"), idh_df, codel_df, mgmt_df)
  mol_marks_all <- rbind(mol_marks_all, mol_marks)
  
  #subtypes
  subtypes_tab <- table(new_meta$Expr_class, cancer_fit$Original.Subtype)
  subtypes_df <- t(data.frame(upper=subtypes_tab[1,],
                        lower = subtypes_tab[2,]))
  subtypes_df <- round((subtypes_df/rowSums(subtypes_df))*100,6)

  subtypes_marks <- data.frame(Gene=rep(gene,2), Expression=c("High", "Low"), subtypes_df)
  subtypes_all <- rbind(subtypes_all, subtypes_marks)
}
subtypes_all<- subtypes_all[-1,]
mol_marks_all<- mol_marks_all[-1,]

rownames(subtypes_all) <- paste0(subtypes_all$Gene, "_", subtypes_all$Expression)
subtype_heat <- t(subtypes_all[,3:10])
colnames(subtype_heat) <- gsub("_.*", "", colnames(subtype_heat))

###quick ggplots
g_wt = ggplot(mol_marks_all, aes(x=Expression, y=IDH_wt, fill=Expression)) + 
  geom_col(position='dodge') + ylab("Percentage of Samples with IDH wild type") +
  scale_fill_manual(values = c("navyblue", "dodgerblue2"))
pdf("meta_plots/idh_wt_association.pdf", height = 7, width = 12)
g_wt + facet_wrap(~Gene)
dev.off()

g_mut = ggplot(mol_marks_all, aes(x=Expression, y=IDH_mutant, fill=Expression)) + 
  geom_col(position='dodge') + ylab("Percentage of Samples with IDH mutation") +
  scale_fill_manual(values = c("navyblue", "dodgerblue2"))
pdf("meta_plots/idh_mut_association.pdf", height = 7, width = 12)
g_mut + facet_wrap(~Gene)
dev.off()
##

ha <- columnAnnotation(df=data.frame(Expression=gsub(".*_", "", rownames(subtypes_all))
                                     ),
                       col = list(Expression=c(High='black', Low="grey50")))

ll <- columnAnnotation(labels = anno_text(colnames(subtype_heat), gp=gpar(fontsize=8),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(1.2,"cm")), 
                       height = unit(1,"cm"))

pdf("meta_plots/subtypes_plus_exp.pdf", height = 4, width = 10)
Heatmap(subtype_heat, name="Percentage", 
        top_annotation = ha, bottom_annotation = ll,
        bottom_annotation_height = unit(1.2,"cm"),
        show_column_names = F,
        col = colorRamp2(c(0, 80), c("darkslategray2", "navyblue")))
dev.off()

##mol marks
rownames(mol_marks_all) <- paste0(mol_marks_all$Gene, "_", mol_marks_all$Expression)
mol_marks_heat <- t(mol_marks_all[,3:8])
colnames(mol_marks_heat) <- gsub("_.*", "", colnames(mol_marks_heat))

ha <- columnAnnotation(df=data.frame(Expression=gsub(".*_", "", rownames(mol_marks_all))
),
col = list(Expression=c(High='black', Low="grey50")))

ll <- columnAnnotation(labels = anno_text(colnames(mol_marks_heat), gp=gpar(fontsize=8),
                                          which = "column", rot=60, 
                                          just = 'right', offset=unit(1.2,"cm")), 
                       height = unit(1,"cm"))

pdf("meta_plots/mol_marks_plus_exp.pdf", height = 4, width = 10)
Heatmap(mol_marks_heat, name="Percentage", 
        #split=c("IDH", "IDH", "codel", "codel", "Methylation", "Methylation"),
        top_annotation = ha, bottom_annotation = ll,
        bottom_annotation_height = unit(1.2,"cm"),
        show_column_names = F,
        col = colorRamp2(c(0, 80), c("darkslategray2", "navyblue")))
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


### option 2 DE genes between groups

library(edgeR)
design <- model.matrix(~test, cancer_fit)
edgeR_cds = DGEList(cancer_df[,keep_ids])
edgeR_cds = calcNormFactors( edgeR_cds )
edgeR_cds = estimateGLMCommonDisp(edgeR_cds, design)

contrast_list <- makeContrasts(Mutant_0 - Mutant_1,
                               Mutant_0 - WT_0,
                               Mutant_0 - WT_1,
                               Mutant_1 - WT_0,
                               Mutant_1 - WT_1,
                               WT_0 - WT_1,
                               levels = levels(cancer_fit$test))

fit <- glmFit(edgeR_cds, design)
res_df <- data.frame(logFC=NA, logCPM=NA, LR=NA, PValue=NA, padj=NA, gene=NA, comp=NA)
for(i in 1:6){
res <- glmLRT(fit, contrast=contrast_list[,i])$table
pval = res$PValue
padj = p.adjust( pval, method="BH")
res = data.frame(res, padj=padj, gene=rownames(res))
res$comp <- rep(colnames(contrast_list)[i],nrow(res))
res_df <- rbind(res_df, res)
}
res_df <- res_df[-1,]

glist=list()
for(i in 1:20){
  curr_gene = colnames(cancer_fit)[i]

  pvals <- subset(res_df, gene == curr_gene)$padj
  gmax <- max(log2(cancer_fit[,curr_gene]+1))
  annot_df <- data.frame(x_min = c(1,1,1,2,2,3), x_max = c(2,3,4,3,4,4), 
                         y_pos = c(gmax+1, gmax+3,gmax+5, gmax+2,gmax+4,gmax+1),
                         pvals = pvals)
  annot_df$sym <- annot_df$pvals
  annot_df$sym[annot_df$sym<0.05]<-"*"
  annot_df$sym[annot_df$sym<0.01]<-"**"
  annot_df$sym[annot_df$sym<0.001]<-"**"
  annot_df$sym[annot_df$sym>0.05]<-'n.s.'
  new_labs <- c("Mutant_0" = "Mutant\nAlive\nn=348",
                "Mutant_1" = "Mutant\nDeceased\nn=81",
                "WT_0" = "Mutant\nAlive\nn=70",
                "WT_1" = "Mutant\nAlive\nn=166")
  
  gg = ggplot(cancer_fit, aes(test, log2(cancer_fit[,i]+1), fill=test))+
    geom_boxplot(alpha=0.2, outlier.alpha = 0) +
    geom_jitter(width = 0.2, aes(color=test)) +
    scale_x_discrete(labels=new_labs)+
    ylab("Log2 Expression") + ggtitle(curr_gene) +
    theme(legend.position = 'none')
  glist[[i]]<-ggplotGrob(gg + geom_signif(annotations = c(annot_df$sym), y_position = c(annot_df$y_pos),
                               xmin=c(annot_df$x_min), xmax=c(annot_df$x_max),textsize=5))
}

library(cowplot)

pdf("meta_plots/idh_plots.pdf", height=25, width = 18)
do.call("grid.arrange", c(glist, ncol=4))
dev.off()

asf1b_idh <- glist[[18]]
grin3a_idh <- glist[[8]]

#
#
#
####
library(edgeR)
design <- model.matrix(~IDH.status, cancer_fit)
edgeR_cds = DGEList(cancer_df[,keep_ids])
edgeR_cds = calcNormFactors( edgeR_cds )
edgeR_cds = estimateGLMCommonDisp(edgeR_cds, design)

contrast_list <- makeContrasts(Mutant - WT,
                               levels = levels(cancer_fit$IDH.status))

fit <- glmFit(edgeR_cds, design)
res_df <- data.frame(logFC=NA, logCPM=NA, LR=NA, PValue=NA, padj=NA, gene=NA, comp=NA)
for(i in 1:1){
  res <- glmLRT(fit, contrast=contrast_list[,i])$table
  pval = res$PValue
  padj = p.adjust( pval, method="BH")
  res = data.frame(res, padj=padj, gene=rownames(res))
  res$comp <- rep(colnames(contrast_list)[i],nrow(res))
  res_df <- rbind(res_df, res)
}
res_df <- res_df[-1,]

glist=list()
for(i in 1:20){
  curr_gene = colnames(cancer_fit)[i]
  
  pvals <- subset(res_df, gene == curr_gene)$padj
  gmax <- max(log2(cancer_fit[,curr_gene]+1))
  annot_df <- data.frame(x_min = c(1), x_max = c(2), 
                         y_pos = c(gmax+0.01),
                         pvals = pvals)
  annot_df$sym <- annot_df$pvals
  annot_df$sym[annot_df$sym<0.05]<-"*"
  annot_df$sym[annot_df$sym<0.01]<-"**"
  annot_df$sym[annot_df$sym<0.001]<-"***"
  annot_df$sym[annot_df$sym>0.05]<-'n.s.'
  new_labs <- c("Mutant" = "Mutant\n(n=429)",
                "WT" = "WT\n(n=236)")
  
  gg = ggplot(cancer_fit, aes(IDH.status, log2(cancer_fit[,i]+1), fill=IDH.status))+
    geom_boxplot(alpha=0.2, outlier.alpha = 0) +
    geom_jitter(width = 0.2, aes(color=IDH.status)) +
    scale_x_discrete(labels=new_labs)+
    ylab("Log2 Expression") + ggtitle(curr_gene) +
    theme(legend.position = 'none')
  glist[[i]]<-ggplotGrob(gg + geom_signif(annotations = c(annot_df$sym), y_position = c(annot_df$y_pos),
                                          xmin=c(annot_df$x_min), xmax=c(annot_df$x_max),textsize=5))
}

library(cowplot)

pdf("meta_plots/idh_no_surv_plots.pdf", height=20, width = 18)
do.call("grid.arrange", c(glist, ncol=5))
dev.off()

##

#
## subtype?

cancer_fit <- subset(cancer_full, !is.na(Original.Subtype))
cancer_fit$test <- as.factor(make.names(paste0(cancer_fit$Original.Subtype,"_", cancer_fit$vitalstatus)))
keep_ids <- gsub("-", "\\.", as.character(cancer_fit$id))

library(edgeR)
design <- model.matrix(~test, cancer_fit)
edgeR_cds = DGEList(cancer_df[,keep_ids])
edgeR_cds = calcNormFactors( edgeR_cds )
edgeR_cds = estimateGLMCommonDisp(edgeR_cds, design)

contrast_list <- makeContrasts(Classical_0 - Classical_1,
                               G.CIMP_0 - G.CIMP_1,
                               IDHmut.codel_0 - IDHmut.codel_1,
                               IDHmut.non.codel_0 - IDHmut.non.codel_1,
                               IDHwt_0 - IDHwt_1,
                               Mesenchymal_0 - Mesenchymal_1,
                               Neural_0 - Neural_1,
                               Proneural_0 - Proneural_1,
                               levels = levels(cancer_fit$test))

fit <- glmFit(edgeR_cds, design)
res_df <- data.frame(logFC=NA, logCPM=NA, LR=NA, PValue=NA, padj=NA, gene=NA, comp=NA)
for(i in 1:8){
  res <- glmLRT(fit, contrast=contrast_list[,i])$table
  pval = res$PValue
  padj = p.adjust( pval, method="BH")
  res = data.frame(res, padj=padj, gene=rownames(res))
  res$comp <- rep(colnames(contrast_list)[i],nrow(res))
  res_df <- rbind(res_df, res)
}
res_df <- res_df[-1,]

glist=list()
for(i in 1:20){
  curr_gene = colnames(cancer_fit)[i]
  
  pvals <- subset(res_df, gene == curr_gene)$padj
  gmax <- max(log2(cancer_fit[,curr_gene]+1))
  annot_df <- data.frame(x_min = c(1,3,5,7,9,11,13,15), x_max = c(2,4,6,8,10,12,14,16), 
                         y_pos = c(gmax+0.6,gmax+0.6,gmax+0.6,gmax+0.6,gmax+0.6,gmax+0.6,gmax+0.6,gmax+0.6),
                         pvals = pvals)
  annot_df$sym <- annot_df$pvals
  annot_df$sym[annot_df$sym<0.05]<-"*"
  annot_df$sym[annot_df$sym<0.01]<-"**"
  annot_df$sym[annot_df$sym<0.001]<-"**"
  annot_df$sym[annot_df$sym>0.05]<-'n.s.'
  new_labs <- c("Classical_0" = "Alive\nn=8",
                "Classical_1" = "Deceased\nn=32",
                "G.CIMP_0" = "Alive\nn=5",
                "G.CIMP_1" = "Deceased\nn=3",
                "IDHmut.codel_0" = "Alive\nn=147",
                "IDHmut.codel_1" = "Deceased\nn=22",
                "IDHmut.non.codel_0" = "Alive\nn=197",
                "IDHmut.non.codel_1" = "Deceased\nn=51",
                "IDHwt_0" = "Alive\nn=44",
                "IDHwt_1" = "Deceased\nn=51",
                "Mesenchymal_0" = "Alive\nn=11",
                "Mesenchymal_1" = "Deceased\nn=40",
                "Neural_0" = "Alive\nn=3",
                "Neural_1" = "Deceased\nn=25",
                "Proneural_0" = "Alive\nn=3",
                "Proneural_1" = "Deceased\nn=27")
  
  plotter <- cancer_fit
  plotter$vitalstatus[plotter$vitalstatus==0] <- "Alive"
  plotter$vitalstatus[plotter$vitalstatus==1] <- "Deceased"
  
  liner <- data.frame(xpos=seq(2.5,14.5,by=2))
  texter <- data.frame(types = levels(cancer_fit$Original.Subtype), 
                       xpos = seq(1.5,15.5,by=2), ypos = rep(gmax+2,8))
  
  gg = ggplot(plotter, aes(test, log2(plotter[,i]+1)))+
    geom_boxplot(alpha=0.2, outlier.alpha = 0, aes(fill=vitalstatus)) +
    geom_jitter(width = 0.2, aes(color=vitalstatus)) +
    geom_vline(data=liner, mapping=aes(xintercept=xpos), color="grey50", alpha=0.5, linetype=2)+
    geom_text(data=texter, mapping=aes(x=xpos, y=ypos, label=types))+
    scale_x_discrete(labels=new_labs)+
    ylab("Log2 Expression") + ggtitle(curr_gene) +
    theme(legend.position = 'none')
  glist[[i]]<-ggplotGrob(gg + geom_signif(annotations = c(annot_df$sym), y_position = c(annot_df$y_pos),
                                          xmin=c(annot_df$x_min), xmax=c(annot_df$x_max),textsize=5))
}

library(cowplot)

pdf("meta_plots/subtype_plots.pdf", height=35, width = 25)
do.call("grid.arrange", c(glist, ncol=2, nrow=10))
dev.off()

asf1b_subtype <- glist[[18]]
grin3a_subtype <- glist[[8]]

##
#
#
## what about a big ole heatmap of expression levels and metadata for all glioma?


cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_metadata <- "metadata/Glioma_metadata.txt"
cancer_df <- read.table(cancer_files, header=T)
rownames(cancer_df) <- cancer_df$GeneSymbol

##what genes do u want

### load best model

bestMod <- readRDS("curated_results/lasso_model_7_5yrs_min.max_norm.rds")
bestModParam <- unlist(bestMod$bestTune)
selectedIndices <- bestMod$pred$mtry == bestModParam
huh <- bestMod$pred[selectedIndices, ]


# Best features bars
bestImp <- varImp(bestMod)$importance
bestImp$Genes <- rownames(bestImp)
bestImp <- bestImp[order(bestImp$Overall, decreasing = T),]
bestImpHigh <- bestImp[1:50,]
#bestImpHigh$Genes <- factor(bestImpHigh$Genes, levels = rev(bestImpHigh$Genes))

#merge with cell type + stim colors
basisSource <- read.table("curated_basis12_source.txt",
                          header=T, sep='\t')
rownames(basisSource) <- basisSource$GeneSymbol
basisSource2 <- subset(basisSource, GeneSymbol %in% bestImpHigh$Genes)
#basisSource2 <- basisSource2[as.character(bestImpHigh$Genes),]
bestImpHigh <- bestImpHigh[as.character(basisSource2$GeneSymbol),] #for alt fig

bestImpHigh2 <- cbind(bestImpHigh, basisSource2)
bestImpHigh2$Source <- factor(bestImpHigh2$Source, levels = unique(bestImpHigh2$Source))
bestImpHigh2$Genes <- factor(bestImpHigh2$Genes, levels = rev(bestImpHigh2$Genes))

#
#
#
feats <- as.character(bestImpHigh2$Genes)
#feats <- basisSource$GeneSymbol

cancer_full <- subset(cancer_df, GeneSymbol %in% feats)
colnames(cancer_full) <- gsub("\\.", "-", colnames(cancer_full))
colnames(cancer_full) <- make.unique(substring(colnames(cancer_full), 0, 12))
rownames(cancer_full) <- cancer_full$GeneSymbol
cancer_full <- cancer_full[feats,]
cancer_full <- cancer_full[,-1]

cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))

library(ComplexHeatmap)
library(circlize)

cscale <- t(scale(t(cancer_full)))
keep_rows <- which(!is.na(cscale[,1]))

ha <- columnAnnotation(df=cancer_meta[,c(3,29,38,57,59,61)])

full_heat_glioma <- Heatmap(cscale[keep_rows,],
                         top_annotation = ha,
                         heatmap_legend_param = list(title = "Scaled Expression"),#,
                                                     #title_gp = gpar(fontsize = 15),
                                                     #legend_height = unit(3, "cm"),
                                                     #labels_gp = gpar(fontsize = 14)),
                         #split = n_genes$Source,
                         col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                         cluster_rows = T, cluster_columns = T,
                         column_title = "",
                         column_title_gp = gpar(fontsize = 15),
                         show_column_names = F, show_row_names = T,
                         row_names_gp = gpar(fontsize = 8))

pdf("meta_plots/full_meta_heat_GeneExp_top50.pdf", height = 8, width = 10)
full_heat_glioma
dev.off()

#
#
#
###
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
#### Fig 5 E-H
### Survival curves by gene for TCGA and Precog!

load("ind_test/surv_test_GSE_quart.RData")
precog_curves <- pp
load("ind_test/surv_test_TCGA_quart.RData")
tcga_curves <- pp


asf1b_p <- precog_curves[[12]]
asf1b_t <- tcga_curves[[102]]

rbbp6_p <- precog_curves[[22]]
rbbp6_t <- tcga_curves[[54]]

plscr1_t <- tcga_curves[[55]]
slc1a4_t <- tcga_curves[[25]]
grin3a_t <- tcga_curves[[47]]
eri1_t <- tcga_curves[[107]]
eepd1 <- tcga_curves[[21]]

apobec_t <- tcga_curves[[37]]
apobec_p <- precog_curves[[24]]

ptgir_t <- tcga_curves[[78]]
ptgir_p <- precog_curves[[19]]

figE_H <- arrange_ggsurvplots(list(asf1b_t, slc1a4_t, 
                                   plscr1_t, grin3a_t,
                                   eri1_t, eepd1), 
                            print = F, nrow = 2, ncol = 3)


pdf("fig/fig_5E_H.pdf", height = 8, width = 12)
figE_H
dev.off()


supp <- arrange_ggsurvplots(list(asf1b_t, rbbp6_t, apobec_t, ptgir_t,
                                 asf1b_p, rbbp6_p, apobec_p, ptgir_p), 
                              print = F, nrow = 4, ncol = 2)
pdf("suppl_curves.pdf", height = 15, width = 10)
supp
dev.off()



#
#


