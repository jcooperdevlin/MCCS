
##############

############## Building a Machine Learning Model to predict cancer survival from our & other gene signatures

############## 2.16.19

args = commandArgs(trailingOnly=TRUE)

library(caret)
#library(preprocessCore)
#library(data.table)
#library(reshape)
#library(circlize)
#library(TCGA2STAT)
library(ROCR)
library(plotROC)
library(randomForest)
library(ComplexHeatmap)
library(circlize)

normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}


setwd("/Users/devlij03/Google Drive/Png/cooper/PNAS_manuscript/int/glioma/")

#
####### for each cancer

###### cancer list?
#cancer_files <- list.files("genesymbols", full.names = T)
#cancer_files <- "genesymbols/Glioma_genesymbols.txt"
cancer_files <- args[1]
#cancer_metadata <- list.files("metadata", full.names=T)
cancer_metadata <- args[2]
#cancer_metadata <- "metadata/Glioma_metadata.txt"
sig_matrix <- args[3]
#sig_matrix <- c("curated_basis12.txt")
out_dir <- args[4]
CV <- 7
surv_num <- as.numeric(args[5])
f = 2

trains_list <- c("glioma_train_list2yr.txt","glioma_train_list5yr.txt","glioma_train_list10yr.txt")
tests_list <- c("glioma_test_list2yr.txt","glioma_test_list5yr.txt","glioma_test_list10yr.txt")

results <- data.frame(cancer_type = NA, test = NA, norm = NA, Survival = NA,
                      trainSetSize = NA, TrainAccuracy = NA, TrainAUC=NA,
                      testSetSize=NA, truePos = NA, falsePos = NA, falseNeg = NA, trueNeg = NA,
                      TestAccuracy = NA, Sensitivity = NA, Specificity = NA, Precision = NA,
                      Recall = NA, BalancedAccuracy = NA, TestAUC = NA)

current_cancer <- gsub("_genesymbols\\.txt", "", basename(cancer_files))
message <- paste0("Working on ", current_cancer, "...")
print(message)

### load cancer data and subset for genes in our signature!
cancer_df <- read.table(cancer_files, header=T)

X <- read.table(sig_matrix,header=T,sep="\t")
basis <- unique(as.character(X$GeneSymbol))

cancer_df <- subset(cancer_df, GeneSymbol %in% basis)
colnames(cancer_df) <- gsub("\\.", "-", colnames(cancer_df))
colnames(cancer_df) <- make.unique(substring(colnames(cancer_df), 0, 12))
rownames(cancer_df) <- cancer_df$GeneSymbol
cancer_df <- cancer_df[,-1]

cancer_meta <- read.table(cancer_metadata, header=T, sep='\t')
rownames(cancer_meta) <- make.unique(substring(cancer_meta$id, 0, 12))

cancer_meta$orig_vital <- cancer_meta$vitalstatus
new_vital <- cancer_meta$vitalstatus
new_vital[cancer_meta$daystodeath > 365*surv_num] <- 0
cancer_meta$newvital <- new_vital
###
# fix class imbalance
#cancer0 <- subset(cancer_meta, new_vital == 0)
#cancer1 <- subset(cancer_meta, new_vital == 1)

#keep0 <- as.character(sample(cancer0$id, nrow(cancer1)))
#cancer0_small <- cancer0[keep0,]

#cancer_meta_new <- rbind(cancer0_small, cancer1)
#cancer_meta_new <- cancer_meta_new[sample(1:nrow(cancer_meta_new)),]

#cancer_df_new <- cancer_df[,rownames(cancer_meta_new)]

### perfect class balance
## or not
cancer_meta <- cancer_meta
cancer_df <- cancer_df


feat_heat <- data.frame(Genes = rownames(cancer_df))
glist <- list()
counter = 1


##
# split into test and train 80/20

#test_group <- as.character(sample(cancer_meta$id, ceiling(nrow(cancer_meta)*0.2)))
#train_group <- colnames(cancer_df)[!colnames(cancer_df) %in% test_group]

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


#
#
#
# Let's machine learn!!
#feat_heat <- data.frame(Genes = rownames(data_train_list[[1]]))
#glist <- list()
#counter = 1
norm_names <- c("none", "min.max_norm")

data.traink <- data_train_list[[f]]
train_class <- make.names(train_class)
data_run1 <- data.frame(t(data.traink), condition=train_class)
data_test1 <- data_test_list[[f]]
test_class <- make.names(test_class)
data_test1 <- data.frame(t(data_test1), condition=test_class)


### lasso
# add lasso

added_ctrl <- trainControl(method = "repeatedcv", number = CV,
                           repeats = 10, verboseIter = F,classProbs = F,
                           savePredictions = T)

message <- paste0("Working on ", current_cancer, "...", "Lasso")
print(message)

rlabels=as.vector(data_run1$condition)
rlabels <- gsub("X0", 0, rlabels)
rlabels <- gsub("X1", 1, rlabels)
rlabels <- as.numeric(rlabels)
data_run1$condition <- rlabels

tlabels=as.vector(data_test1$condition)
tlabels <- gsub("X0", 0, tlabels)
tlabels <- gsub("X1", 1, tlabels)
tlabels <- as.numeric(tlabels)
data_test1$condition <- tlabels

lasso_test <- train(condition ~ ., 
                    data = data_run1, tuneLength = 50,
                    method = "lasso",
                    trControl = added_ctrl)

saveRDS(lasso_test, paste0(out_dir, "lasso_model_", CV,"_", surv_num, "yrs_", norm_names[f], ".rds"))

## features
lasso_imp <- data.frame(Overall=varImp(lasso_test)$importance$Overall)
feat_heat <- cbind(feat_heat, lasso_imp)

#gather results
#train results
bestModParam <- unlist(lasso_test$bestTune)
selectedIndices <- lasso_test$pred$fraction == bestModParam
huh <- lasso_test$pred[selectedIndices, ]

g <- ggplot(huh, aes(m=pred, d=obs)) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc() +
  ggtitle(paste0("lasso_model_", CV, "_", norm_names[f], "_", surv_num))
glist[[counter]] <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
counter = counter+1

cv_auc <- round((calc_auc(g))$AUC, 4)

namer <- paste0(out_dir, "Lasso_Train_", CV, "_", surv_num, "_yrs_", norm_names[f], "_AUC.txt")
write.table(huh, namer, sep='\t', row.names = F, quote = F)
#

preds <- as.numeric(as.character(predict(lasso_test, newdata = data_run1)))
preds[preds>0.5]<-1
preds[preds<0.5]<-0
cf <- confusionMatrix(data = factor(preds, levels=c(0,1)), 
                      reference = as.factor(data_run1$condition), 
                      mode = "prec_recall")

tacc <- cf$overall[1]
names(tacc) <- "TrainAccuracy"
trainSetSize = length(train_class)
train_stats <- c(trainSetSize, tacc, TrainAUC = cv_auc)
names(train_stats)[1] <- "trainSetSize"


#test results
preds1 <- data.frame(X1=predict(lasso_test, newdata = data_test1))
#preds1$obs <- predict(lasso_test, newdata = data_test1)
preds1$true <- data_test1$condition
preds1$true <- gsub(0, "X0",preds1$true)
preds1$true <- gsub(1, "X1",preds1$true)

g <- ggplot(preds1, aes(m=X1, d=factor(true, levels = c("X0", "X1")))) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc() +
  ggtitle(paste0("Lasso_Test_", CV, "_", norm_names[f], "_", surv_num))
glist[[counter]] <- g + annotate("text", x=0.75, y=0.25, label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
counter = counter+1

namer <- paste0(out_dir, "Lasso_Test_", CV, "_", surv_num, "_yrs_", norm_names[f], "_AUC.txt")
write.table(preds1, namer, sep='\t', row.names = F, quote = F)

preds <- as.numeric(as.character(predict(lasso_test, newdata = data_test1)))
preds[preds>0.5]<-1
preds[preds<0.5]<-0
cf <- confusionMatrix(data = factor(preds, levels=c(0,1)), 
                      reference = as.factor(data_test1$condition), 
                      mode = "prec_recall")

c_mat <- c(truePos=cf$table[1,1], falsePos=cf$table[2,1], 
           falseNeg=cf$table[1,2], trueNeg=cf$table[2,2])

accuracies <- c(cf$overall[1], cf$byClass[1], cf$byClass[2], 
                cf$byClass[5], cf$byClass[6], cf$byClass[11])
names(accuracies)[1] <- "TestAccuracy"
names(accuracies)[6] <- "BalancedAccuracy"
predictions=as.vector(preds)
predictions <- gsub("X0", 0, predictions)
predictions <- gsub("X1", 1, predictions)
predictions <- as.numeric(predictions)

labels=as.vector(data_test1$condition)
labels <- gsub("X0", 0, labels)
labels <- gsub("X1", 1, labels)
labels <- as.numeric(labels)

pred=prediction(predictions,labels)
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC=perf_AUC@y.values[[1]]

testSetSize = length(test_class)
test_stats <- c(testSetSize, c_mat, accuracies, TestAUC=AUC)
names(test_stats)[1] <- "testSetSize"

results_add <- cbind(data.frame(cancer_type = current_cancer, test = "lasso", 
                                norm = norm_names[f], Survival = surv_num),
                     t(data.frame(train_stats)), t(data.frame(test_stats)))
results <- rbind(results, results_add)
  

message <- paste0(current_cancer, "...Complete!")
print(message)


results <- results[-1,]
namer <- gsub("_genesymbols.txt", "", paste0(out_dir, 
                                             current_cancer,
                                             "_", CV, "_", surv_num, "yrs_", norm_names[f], "_ml_results.txt"))
write.table(results, namer, sep='\t', quote=F, row.names=F)

library(cowplot)

namer <- gsub("_genesymbols.txt", "", paste0(out_dir, 
                                             current_cancer,
                                             "_", CV, "_", surv_num, "yrs_", norm_names[f], "_ml_results.pdf"))
pdf(namer, height=8, width = 15)
pp = plot_grid(plotlist = glist, ncol = 4)
print(pp)
dev.off()

namer <- gsub("_genesymbols.txt", "", paste0(out_dir, 
                                             current_cancer,
                                             "_", CV, surv_num, "yrs_", norm_names[f], "_plots.RDS"))
saveRDS(glist, file = namer)
