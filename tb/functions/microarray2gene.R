#### coverter from micro probes to gene symbols ready for cibersort

#source("https://bioconductor.org/biocLite.R")
#biocLite("illuminaHumanv4.db")
library("illuminaHumanv4.db")

microarray2gene <- function(input_df) {
  input_df <- subset(input_df, ID_REF != "")
  rownames(input_df) <- input_df$ID_REF
  
  probeID <- input_df$ID_REF
  
  new_ids <- data.frame(Gene=unlist(mget(x = probeID,envir = illuminaHumanv4SYMBOL, ifnotfound = NA)))
  new_ids <- subset(new_ids, !is.na(Gene))
  
  input_df <- input_df[rownames(new_ids),]
  GeneSymbol <- new_ids$Gene

  input_df_agg <- aggregate(input_df[,-1], list(GeneSymbol), mean)
  colnames(input_df_agg)[1] <- "GeneSymbol"
  
  return(input_df_agg)
}