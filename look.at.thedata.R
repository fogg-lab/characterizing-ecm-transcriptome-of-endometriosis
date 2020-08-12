library(tidyverse);
##  2020-08-09

## Look at the DE analysis results 

obj = read_tsv("../analysis/unified_cervical_data_unfiltered_DESeq_results.tsv")

## Read in data
setwd("D:/Box Sync/Fogg-Lab-RNA-Seq-analysis/R");

counts = read_tsv("../data/unified_cervical_data/counts.tsv") %>%  select(-"Entrez_Gene_Id") %>% column_to_rownames (var = "Hugo_Symbol");
y = as.matrix(counts);

## Look at genes with greatest log2FoldChange
or = order(obj$log2FoldChange, decreasing = TRUE)
top = or[1:10]

obj[top,]

par(mfrow=c(2,5))
for (i in top) {
  print(y[i, 12:13]);
  boxplot(log(y[i,]) ~ coldata$condition * coldata$data_source)
}


## Look at genes with least p-values
or = order(obj$pvalue)
top = or[1:10]
top = or[11:20]

obj[top,]

par(mfrow=c(2,5))
for (i in top) {
  print(y[i, 12:13]);
  boxplot(log(y[i,]) ~ coldata$condition * coldata$data_source)
}



 