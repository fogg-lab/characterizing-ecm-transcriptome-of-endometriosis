**Description:** Mean plots to compare expression of 49 genes (from the top 5 KEGG pathways) between mild and severe endometriosis.

The data files (`l2fc_{stage}_{direction_of_dysregulation}.tsv`) contain the $log_2$ fold change values to plot. Only the first column is actual data (the means), the other two are mean +/- std. The figure generation code uses these two dummy values around the mean infer what the standard deviation is. This is done for programming convenience.
