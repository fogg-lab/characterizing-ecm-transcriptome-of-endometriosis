library(org.Hs.eg.db)
library(ggplot2)
library(readr)
library(car)

# Load clinical and expression data
clinical <- read.csv("../data_prep/metadata/stagewise_coldata.tsv", sep = "\t")
expression <- read.csv("../data/all/all_phases_all_matrisome_counts.tsv", sep = "\t", row.names = 1)

# Load DEG lists
deg_proliferative <- scan("../DEMGs/stage_significant_demgs/proliferative_demg.txt", what = "", sep = "\n")
deg_secretory_early <- scan("../DEMGs/stage_significant_demgs/secretory_early_demg.txt", what = "", sep = "\n")
deg_secretory_mid <- scan("../DEMGs/stage_significant_demgs/secretory_mid_demg.txt", what = "", sep = "\n")

# Load enrichment results
go_results <- read.csv("../enrichment_analysis/stage_significant_results/go_simplified_results.tsv", sep = "\t")
kegg_results <- read.csv("../enrichment_analysis/stage_significant_results/kegg_results.tsv", sep = "\t")

# Function to perform ANOVA and return significant genes
perform_anova <- function(expression, clinical, deg_list, phase, go_results, kegg_results) {
  phase_samples <- clinical$sample_name[clinical$phase == phase]
  phase_expression <- expression[, phase_samples]
  phase_condition <- clinical$condition[clinical$phase == phase]

  genes_and_pvalues <- sapply(deg_list, function(gene) {
    if(gene %in% rownames(phase_expression)){
      gene_expression <- data.frame(Expression = as.numeric(phase_expression[gene, ]), Condition = phase_condition)
      aov_results <- aov(Expression ~ Condition, data = gene_expression)
      p_value <- summary(aov_results)[[1]][["Pr(>F)"]][1]

      if (!is.na(p_value) && p_value < 0.05) {
        return(p_value)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })

    # Remove NA p-values and genes
    genes_and_pvalues <- genes_and_pvalues[!is.na(genes_and_pvalues)]

    # Order the genes by p-values
    genes_and_pvalues <- genes_and_pvalues[order(genes_and_pvalues)]

    # Extract genes from top 5 GO categories
    top_go_genes <- unlist(strsplit(as.character(head(go_results$geneID, 5)), split = "/"))

    # Extract genes from top 5 KEGG categories
    top_kegg_entrez <- unlist(strsplit(as.character(head(kegg_results$geneID, 5)), split = "/"))

    # Convert Entrez IDs to gene symbols
    top_kegg_genes <- mapIds(org.Hs.eg.db, keys = top_kegg_entrez, column = "SYMBOL", keytype = "ENTREZID")

    # Filter significant genes to only include those in top GO or KEGG categories
    go_significant_genes <- genes_and_pvalues[names(genes_and_pvalues) %in% top_go_genes]
    kegg_significant_genes <- genes_and_pvalues[names(genes_and_pvalues) %in% top_kegg_genes]

    go_significant_genes <- data.frame(symbol = names(go_significant_genes), p_value = go_significant_genes, row.names = NULL)
    kegg_significant_genes <- data.frame(symbol = names(kegg_significant_genes), p_value = kegg_significant_genes, row.names = NULL)

    return(list(go = go_significant_genes, kegg = kegg_significant_genes))
}

# Perform ANOVA and save results
results_proliferative <- perform_anova(expression, clinical, deg_proliferative, "proliferative", go_results, kegg_results)
write.table(results_proliferative$go, file = "./results_proliferative_go.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_proliferative$kegg, file = "./results_proliferative_kegg.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

results_secretory_early <- perform_anova(expression, clinical, deg_secretory_early, "early_secretory", go_results, kegg_results)
write.table(results_secretory_early$go, file = "./results_secretory_early_go.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_secretory_early$kegg, file = "./results_secretory_early_kegg.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

results_secretory_mid <- perform_anova(expression, clinical, deg_secretory_mid, "mid_secretory", go_results, kegg_results)
write.table(results_secretory_mid$go, file = "./results_secretory_mid_go.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_secretory_mid$kegg, file = "./results_secretory_mid_kegg.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
