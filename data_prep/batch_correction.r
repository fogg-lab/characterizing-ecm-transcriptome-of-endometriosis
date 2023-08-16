library(readr)
library(dplyr)
library(sva)

args = commandArgs(trailingOnly = TRUE)

# Paths
expression_file = args[1]
clinical_file = args[2]

# Read in the expression matrix and clinical data table
expression_df <- read_tsv(expression_file)
clinical_df <- read_tsv(clinical_file)

# Prep batch correction
expression_matrix <- as.matrix(expression_df[-1])
rownames(expression_matrix) <- expression_df$hgnc_symbol
model_m <- model.matrix(~ condition, data = clinical_df)
batch <- clinical_df$batch

# Get symbol column for adding to output later
symbol <- expression_df$hgnc_symbol
symbols <- data.frame(symbol)

# Batch correction
bc_expression_matrix <- ComBat(expression_matrix, batch = batch, mod = model_m, ref.batch = 1)

# Add symbol column and write the dataframe to file
bc_expression_matrix <- data.frame(bc_expression_matrix)
result <- bind_cols(symbols, bc_expression_matrix)
write_tsv(result, file.path(expression_file))
