library(readr)
library(dplyr)
library(sva)

args = commandArgs(trailingOnly = TRUE)

# Paths
countspath = args[1]
coldatapath = args[2]

# Read in the sample data and counts
counts <- read_tsv(countspath)
metadata <- read_tsv(coldatapath)

# Prep batch correction
rma_expr <- as.matrix(counts[-1])
rownames(rma_expr) <- counts$hgnc_symbol
model_m <- model.matrix(~ condition, data = metadata)
batch <- metadata$batch

# Get symbol column for adding to output later
symbol <- counts$hgnc_symbol
symbols <- data.frame(symbol)

# Batch correction
bc_rma_expr <- ComBat(rma_expr, batch = batch, mod = model_m, ref.batch = 1)

# Add symbol column and write the dataframe to file
bc_rma_expr <- data.frame(bc_rma_expr)
result <- bind_cols(symbols, bc_rma_expr)
write_tsv(result, file.path(countspath))
