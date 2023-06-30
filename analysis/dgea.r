library(readr)
library(dplyr)
library(tibble)
library(limma)
library(yaml)
library(ggrepel)
library(jsonlite)

# Parse command line arguments for file paths
args = commandArgs(trailingOnly = TRUE)
counts_filepath <- args[1]
coldata_filepath <- args[2]
config_filepath <- args[3]
filter_filepath <- args[4]
output_dir <- args[5]

# Read data from files
counts_df <- read_tsv(counts_filepath, col_types=cols())
coldata_df <- read_tsv(coldata_filepath, col_types=cols())

# Read configuration settings from YAML file
config_yml <- read_yaml(config_filepath)

# Extract config parameters
min_expr <- config_yml$min_expr
min_prop <- config_yml$min_prop
padj_thresh <- config_yml$padj_thresh
adj_method <- config_yml$adj_method
contrast_level <- config_yml$contrast_level
reference_level <- config_yml$reference_level
use_qual_weights <- config_yml$use_qual_weights

# Convert column names to lower case
names(coldata_df) <- tolower(names(coldata_df))

# Apply factor levels to the 'condition' column
coldata_df <- coldata_df %>%
    mutate(condition = factor(condition, levels = c(reference_level, contrast_level)))

# Filter count data based on the minimum expression and proportion criteria
filt_counts_df <- counts_df %>%
    filter(rowSums(.[-1] > min_expr) / (ncol(.) - 1) >= min_prop)

# Convert filtered count data to a matrix
filt_expr <- filt_counts_df %>%
    column_to_rownames("symbol") %>%
    as.matrix()

# Create a design matrix for differential expression analysis
design <- model.matrix(~ condition, data = coldata_df)
rownames(design) <- coldata_df$sample_name

# Check if column names of the expression data match the row names of the design matrix
all(colnames(filt_expr) == rownames(design))

# Calculate weights if the flag for using quality weights is set
if (use_qual_weights) {
    qual_weights <- arrayWeights(filt_expr, design = design)
} else {
    qual_weights <- NULL
}

# Fit the linear model
lm_fit <- lmFit(filt_expr, design = design, weights = qual_weights)

# Apply empirical Bayes smoothing to the standard errors
bayes_fit <- eBayes(lm_fit)

# Get the column names of the coefficients in the fitted model
bayes_fit$coefficients %>% colnames()

# Extract DEGs and rename columns
fit_de_res_df <- topTable(bayes_fit, coef = paste0("condition", contrast_level), number = nrow(filt_counts_df), adjust.method = adj_method, p.value = padj_thresh) %>%
    rename(lfc = logFC, ave_expr = AveExpr, pval = P.Value, padj = adj.P.Val) %>%
    as_tibble(rownames = "symbol")

# Rename columns
colnames(fit_de_res_df) <- c("symbol", "l2fc", "base_avg", "test_stat", "pval", "padj", "B")

# Write the result to a tsv file
write_tsv(fit_de_res_df, file.path(output_dir, "output.tsv"))

# Filter the results if a gene filter list is provided
if (file.info(filter_filepath)$size != 0) {
    filter_json <- fromJSON(filter_filepath)
    filter_list <- filter_json$symbols
    filtered_df <- fit_de_res_df[fit_de_res_df$symbol %in% filter_list,]
    write_tsv(filtered_df, file.path(output_dir, "filter_output.tsv"))
}
