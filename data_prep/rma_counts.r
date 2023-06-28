library(affy)
library(Biobase)

args <- commandArgs(trailingOnly=TRUE)

data_dir <- args[1]
output_dir <- args[2]

series_dirs <- list.dirs(data_dir, recursive=FALSE)

saveCounts <- function(expr_df, series) {
    output_file = paste0(output_dir, "/", series, "_counts_unmapped.tsv")
    write.table(expr_df, file=output_file, sep="\t", quote=F, row.names=FALSE)
}

writeLines("Preprocessing data in the following directories:")
writeLines(series_dirs)
writeLines("------------------------------------------------\n")

for (series_dir in series_dirs) {
    series <- basename(series_dir)

    # Create output directory if it doesn't exist
    dir.create(output_dir, recursive=FALSE, showWarnings=F)

    # Read in CEL files
    raw_data <- ReadAffy(celfile.path = series_dir)

    # Perform RMA normalization
    rma_data <- rma(raw_data)

    # Extract expression matrix from the ExpressionSet
    expr_matrix <- exprs(rma_data)

    # Convert to data frame
    expr_df <- as.data.frame(expr_matrix)

    # Add probe column
    expr_df <- data.frame(probe = rownames(expr_df), expr_df, row.names = NULL)

    # Write expression matrix to file
    saveCounts(expr_df, series)
}
