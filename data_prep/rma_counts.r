library(affy)
library(Biobase)
library(GEOquery)

args <- commandArgs(trailingOnly=TRUE)

data_dir <- args[1]
output_dir <- args[2]

series_dirs <- list.dirs(data_dir, recursive=FALSE)

saveCounts <- function(expr_df, series) {
    output_file = paste0(output_dir, "/", series, "_counts_unmapped.tsv")
    write.table(counts, file=output_file, sep="\t", quote=F, row.names=FALSE)
}

saveColdata <- function(coldata, series) {
    # Keep only "sample_name" and "condition" columns
    coldata <- coldata[, c("sample_name", "condition")]

    output_file = paste0(output_dir, "/", series, "_coldata.tsv")
    write.table(coldata, file=output_file, sep="\t", quote=F, row.names = FALSE)
}

processColdata <- function(coldata) {
    condition <- character(nrow(coldata))

    for(i in 1:nrow(coldata)) {
        sample_row <- tolower(paste(coldata[i, ], collapse=" "))
        if(grepl("non-endometriosis|normal patient|disease status: healthy", sample_row)) {
            condition[i] <- "healthy"
        } else if(grepl("endometriosis severity|patient with\\s*endometriosis", sample_row)) {
            condition[i] <- "endometriosis"
        } else {
            condition[i] <- "healthy"
        }
    }

    return(data.frame(sample_name = row.names(coldata), condition = condition))
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

    # Prepare coldata
    out = tryCatch({
        # Get GEO series and corresponding platform
        gset <- getGEO(series, GSEMatrix=TRUE, getGPL=TRUE)
        if (length(gset) > 1) {
            idx <- grep(pattern = series, x = attr(gset, "names"))
        } else idx <- 1
        gset <- gset[[idx]]

        # Get and process coldata
        coldata <- gset@phenoData@data
        coldata <- processColdata(coldata)

        if (!is.null(coldata)) {
            # Write coldata to file
            saveColdata(coldata, series)
        }
    }, error = function(e) {
        print(paste('error:', e))
    })
}
