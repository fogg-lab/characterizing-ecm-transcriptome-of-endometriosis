library(tidyverse)
library(DESeq2)
library(BiocParallel)


data_dir <- "../../../../../mnt/d/unified_TCGA_GTEx"
# dsets <- c("unified_cervical_data", "unified_uterine_data")
dsets <- c("unified_cervical_data")
dset_paths <- unlist(map(dsets, function(d) paste0(data_dir, "/", d)))
padj_thresh <- 0.05


set_bioc_parallel <- function() {
    n_cores <- parallel::detectCores()
    BiocParallel::register(BiocParallel::MulticoreParam(n_cores))
}


# Filter out rows (genes) which are not expressed in any sample
filter_unexpressed <- function(d) {
    not_expr_mask <- rowSums(DESeq2::counts(d)) == 0
    return(d[!not_expr_mask, ])
}


set_bioc_parallel()

for (i in seq_len(length(dset_paths))) {
    # Load data
    # Need to round the counts so that they are integers. Not optimal, but suggested method
    # See: https://support.bioconductor.org/p/94003/#94028
    counts <- read_tsv(paste0(dset_paths[i], "/counts.tsv")) %>%
        select(-"Entrez_Gene_Id") %>%
        mutate_if(is.numeric, round, 0) %>%
        column_to_rownames(var = "Hugo_Symbol")
    coldata <- read_tsv(paste0(dset_paths[i], "/coldata.tsv")) %>%
        column_to_rownames(var = "sample_name")

    # coldata and counts must line up for DESeq2 to work
    # See: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input
    stopifnot(all(rownames(coldata) == colnames(counts)))

    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~ condition
    )
    dds_filtered <- filter_unexpressed(dds)
    dds_seq <- DESeq(dds_filtered, parallel = TRUE)

    res <- results(
        dds_seq,
        contrast = c("condition", "tumor", "healthy"),
        pAdjustMethod = "BH",
        parallel = TRUE
    )
    res_df <- as_tibble(res, rownames = "geneID") %>%
        drop_na() %>%
        filter(padj < padj_thresh)

    DEG_dest <- paste0(
        data_dir,
        "/analysis/",
        paste0(
            dsets[i],
            "_DGE_padj",
            as.character(padj_thresh),
            ".tsv"
        )
    )
    write_tsv(res_df, DEG_dest)
}
