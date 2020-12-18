library(tidyverse)
library(DESeq2)
library(BiocParallel)

# Custom package
library(rutils)


# Set number of cores to use
n_cores <- detectCores() - 2
BiocParallel::register(MulticoreParam(n_cores))


# Define constants and load data
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
dset_paths <- unlist(map(dsets, function(d) paste0(dirs$data_dir, "/", d)))
matrisome_list <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

padj_thresh <- 0.05
expr_thresh <- 0


# Helper functions
run_DESeq_and_get_results <- function(dds) {
    dds_seq <- DESeq(dds, parallel = TRUE)
    res <- results(
        dds_seq,
        contrast = c("condition", "tumor", "healthy"),
        pAdjustMethod = "BH",
        alpha = padj_thresh,
        parallel = TRUE
    )
    return(res)
}

for (dset_idx in 1:3) {
    # Read in data
    counts_df <- read_tsv(paste0(dset_paths[dset_idx], "/counts.tsv")) %>%
        mutate_if(is.numeric, round, 0) %>%
        dplyr::select(-Entrez_Gene_Id) %>%
        dplyr::rename(geneID = Hugo_Symbol)
    coldata_df <- read_tsv(paste0(dset_paths[dset_idx], "/coldata.tsv"))
    all(coldata_df$sample_name == colnames(counts_df[, -1]))

    matrisome_genes_df <- rutils::load_matrisome_df(matrisome_list) %>%
        dplyr::select(gene_symbol)

    # Population breakdown
    sum(coldata_df$condition == "healthy")
    sum(coldata_df$condition == "tumor")

    # Pre-filter genes
    all_counts_res_df <- rutils::get_unified_thresh_results_for_all(
        counts_df,
        coldata_df,
        c("GTEX", "TCGA_healthy", "TCGA_tumor"),
        thresh = 0
    )

    sufficiently_expr_genes_df <- all_counts_res_df %>%
        dplyr::filter(tot_over_thresh_prop > 1/3)

    nrow(sufficiently_expr_genes_df)
    # Proportion of genes which will be kept
    nrow(sufficiently_expr_genes_df) / nrow(counts_df)

    filtered_counts_df <- counts_df %>%
        dplyr::filter(geneID %in% sufficiently_expr_genes_df$geneID)

    # Run DGE analysis
    dds <- DESeqDataSetFromMatrix(
        countData = filtered_counts_df %>% column_to_rownames(var = "geneID"),
        colData = coldata_df %>% column_to_rownames(var = "sample_name"),
        design = ~ condition
    )
    dge_res <- run_DESeq_and_get_results(dds)
    dge_res_df <- as_tibble(dge_res, rownames = "geneID") %>%
        dplyr::mutate(qval = WGCNA::qvalue(pvalue)$qvalues)

    # Save results
    write_tsv(dge_res_df, paste0(dirs$analysis_dir, "/", dsets[dset_idx], "_DESeq_results.tsv"))
}
