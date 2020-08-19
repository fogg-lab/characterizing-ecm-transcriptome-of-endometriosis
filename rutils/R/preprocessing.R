library(tidyverse)
library(HDF5Array)

load_matrisome_df <- function(matrisome_list_file) {
    matrisome_df <- readr::read_tsv(matrisome_list_file, quote = "")
    colnames(matrisome_df) <- purrr::map(sub(" ", "_", colnames(matrisome_df)), tolower)
    matrisome_df <- select(matrisome_df, gene_symbol, everything()) %>%
        dplyr::filter(division != "Retired")    # Ignore "Retired" matrisome genes
    return(matrisome_df)
}


load_and_combine_count_matrix_data <- function(count_files, coldata_files, count_join_symbols) {
    count_df_ls <- purrr::map(.x = count_files, .f = readr::read_tsv)
    coldata_df_ls <- purrr::map(.x = coldata_files, .f = readr::read_tsv)

    coldata_df <- dplyr::bind_rows(coldata_df_ls)
    count_df <- purrr::reduce(count_df_ls, dplyr::inner_join, by = count_join_symbols)
    return(list(counts_df = count_df, coldata_df = coldata_df))
}


balanced_group_sample <- function(counts, coldata, centroids, groups, n, group_col, sample_col) {
    samples <- list()
    for (group in groups) {
        group_mask <- coldata[[group_col]] == group
        group_counts <- counts[, coldata[[sample_col]][group_mask]]
        centroid <- centroids[[group]]
        res <- find_n_closest(group_counts, centroid, n, sample_col) %>%
            mutate("group" = group)
        print(dim(res))
        samples[[group]] <- res
    }
    return(bind_rows(samples) %>% dplyr::arrange_at(sample_col))
}


load_RSE_objects <- function(dir, projects, prefixes) {
    data_ls <- list()
    for (i in seq_len(length(projects))) {
        data_ls[[projects[i]]] <- loadHDF5SummarizedExperiment(dir = dir, prefix = prefixes[i])
    }
    return(data_ls)
}
