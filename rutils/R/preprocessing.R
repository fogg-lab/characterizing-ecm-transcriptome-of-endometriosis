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


load_survival_df <- function(survival_data_path, event_code) {
    survival_df <- read_tsv(survival_data_path) %>%
        mutate(vital_status_num = case_when(
            vital_status == "Dead" ~ event_code[["Dead"]],
            vital_status == "Alive" ~ event_code[["Alive"]]
        )) %>%
        dplyr::select(sample_name, vital_status_num, everything(), -vital_status) %>%
        dplyr::rename(vital_status = vital_status_num)
    return(survival_df)
}


to_one_hot <- function(df, col) {
    one_hot <- model.matrix(
        as.formula(paste0("~ ", col, " - 1")),    # We do not want the intercept
        model.frame(~ ., df[col], na.action = na.pass)
    )
    # Don't want white space
    colnames(one_hot) <- gsub(" ", "_", colnames(one_hot))
    # model.matrix() will prepend original column name to each one-hot column
    colnames(one_hot) <- gsub(col, paste0(col, "_"), colnames(one_hot))
    return(tibble::as_tibble(one_hot))
}
