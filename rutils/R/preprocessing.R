library(tidyverse)

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