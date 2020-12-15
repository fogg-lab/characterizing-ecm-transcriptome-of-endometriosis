library(tidyverse)
library(HDF5Array)

figo_map_df <- tibble(
    roman_num = c("I", "II", "III", "IV"),
    figo_chr = c('figo_stage_1', 'figo_stage_2', 'figo_stage_3', 'figo_stage_4'),
    figo_num = c(1, 2, 3, 4)
)

load_matrisome_df <- function(matrisome_list_file) {
    matrisome_df <- readr::read_tsv(matrisome_list_file, quote = "")
    colnames(matrisome_df) <- purrr::map(sub(" ", "_", colnames(matrisome_df)), tolower)
    matrisome_df <- select(matrisome_df, gene_symbol, everything()) %>%
        dplyr::filter(division != "Retired")    # Ignore "Retired" matrisome genes
    return(matrisome_df)
}


load_survival_df <- function(survival_data_file, event_code) {
    survival_df <- read_tsv(survival_data_file) %>%
        mutate(vital_status_num = case_when(
            vital_status == "Dead" ~ event_code[["Dead"]],
            vital_status == "Alive" ~ event_code[["Alive"]]
        )) %>%
        dplyr::select(sample_name, vital_status_num, everything(), -vital_status) %>%
        dplyr::rename(vital_status = vital_status_num)
    return(survival_df)
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
        samples[[group]] <- res
    }
    return(bind_rows(samples) %>% dplyr::arrange_at(sample_col))
}


get_hm_dfs <- function(counts_df, coldata_df, hm_sample_meta_df, drop_low_sd = TRUE, eps = 1e-6) {
    hm_sample_counts_df <- counts_df[, c("external_gene_name", hm_sample_meta_df$sample_name)] %>%
        column_to_rownames(var = "external_gene_name")
    hm_sample_coldata_df <- coldata_df %>%
        dplyr::filter(sample_name %in% hm_sample_meta_df$sample_name) %>%
        arrange(match(sample_name, hm_sample_meta_df$sample_name)) %>%
        column_to_rownames(var = "sample_name") # Heatmap needs row names

    if (drop_low_sd) {
        gene_sd_mask <- rowSds(as.matrix(hm_sample_counts_df)) > eps
        hm_sample_counts_df <- hm_sample_counts_df[gene_sd_mask, ]
    }

    hm_sample_ls <- list(
        coldata = hm_sample_coldata_df,
        counts = hm_sample_counts_df
    )
    stopifnot(all(colnames(hm_sample_counts_df) == rownames(hm_sample_coldata_df)))
    return(hm_sample_ls)
}


get_hm_clusters <- function(counts_df, col_cor = "spearman", row_cor = "pearson", col_metric = "complete", row_metric = "complete") {
    col_dist <- as.dist(1 - cor(as.matrix(counts_df), method = col_cor))
    row_dist <- as.dist(1 - cor(t(as.matrix(counts_df)), method = row_cor))

    return(list(
        col = hclust(col_dist, method = col_metric),
        row = hclust(row_dist, method = row_metric)
    ))
}


load_RSE_objects <- function(dir, projects, prefixes) {
    data_ls <- list()
    for (i in seq_len(length(projects))) {
        data_ls[[projects[i]]] <- loadHDF5SummarizedExperiment(dir = dir, prefix = prefixes[i])
    }
    return(data_ls)
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


decode_figo_stage <- function(df, to = "num") {
    if (str_sub(to, 1, 1) == "n") {
        drop_col <- "figo_chr"
        keep_col <- "figo_num"
    }
    else if (str_sub(to, 1, 1) == "c") {
        drop_col <- "figo_num"
        keep_col <- "figo_chr"
    }

    new_df <- df %>%
        dplyr::mutate(
            figo_rn = str_extract(figo_stage, "IV|III|II|I")
        ) %>%
        dplyr::inner_join(figo_map_df, by = c("figo_rn" = "roman_num")) %>%
        dplyr::select(-c("figo_rn", "figo_stage", drop_col)) %>%
        dplyr::rename("figo_stage" = keep_col)
    return(new_df)
}


get_unified_group_samples <- function(counts_df, coldata_df, sample_group) {
    if (sample_group == "GTEX") {
        group_samples <- (coldata_df %>%
            dplyr::filter(data_source == "GTEx"))$sample_name
    }
    else if (sample_group == "TCGA_healthy") {
        group_samples <- (coldata_df %>%
            dplyr::filter(condition == "healthy" & data_source == "TCGA"))$sample_name
    }
    else if (sample_group == "TCGA_tumor") {
        group_samples <- (coldata_df %>%
            dplyr::filter(condition == "tumor"))$sample_name
    }
    
    return(counts_df %>% dplyr::select("geneID", all_of(group_samples)))
}


get_unified_thresh_results <- function(group_df, thresh, group_name) {
    over_thresh_str <- paste0(group_name, "_over_thresh")
    over_thresh_prop_str <- paste0(group_name, "_over_thresh_prop")

    res_df <- group_df %>%
        mutate(over_thresh = rowSums(. [, -1] > thresh)) %>%
        mutate(over_thresh_prop = over_thresh / (ncol(.) - 2)) %>%
        dplyr::rename(!!over_thresh_str := over_thresh) %>%
        dplyr::rename(!!over_thresh_prop_str := over_thresh_prop) %>%
        dplyr::select(matches(c("geneID", over_thresh_str, over_thresh_prop_str)))
    return(res_df)
}


get_unified_thresh_results_for_all <- function(counts_df, coldata_df, group_names, thresh) {
    df_list <- list()
    for (gn in group_names) {
        group_counts_df <- get_unified_group_samples(counts_df, coldata_df, gn)
        thresh_res_df <- get_unified_thresh_results(group_counts_df, thresh, gn)
        df_list[[gn]] <- thresh_res_df
    }
    final_df <- df_list %>%
        purrr::reduce(inner_join, by = "geneID") %>%
        mutate(
            tot_over_thresh = rowSums(select(., paste0(group_names, "_over_thresh"))),
            tot_over_thresh_prop = tot_over_thresh / (ncol(counts_df) - 1)    # subtract 1 because we assume gene names is a column
        )
    return(final_df)
}


transpose_df <- function(df, future_colnames_col, previous_colnames_col) {
    temp_df <- as.data.frame(df)
    rownames(temp_df) <- df[[future_colnames_col]]
    temp_df <- temp_df %>% dplyr::select(-(!!future_colnames_col))
    t(temp_df) %>% as_tibble(rownames = previous_colnames_col)
}


filter_outliers_IQR <- function(df, filter_col, coef) {
    quantiles <- quantile(df[[filter_col]])
    q_1 <- quantiles[2]
    q_3 <- quantiles[4]
    iqr <- q_3 - q_1
    min_thresh <- q_1 - coef * iqr
    max_thresh <- q_3 + coef * iqr

    filtered_df <- df %>%
        dplyr::mutate(outlier_status = (!!as.name(filter_col) < min_thresh) | (max_thresh < !!as.name(filter_col))) %>%
        dplyr::filter(outlier_status == FALSE) %>%
        dplyr::select(-outlier_status)

    return(filtered_df)
}


load_matrisome_norm_counts <- function(dset_path, event_code, keep_clinical = NULL, keep_conditions = c("healthy", "tumor"), drop_unexpressed = TRUE) {
    # Load clinical data
    clinical_df <- load_survival_df(paste0(dset_path, "/survival_data.tsv"), event_code) %>%
        # {if (!is.null(keep_clinical)) dplyr::select(., one_of(keep_clinical)) else .} %>%
        dplyr::select(one_of(keep_clinical)) %>%
        dplyr::filter(rowSums(is.na(.)) == 0)
    coldata_df <- read_tsv(paste0(dset_path, "/coldata.tsv")) %>%
        dplyr::filter(condition %in% keep_conditions)
    
    if (drop_unexpressed) {
        unexpressed_genes <- read_tsv(paste0(dset_path, "/matrisome_counts.tsv")) %>%
            dplyr::select(one_of("geneID", coldata_df$sample_name)) %>%
            dplyr::filter(rowSums(.[, -1]) == 0) %>%
            dplyr::pull(geneID)
    } else {
        unexpressed_genes <- c()
    }

    if(is.null(keep_clinical)) {
        keeper_samples <- coldata_df$sample_name
    } else {
        keeper_samples <- clinical_df$sample_name
    }

    matrisome_norm_counts_df <- read_tsv(paste0(dset_path, "/norm_matrisome_counts.tsv")) %>%
        dplyr::select(one_of(c("geneID", coldata_df$sample_name))) %>%
        dplyr::filter(!(geneID %in% unexpressed_genes)) %>%    # Drop unexpressed genes (if desired)
        transpose_df(future_colnames_col = "geneID", previous_colnames_col = "sample_name") %>%
        dplyr::filter(sample_name %in% keeper_samples) %>%   # Drop samples w/ missing clinical data
        inner_join(coldata_df, by = "sample_name") %>%
        dplyr::select(sample_name, condition, data_source, everything())
    return(matrisome_norm_counts_df)
}
