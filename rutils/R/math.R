library(tidyverse)
library(Biobase)
library(WGCNA)


dist_L1 <- function(x, y) {
    return(sum(abs(x - y)))
}


get_z_score <- function(x) {
    return((x - mean(x)) / sd(x))
}


find_n_closest <- function(counts, centroid, n, sample_colname) {
    res <- apply(counts, 2, function(x) { dist_L1(x, centroid) })
    return(
        res %>%
            tibble::as_tibble(rownames = sample_colname) %>%
            top_n(-n, wt = value)
    )
}

get_group_centroids <- function(counts, coldata, groups, group_col, sample_col) {
    centroids <- list()
    for (group in groups) {
        group_mask <- coldata[[group_col]] == group
        group_samples <- coldata[[sample_col]][group_mask]
        # Use medians as centroids since we know there are many outliers
        centroid <- rowMedians(as.matrix(counts[, group_samples]))
        centroids[[group]] <- centroid
    }
    return(tibble::as_tibble(centroids))
}


colwise_cor_test <- function(df, indep_vars, dep_var, var = "var", method = "pearson") {
    cors <- c()
    pvals <- c()
    ns <- c()
    for (iv in indep_vars) {
        ct_res <- cor.test(df[[iv]], df[[dep_var]], method = method)
        cors <- c(cors, ct_res$estimate[[1]])
        pvals <- c(pvals, ct_res$p.value[[1]])
        ns <- c(ns, length(df[[iv]]))
    }
    return(
        tibble(!!var := indep_vars, "cor" = cors, "pval" = pvals, "n" = ns)
    )
}


# Source: Page-Gould, Elizabeth. (2015). Re: What is the formula to calculate the critical value of correlation?. Retrieved from: https://www.researchgate.net/post/What_is_the_formula_to_calculate_the_critical_value_of_correlation/558f05c76307d9c3488b4601/citation/download. 
critical_r <- function(n, alpha = .05) {
  df <- n - 2
  critical_t <- qt(alpha / 2, df, lower.tail = FALSE)
  critical_r <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  return(critical_r)
}


colwise_anova <- function(df, group_var, cols, colnames_col, adjust_method = "BH") {
    pvals <- rep(0, length(cols))
    for (i in seq_len(length(cols))) {
        formula_str <- paste0(cols[i], " ~ ", group_var)
        aov_res <- oneway.test(as.formula(formula_str), data = df)
        pvals[i] <- aov_res$p.value
    }
    # aov_df <- tibble(!!as.name(colnames_col) := cols, "pval" = pvals) %>%
    aov_df <- tibble(!!as.name(colnames_col) := cols) %>%
        dplyr::mutate(pval = pvals) %>%
        dplyr::mutate(padj = p.adjust(pval, method = adjust_method))
    return(aov_df)
}


colwise_t_test <- function(df, group_var, cols, colnames_col, adjust_method = "BH") {
    pvals <- rep(0, length(cols))
    for (i in seq_len(length(cols))) {
        formula_str <- paste0(cols[i], " ~ ", group_var)
        welch_t_res <- t.test(as.formula(formula_str), data = df, )
        pvals[i] <- welch_t_res$p.value
    }
    welch_t_df <- tibble(!!as.name(colnames_col) := cols) %>%
        dplyr::mutate(pval = pvals) %>%
        dplyr::mutate(padj = p.adjust(pval, method = adjust_method))
    return(welch_t_df)
}


get_most_conn_genes <- function(data_expr, module_colors, soft_power, type = "unsigned", conn_vs_hub_thresh = 0.5) {
    hub_ls <- list()
    for (mc in unique(module_colors)) {
        mini_adj <- WGCNA::adjacency(data_expr[, module_colors == mc], power = soft_power, type = type)
        mc_hubs_df <- as_tibble(mini_adj, rownames = "geneID") %>%
            dplyr::select(geneID, everything()) %>%
            dplyr::mutate(conn = rowSums(as.matrix(.[-c(1)]))) %>%
            dplyr::select(geneID, conn) %>%
            dplyr::arrange(desc(conn)) %>%
            # Get connectivity of each gene (as a proportion of hub gene's connectivity)
            dplyr::mutate(conn_vs_hub = 1 - abs((conn - first(conn)) / first(conn))) %>%
            dplyr::filter(conn_vs_hub > conn_vs_hub_thresh)
        hub_ls[[mc]] <- mc_hubs_df
    }
    return(hub_ls)
}
