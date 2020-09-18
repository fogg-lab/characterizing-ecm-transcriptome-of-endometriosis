
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


count_pca_results <- function() {
    return()
}
