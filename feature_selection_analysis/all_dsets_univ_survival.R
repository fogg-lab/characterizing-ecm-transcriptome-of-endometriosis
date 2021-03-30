library(tidyverse)
library(survival)
library(survminer)
library(WGCNA)

# Custom package
library(rutils)


# Define constants
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_path <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

event_code <- list("Alive" = 0, "Dead" = 1)
dep_cols <- c("vital_status", "survival_time")


# Functions
get_high_low <- function(df, col_str, cutoff) {
    col = as.name(col_str)
    df %>%
        mutate(high = !!col > cutoff, high_low = ifelse(high == TRUE, "high", "low")) %>%
        select(-high)
}

test_all_genes_km <- function(count_df, cutoff_df, gene_names) {
    n_genes <- length(gene_names)
    pvals <- rep(1, n_genes)
    for (i in seq_len(n_genes)) {
        gene_i <- gene_names[i]
        simp_df <- count_df %>%
            select(sample_name, vital_status, survival_time, !!as.name(gene_i))
        cutoff <- cutoff_df %>%
            filter(geneID == gene_i) %>%
            pull(cutoff)
        # If there's only one group (no good cutoff), use median
        tryCatch({
            simp_df <- get_high_low(simp_df, gene_i, cutoff)
            km_fit <- survfit(Surv(survival_time, vital_status) ~ high_low, type = "kaplan-meier", data = simp_df)
            km_diff <- survdiff(Surv(survival_time, vital_status) ~ high_low, data = simp_df)
        }, error = function(error_condition) {
            simp_df <- get_high_low(simp_df, gene_i, median(count_df[[gene_i]]))
            km_fit <- survfit(Surv(survival_time, vital_status) ~ high_low, type = "kaplan-meier", data = simp_df)
            km_diff <- survdiff(Surv(survival_time, vital_status) ~ high_low, data = simp_df)
        })
        pvals[i] <- pchisq(km_diff$chisq, length(km_diff$n) - 1, lower.tail = FALSE)
    }
    tibble(geneID = gene_names, km_pval = pvals, km_qval = WGCNA::qvalue(km_pval)$qvalues)
}

test_all_genes_cph <- function(count_df, gene_names) {
    n_genes <- length(gene_names)
    pvals <- rep(1, n_genes)
    coeffs <- rep(0, n_genes)
    for (i in seq_len(n_genes)) {
        gene_i <- gene_names[i]
        cph_fit <- coxph(as.formula(paste0("Surv(survival_time, vital_status) ~ ", gene_i)), data = count_df)
        pvals[i] <- summary(cph_fit)$logtest["pvalue"]
        coeffs[i] <- as.numeric(cph_fit$coefficient)
    }
    tibble(geneID = gene_names, cph_pval = pvals, cph_qval = WGCNA::qvalue(cph_pval)$qvalues, coeff = coeffs)
}

for (dset_idx in 1:3) {
    # Load and filter survival data
    survival_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/survival_data.tsv")
    survival_df <- load_survival_df(survival_path, event_code)
    umsmg_demg_list <- read_lines(paste0(dirs$analysis_dir, "/gene_lists/", unified_dsets[dset_idx], "_umsmg_demg_list.txt"))
    cutoff_df <- read_tsv(paste0(dirs$analysis_dir, "/survival/", unified_dsets[dset_idx], "_expression_cutoffs.tsv"))

    filtered_survival_df <- survival_df %>%
        dplyr::select(one_of(c("sample_name", dep_cols))) %>%
        dplyr::filter(rowSums(is.na(.)) == 0)

    # Load normalized matrisome count data
    norm_matrisome_counts_df <- read_tsv(paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/norm_matrisome_counts.tsv"))
    norm_matrisome_counts_t_df <- norm_matrisome_counts_df %>%
        dplyr::select(c("geneID", filtered_survival_df$sample_name)) %>%
        transpose_df("geneID", "sample_name")
    # Combine survival data and normalized count data
    filtered_joined_df <- filtered_survival_df %>%
        inner_join(norm_matrisome_counts_t_df, by = "sample_name") %>%
        select(one_of("sample_name", "vital_status", "survival_time", umsmg_demg_list)) %>%
        # cannot have survival times of 0 for univariate Cox PH analysis
        dplyr::filter(survival_time > 0)

    gene_names <- colnames(filtered_joined_df)[-(1:3)]
    km_df <- test_all_genes_km(filtered_joined_df, cutoff_df, gene_names)
    cph_df <- test_all_genes_cph(filtered_joined_df, gene_names)
    joined_surv_df <- km_df %>%
        inner_join(cph_df, by = "geneID")
    write_tsv(joined_surv_df, paste0(dirs$analysis_dir, "/survival/", unified_dsets[dset_idx], "_univ_survival_results.tsv"))
}
