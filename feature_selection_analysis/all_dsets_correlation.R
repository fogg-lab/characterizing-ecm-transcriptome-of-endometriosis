library(tidyverse)

# Custom package
library(rutils)

# Define constants and load data
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_path <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

event_code <- list("Alive" = 0, "Dead" = 1)
covariate_cols_no_figo <- c("age_at_diagnosis", "race", "ethnicity")
covariate_cols <- c("figo_stage", covariate_cols_no_figo)
dep_cols <- c("vital_status", "survival_time")


for (dset_idx in 1:3) {
    # Load and filter survival data
    survival_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/survival_data.tsv")
    survival_df <- load_survival_df(survival_path, event_code)
    filtered_survival_df <- survival_df %>%
        dplyr::select(sample_name, vital_status, survival_time) %>%
        dplyr::filter(vital_status == event_code$Dead, rowSums(is.na(.)) == 0)

    # Load normalized matrisome count data
    norm_matrisome_counts_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/norm_matrisome_counts.tsv")
    norm_survival_counts_df <- read_tsv(norm_matrisome_counts_path) %>%
        dplyr::select(c("geneID", filtered_survival_df$sample_name))
    norm_survival_counts_t_df <- norm_survival_counts_df %>%
        column_to_rownames(var = "geneID") %>%
        t() %>%
        as_tibble(rownames = "sample_name") %>%
        inner_join(filtered_survival_df, by = "sample_name") %>%
        dplyr::select(sample_name, survival_time, everything(), -vital_status)

    # Perform correlation  test
    cor_test_df <- colwise_cor_test(
        norm_survival_counts_t_df,
        colnames(norm_survival_counts_t_df)[-c(1:2)],
        "survival_time",
        v = "geneID"
    ) %>%
        dplyr::mutate(padj = p.adjust(pval, method = "BH")) %>%
        dplyr::select(geneID, cor, pval, padj, n) %>%
        na.omit() %>%
        dplyr::mutate(qval = WGCNA::qvalue(pval)$qvalues) %>%
        dplyr::select(geneID, cor, pval, padj, qval, n)

    # Save results
    write_tsv(cor_test_df, paste0(dirs$analysis_dir, "/", unified_dsets[dset_idx], "_cor_results.tsv"))
}