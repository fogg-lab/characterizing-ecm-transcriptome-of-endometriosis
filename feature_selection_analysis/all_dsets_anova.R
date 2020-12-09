library(tidyverse)

# Custom package
library(rutils)


# Define constants
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_path <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

event_code <- list("Alive" = 0, "Dead" = 1)
covariate_cols_no_figo <- c("age_at_diagnosis", "race", "ethnicity")
covariate_cols <- c("figo_stage", covariate_cols_no_figo)
dep_cols <- c("vital_status", "survival_time")

for (dset_idx in 1:3) {
    survival_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/survival_data.tsv")
    survival_df <- load_survival_df(survival_path, event_code)

    # Load and filter survival data
    filtered_survival_df <- survival_df %>%
        decode_figo_stage(to = "c") %>%
        dplyr::select(sample_name, figo_stage, race, ethnicity, age_at_diagnosis) %>% # make sure using same samples as classification models
        dplyr::filter(rowSums(is.na(.)) == 0) %>%
        dplyr::select(sample_name, figo_stage)

    # Load normalized matrisome count data
    norm_matrisome_counts_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/norm_matrisome_counts.tsv")
    norm_survival_counts_t_df <- read_tsv(norm_matrisome_counts_path) %>%
        dplyr::select(c("geneID", filtered_survival_df$sample_name)) %>%
        transpose_df("geneID", "sample_name")

    # Combine survival data and normalized count data
    joined_df <- filtered_survival_df %>%
        inner_join(norm_survival_counts_t_df, by = "sample_name")

    # Some genes contain the '-' symbol, which affects formulae
    colnames(joined_df) <- gsub("-", "_", colnames(joined_df))

    # Perform Welch ANOVA for each gene and apply BH adjustment to p-values
    gene_names <- colnames(joined_df[-c(1:2)])
    waov_df <- tibble(geneID = c(), pval = c())

    waov_df <- colwise_anova(joined_df, "figo_stage", gene_names, "geneID", adjust_method = "BH")

    # Re-sub '-' symbol
    waov_df$geneID <- gsub("_", "-", waov_df$geneID)

    write_tsv(waov_df, paste0(dirs$analysis_dir, "/", unified_dsets[dset_idx], "_welch_anova_results.tsv"))
}
