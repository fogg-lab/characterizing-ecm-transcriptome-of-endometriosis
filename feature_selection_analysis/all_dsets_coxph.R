library(tidyverse)
library(survival)

# Custom package
library(rutils)

# Define constants and load data
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_path <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

event_code <- list("Alive" = 0, "Dead" = 1)
covariate_cols <- c("figo_stage", "age_at_diagnosis", "race", "ethnicity")
dep_cols <- c("vital_status", "survival_time")
figo_map_df <- tibble(
    roman_num = c("I", "II", "III", "IV"),
    figo_code = c('1', '2', '3', '4')
)

cox_null_scores_df <- tibble(score = c("lr_test_pval", "wald_test_pval", "score_test_pval"))

for (dset_idx in 1:3) {
    # Load and filter survival data
    survival_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/survival_data.tsv")
    survival_df <- load_survival_df(survival_path, event_code)
    filtered_survival_df <- survival_df %>%
        dplyr::select(one_of(c("sample_name", dep_cols, covariate_cols))) %>%
        dplyr::filter(rowSums(is.na(.)) == 0) %>%
        dplyr::mutate(
            figo_rn = str_extract(figo_stage, "IV|III|II|I")
        ) %>%
        dplyr::inner_join(figo_map_df, by = c("figo_rn" = "roman_num")) %>%
        dplyr::select(-c(figo_rn, figo_stage)) %>%
        dplyr::rename(figo_stage = figo_code)

    # Load normalized matrisome count data
    norm_matrisome_counts <- read_tsv(paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/norm_matrisome_counts.tsv")) %>%
        column_to_rownames(var = "geneID") %>%
        as.matrix()

    # Match up columns of counts with rows of survival data & only include samples present in survival data
    norm_matrisome_survival_counts <- norm_matrisome_counts[, filtered_survival_df$sample_name]
    all(rownames(t(norm_matrisome_survival_counts)) == filtered_survival_df$sample_name)

    # Combine survival data and normalized count data
    joined_survival_counts_df <- filtered_survival_df %>%
    inner_join(
        as_tibble(t(norm_matrisome_survival_counts), rownames = "sample_name"),
        by = "sample_name"
    )

    # Some genes contain the '-' symbol, which affects formulae
    colnames(joined_survival_counts_df) <- gsub("-", "_", colnames(joined_survival_counts_df))

    # Fit null model
    null_model_formula_chr <- paste0(
        "Surv(survival_time, vital_status) ~ ",
        paste0(covariate_cols, collapse = " + ")
    )
    cox_fit_null <- coxph(
        as.formula(null_model_formula_chr),
        data = joined_survival_counts_df,
        singular.ok = TRUE
    )
    res_cox_fit_null <- summary(cox_fit_null)
    dset <- unified_dsets[dset_idx]
    cox_null_scores_df <- cox_null_scores_df %>%
        add_column(!!dset := c(res_cox_fit_null$logtest[["pvalue"]], res_cox_fit_null$waldtest[["pvalue"]], res_cox_fit_null$sctest[["pvalue"]]))

    # Fit gene models
    genes_of_interest <- colnames(joined_survival_counts_df %>% dplyr::select(-colnames(filtered_survival_df)))
    gene_pvals <- c()
    gene_coeffs <- c()

    for (g in genes_of_interest) {
        gene_model_formula_chr <- paste0(null_model_formula_chr, " + ", g)
        cox_fit_gene <- coxph(
            as.formula(gene_model_formula_chr),
            data = joined_survival_counts_df,
            singular.ok = TRUE
        )
        anova_res <- anova(cox_fit_null, cox_fit_gene, test = "LRT")
        gene_pvals <- c(gene_pvals, anova_res[["P(>|Chi|)"]][2])
        gene_coeffs <- c(gene_coeffs, cox_fit_gene$coefficients[[g]])
    }

    # Re-sub '-' for '_' now that no longer needed for formulae
    cox_regression_df <- tibble("geneID" = gsub("_", "-", genes_of_interest), "gene_pval" = gene_pvals, "gene_coeff" = gene_coeffs)
    sig_cox_regression_df <- cox_regression_df %>%
        dplyr::filter(gene_pval < 0.05)

    # Save results
    write_tsv(cox_regression_df, paste0(dirs$analysis_dir, "/", unified_dsets[dset_idx], "_coxph_results.tsv"))

    # Num. predictive genes
    nrow(sig_cox_regression_df)
    # Prop. matrisome genes which are predictive
    nrow(sig_cox_regression_df) / nrow(norm_matrisome_survival_counts)
    # Genes associated with negative prognosis
    nrow(sig_cox_regression_df %>%
        dplyr::filter(gene_coeff > 0))

    # Genes associated with positive prognosis
    nrow(sig_cox_regression_df %>%
        dplyr::filter(gene_coeff < 0))
}

write_tsv(cox_null_scores_df, paste0(dirs$analysis_dir, "/", "coxph_null_scores.tsv"))
