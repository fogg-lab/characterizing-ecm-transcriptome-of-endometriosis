library(tidyverse)
library(survival)
library(WGCNA)

# Custom package
library(rutils)

# Define constants
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_path <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

event_code <- list("Alive" = 0, "Dead" = 1)
covariate_cols <- c("figo_stage", "age_at_diagnosis", "race", "ethnicity")
dep_cols <- c("vital_status", "survival_time")


cox_null_scores_df <- tibble(score = c("lr_test_pval", "wald_test_pval", "score_test_pval"))

for (dset_idx in 1:3) {
    # Load and filter survival data
    survival_path <- paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/survival_data.tsv")
    survival_df <- load_survival_df(survival_path, event_code)

    filtered_survival_df <- survival_df %>%
        decode_figo_stage(to = "c") %>%
        dplyr::select(one_of(c("sample_name", dep_cols, covariate_cols))) %>%
        dplyr::filter(rowSums(is.na(.)) == 0)

    # Load normalized matrisome count data
    norm_matrisome_counts_df <- read_tsv(paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/norm_matrisome_counts.tsv"))
    norm_matrisome_counts_t_df <- norm_matrisome_counts_df %>%
        dplyr::select(c("geneID", filtered_survival_df$sample_name)) %>%
        transpose_df("geneID", "sample_name")
        # column_to_rownames(var = "geneID") %>%
        # as.matrix()

    # # Match up columns of counts with rows of survival data & only include samples present in survival data
    # norm_matrisome_survival_counts <- norm_matrisome_counts[, filtered_survival_df$sample_name]
    # all(rownames(t(norm_matrisome_survival_counts)) == filtered_survival_df$sample_name)

    # Combine survival data and normalized count data
    joined_survival_counts_df <- filtered_survival_df %>%
    inner_join(norm_matrisome_counts_t_df, by = "sample_name")

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
    cox_regression_df <- tibble("geneID" = gsub("_", "-", genes_of_interest), "gene_pval" = gene_pvals, "gene_coeff" = gene_coeffs) %>%
        dplyr::mutate(gene_qval = WGCNA::qvalue(gene_pval)$qvalues) %>%
        dplyr::mutate(gene_padj = p.adjust(gene_pval, method = "BH"))

    # Save results
    write_tsv(cox_regression_df, paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_coxph_results.tsv"))

}

write_tsv(cox_null_scores_df %>% rutils::transpose_df("score", "dataset"), paste0(dirs$analysis_dir, "/meta/", "coxph_null_scores.tsv"))
