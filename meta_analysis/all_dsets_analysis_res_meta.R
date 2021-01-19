library(tidyverse)

# Custom package
library(rutils)


# Helper functions
get_consensus_col <- function(df, n, thresh = 0) {
    consensus_df <- df %>%
        dplyr::select(geneID, contains("mean")) %>%
        dplyr::mutate(consensus = rowSums(.[-1] > thresh) == n)
    return(df %>% dplyr::mutate(consensus = consensus_df$consensus))
}


# Define constants
dirs <- rutils::get_dev_directories(dev_paths_file = "../dev_paths.txt")
projects <- c("TCGA-CESC", "TCGA-UCS", "TCGA-UCEC", "TCGA-OV")
unified_dsets <- c("unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data")
matrisome_list <- paste0(dirs$data_dir, "/matrisome/matrisome_hs_masterlist.tsv")

p_thresh = 0.05
lfc_thresh = log2(2)
coxph_coeff_thresh = 0.0
mi_thresh = 0.0
consensus_thresh = 0.0
consensus_n = 5

matrisome_df <- rutils::load_matrisome_df(matrisome_list) %>%
    dplyr::select(gene_symbol, division, category)

#region Init meta DFs
deg_meta_df <- tibble(
    "item" = c("n_degs", "deg_prop", "n_degs_up", "n_degs_down", "n_matrisome_degs", "matrisome_deg_prop", "n_matrisome_degs_up", "n_matrisome_degs_down")
)
coxph_meta_df <- tibble(
    "item" = c("coxph_null_sig", "n_coxph_sig", "prop_coxph_sig", "n_coxph_sig_protective", "n_coxph_sig_harmful")
)
cor_meta_df <- tibble(
    "item" = c("n_cor", "n_cor_down", "n_cor_up")
)
mi_survival_meta_df <- tibble(
    "item" = c("n_mi")
)
mse_gbr_meta_df <- tibble(
    "item" = c("mse_gbr_avg", "mse_gbr_imp", "mse_gbr_pct_imp", "n_mse_gbr_consensus_genes", "mse_baseline")
)
mse_rfr_meta_df <- tibble(
    "item" = c("mse_rfr_avg", "mse_rfr_imp", "mse_rfr_pct_imp", "n_mse_rfr_consensus_genes", "mse_baseline")
)
anova_meta_df <- tibble(
    "item" = c("n_sig")
)
mi_figo_meta_df <- tibble(
    "item" = c("n_mi")
)
f1_gbc_meta_df <- tibble(
    "item" = c("f1_gbc_avg", "f1_gbc_imp", "f1_gbc_pct_imp", "n_f1_gbc_consensus_genes", "baseline")
)
f1_rfc_meta_df <- tibble(
    "item" = c("f1_rfc_avg", "f1_rfc_imp", "f1_rfc_pct_imp", "n_f1_rfc_consensus_genes", "baseline")
)
f1_l1_lr_meta_df <- tibble(
    "item" = c("f1_l1_lr_avg", "f1_l1_lr_imp", "f1_l1_lr_pct_imp", "n_f1_l1_lr_consensus_genes", "baseline")
)
#endregion

for (dset_idx in 1:3) {
    norm_counts_df <- read_tsv(paste0(dirs$data_dir, "/", unified_dsets[dset_idx], "/", "norm_counts.tsv"))

    #region DEG
    DESeq_results_df <- read_tsv(paste0(dirs$analysis_dir, "/deg/", unified_dsets[dset_idx], "_DESeq_results.tsv"))
    filtered_DESeq_results_df <- DESeq_results_df %>%
        dplyr::filter(abs(log2FoldChange) > lfc_thresh, qval < p_thresh)
    filtered_matrisome_DESeq_results_df <- filtered_DESeq_results_df %>%
        dplyr::filter(geneID %in% matrisome_df$gene_symbol)

    n_degs <- nrow(filtered_DESeq_results_df)
    deg_prop <- nrow(filtered_DESeq_results_df) / nrow(norm_counts_df)
    n_degs_up <- nrow(filtered_DESeq_results_df %>% dplyr::filter(log2FoldChange > 0))
    n_degs_down <- nrow(filtered_DESeq_results_df %>% dplyr::filter(log2FoldChange < 0))
    n_matrisome_degs <- nrow(filtered_matrisome_DESeq_results_df)
    matrisome_deg_prop <- nrow(filtered_matrisome_DESeq_results_df) / nrow(matrisome_df)
    n_matrisome_degs_up <- nrow(filtered_matrisome_DESeq_results_df %>% dplyr::filter(log2FoldChange > 0))
    n_matrisome_degs_down <- nrow(filtered_matrisome_DESeq_results_df %>% dplyr::filter(log2FoldChange < 0))

    dset <- unified_dsets[dset_idx]
    deg_meta_df <- deg_meta_df %>%
        add_column(!!dset := c(n_degs, deg_prop, n_degs_up, n_degs_down, n_matrisome_degs, matrisome_deg_prop, n_matrisome_degs_up, n_matrisome_degs_down))
    #endregion

    #region Cox PH
    coxph_null_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/meta/", "coxph_null_scores.tsv"))
    coxph_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_coxph_results.tsv"))
    filtered_coxph_results_df <- coxph_results_df %>%
        dplyr::filter(gene_qval < p_thresh)
    
    coxph_null_sig <- (coxph_null_scores_df %>% dplyr::filter(dataset == unified_dsets[dset_idx]))$lr_test_pval < p_thresh
    n_coxph_sig <- nrow(filtered_coxph_results_df)
    prop_coxph_sig <- n_coxph_sig / nrow(matrisome_df)
    n_coxph_sig_protective <- nrow(filtered_coxph_results_df %>% dplyr::filter(gene_coeff < 0))
    n_coxph_sig_harmful <- nrow(filtered_coxph_results_df %>% dplyr::filter(gene_coeff > 0))

    dset <- unified_dsets[dset_idx]
    coxph_meta_df <- coxph_meta_df %>%
        add_column(!!dset := c(coxph_null_sig, n_coxph_sig, prop_coxph_sig, n_coxph_sig_protective, n_coxph_sig_harmful))
    #endregion

    #region Regression
    ## Correlation
    cor_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_cor_results.tsv"))
    filtered_cor_results_df <- cor_results_df %>% dplyr::filter(qval < p_thresh)

    n_cor <- nrow(filtered_cor_results_df)
    n_cor_down <- nrow(filtered_cor_results_df %>% dplyr::filter(cor < 0))
    n_cor_up <- nrow(filtered_cor_results_df %>% dplyr::filter(cor > 0))

    dset <- unified_dsets[dset_idx]
    cor_meta_df <- cor_meta_df %>%
        add_column(!!dset := c(n_cor, n_cor_down, n_cor_up))

    ## MI (survival)
    mi_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_MI_survival_results.tsv"))
    filtered_mi_results_df <- mi_results_df %>% dplyr::filter(MI_est_median > 0)

    n_mi <- nrow(filtered_mi_results_df)

    dset <- unified_dsets[dset_idx]
    mi_survival_meta_df <- mi_survival_meta_df %>%
        add_column(!!dset := c(n_mi))

    ## Baselines
    reg_baselines_df <- read_tsv(paste0(dirs$analysis_dir, "/meta/", "reg_baselines.tsv"))
    mae_baseline = (reg_baselines_df %>% filter(dataset == unified_dsets[dset_idx]))$L1
    mse_baseline = (reg_baselines_df %>% filter(dataset == unified_dsets[dset_idx]))$L2
    ev_baseline = (reg_baselines_df %>% filter(dataset == unified_dsets[dset_idx]))$explained_variance

    ## GBR (MSE)
    mse_gbr_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_mse_gbr_ref_scores.tsv"))
    mse_gbr_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_mse_gbr_results.tsv"))
    mse_gbr_results_df <- get_consensus_col(mse_gbr_results_df, n = consensus_n, thresh = consensus_thresh)

    # mult. by -1 to convert neg MSE -> MSE
    mse_gbr_avg <- mean(-mse_gbr_scores_df$ref_score)
    mse_gbr_imp <- mse_gbr_avg < mse_baseline
    # mult. by -1 since pct. "improvement" is pct. reduction in this case
    mse_gbr_pct_imp <- - (mse_gbr_avg - mse_baseline) / mse_baseline * 100
    n_mse_gbr_consensus_genes <- nrow(mse_gbr_results_df %>% dplyr::filter(consensus == TRUE))

    dset <- unified_dsets[dset_idx]
    mse_gbr_meta_df <- mse_gbr_meta_df %>%
        add_column(!!dset := c(mse_gbr_avg, mse_gbr_imp, mse_gbr_pct_imp, n_mse_gbr_consensus_genes, mse_baseline))

    ## RFR (MSE)
    mse_rfr_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_mse_rfr_ref_scores.tsv"))
    mse_rfr_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_mse_rfr_results.tsv"))
    mse_rfr_results_df <- get_consensus_col(mse_rfr_results_df, n = consensus_n, thresh = consensus_thresh)

    # mult. by -1 to convert neg MSE -> MSE
    mse_rfr_avg <- mean(-mse_rfr_scores_df$ref_score)
    mse_rfr_imp <- mse_rfr_avg < mse_baseline
    # mult. by -1 since pct. "improvement" is pct. reduction in this case
    mse_rfr_pct_imp <- - (mse_rfr_avg - mse_baseline) / mse_baseline * 100
    n_mse_rfr_consensus_genes <- nrow(mse_rfr_results_df %>% dplyr::filter(consensus == TRUE))

    dset <- unified_dsets[dset_idx]
    mse_rfr_meta_df <- mse_rfr_meta_df %>%
        add_column(!!dset := c(mse_rfr_avg, mse_rfr_imp, mse_rfr_pct_imp, n_mse_rfr_consensus_genes, mse_baseline))
    #endregion

    #region Classification
    ## ANOVA (FIGO)
    anova_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_welch_anova_results.tsv"))
    filtered_anova_results_df <- anova_results_df %>% dplyr::filter(qval < p_thresh)

    n_anov_sig <- nrow(filtered_anova_results_df)

    dset <- unified_dsets[dset_idx]
    anova_meta_df <- anova_meta_df %>%
        add_column(!!dset := c(n_anov_sig))

    ## MI (FIGO)
    mi_figo_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_MI_figo_results.tsv"))
    filtered_mi_figo_results_df <- mi_figo_results_df %>% dplyr::filter(MI_est_median > 0)

    n_mi <- nrow(filtered_mi_figo_results_df)

    dset <- unified_dsets[dset_idx]
    mi_figo_meta_df <- mi_figo_meta_df %>%
        add_column(!!dset := c(n_mi))

    ## Baselines
    cls_baselines_df <- read_tsv(paste0(dirs$analysis_dir, "/meta/", "cls_baselines.tsv"))
    f1_macro_majority_baseline <- (cls_baselines_df %>% filter(dataset == unified_dsets[dset_idx]))$f1_macro_majority
    f1_macro_MC_baseline <- (cls_baselines_df %>% filter(dataset == unified_dsets[dset_idx]))$f1_macro_MC

    ## GBC (F1)
    f1_gbc_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_gbc_ref_scores.tsv"))
    f1_gbc_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_gbc_results.tsv"))
    f1_gbc_results_df <- get_consensus_col(f1_gbc_results_df, n = consensus_n, thresh = consensus_thresh)

    f1_gbc_avg <- mean(f1_gbc_scores_df$ref_score)
    f1_gbc_imp <- f1_gbc_avg > f1_macro_MC_baseline
    f1_gbc_pct_imp <- (f1_gbc_avg - f1_macro_MC_baseline) / f1_macro_MC_baseline * 100
    n_f1_gbc_consensus_genes <- nrow(f1_gbc_results_df %>% dplyr::filter(consensus == TRUE))

    dset <- unified_dsets[dset_idx]
    f1_gbc_meta_df <- f1_gbc_meta_df %>%
        add_column(!!dset := c(f1_gbc_avg, f1_gbc_imp, f1_gbc_pct_imp, n_f1_gbc_consensus_genes, f1_macro_MC_baseline))

    ## RFC (F1)
    f1_rfc_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_rfc_ref_scores.tsv"))
    f1_rfc_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_rfc_results.tsv"))
    f1_rfc_results_df <- get_consensus_col(f1_rfc_results_df, n = consensus_n, thresh = consensus_thresh)

    f1_rfc_avg <- mean(f1_rfc_scores_df$ref_score)
    f1_rfc_imp <- f1_rfc_avg > f1_macro_MC_baseline
    f1_rfc_pct_imp <- (f1_rfc_avg - f1_macro_MC_baseline) / f1_macro_MC_baseline * 100
    n_f1_rfc_consensus_genes <- nrow(f1_rfc_results_df %>% dplyr::filter(consensus == TRUE))

    dset <- unified_dsets[dset_idx]
    f1_rfc_meta_df <- f1_rfc_meta_df %>%
        add_column(!!dset := c(f1_rfc_avg, f1_rfc_imp, f1_rfc_pct_imp, n_f1_rfc_consensus_genes, f1_macro_MC_baseline))

    ## L1 LR (F1)
    f1_l1_lr_scores_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_l1_lr_ref_scores.tsv"))
    f1_l1_lr_results_df <- read_tsv(paste0(dirs$analysis_dir, "/feature_selection/", unified_dsets[dset_idx], "_l1_lr_results.tsv"))
    f1_l1_lr_results_df <- get_consensus_col(f1_l1_lr_results_df, n = consensus_n, thresh = consensus_thresh)

    f1_l1_lr_avg <- mean(f1_l1_lr_scores_df$ref_score)
    f1_l1_lr_imp <- f1_l1_lr_avg > f1_macro_MC_baseline
    f1_l1_lr_pct_imp <- (f1_l1_lr_avg - f1_macro_MC_baseline) / f1_macro_MC_baseline * 100
    n_f1_l1_lr_consensus_genes <- nrow(f1_l1_lr_results_df %>% dplyr::filter(consensus == TRUE))

    dset <- unified_dsets[dset_idx]
    f1_l1_lr_meta_df <- f1_l1_lr_meta_df %>%
        add_column(!!dset := c(f1_l1_lr_avg, f1_l1_lr_imp, f1_l1_lr_pct_imp, n_f1_l1_lr_consensus_genes, f1_macro_MC_baseline))
    #endregion
}

#region Write out meta DFs
write_tsv(deg_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_deg_meta.tsv"))
write_tsv(coxph_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_coxph_meta.tsv"))
write_tsv(cor_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_cor_meta.tsv"))
write_tsv(mi_survival_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_mi_survival_meta.tsv"))
write_tsv(mse_gbr_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_mse_gbr_meta.tsv"))
write_tsv(mse_rfr_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_mse_rfr_meta.tsv"))
write_tsv(anova_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_anova_meta.tsv"))
write_tsv(mi_figo_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_mi_figo_meta.tsv"))
write_tsv(f1_gbc_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_f1_gbc_meta.tsv"))
write_tsv(f1_rfc_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_f1_rfc_meta.tsv"))
write_tsv(f1_l1_lr_meta_df %>% rutils::transpose_df("item", "dataset"), paste0(dirs$analysis_dir, "/meta/", "fs_f1_l1_lr_meta.tsv"))
#endregion
