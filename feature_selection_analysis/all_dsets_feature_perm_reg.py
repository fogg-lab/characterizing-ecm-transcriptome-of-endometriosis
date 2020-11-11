import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, explained_variance_score
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt
import utils.feature_selection as feat_sel


# Helper functions
def collect_feature_perm_results(models, x_df, y_df, r, gene_cols, score, verbose=True):
    all_mean_perm_results = []
    all_ref_scores = []
    all_perm_res_dfs = []
    
    for i, m in enumerate(models):
        if verbose:
            print(f"Running feature perm for model {i}")
        perm_results, ref_scores = opt.cv_permutation_importance(m, x_df, y_df, score, k=5, random_state=r)
        perm_importances = np.concatenate([r.importances for r in perm_results], axis=1)
        perm_importance_means = np.mean(perm_importances, axis=1)
        
        all_mean_perm_results.append(perm_importance_means)
        all_ref_scores.append(ref_scores)
        
        res_df = feat_sel.gather_perm_res(x_df, perm_importance_means, np.mean(ref_scores), gene_cols)
        res_df = res_df.rename(columns={"mean_imp": f"mean_imp_{i}", "score_pct_improvement": f"score_pct_improvement_{i}"})
        all_perm_res_dfs.append(res_df)
    
    return all_mean_perm_results, all_ref_scores, all_perm_res_dfs


def merge_perm_results(perm_res_dfs):
    merge_df = perm_res_dfs[0]
    for i in range(1, len(perm_res_dfs)):
        merge_df = merge_df.merge(perm_res_dfs[i], on = "geneID", how = "inner")
    merge_df = (
        merge_df.assign(consensus_imp_mean = merge_df.filter(regex="mean_imp").mean(axis=1))
            .assign(consensus_imp_std = merge_df.filter(regex="mean_imp").std(axis=1))
    )
    merge_df = merge_df.assign(consensus_imp_cv = merge_df.consensus_imp_std / merge_df.consensus_imp_mean)
    merge_df["consensus_vote"] = (merge_df.set_index("geneID").filter(regex="mean_imp", axis=1) > 0).all(axis=1).values
    return merge_df


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
matrisome_list = f"{dirs.data_dir}/matrisome/matrisome_hs_masterlist.tsv"
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["figo_stage", "age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["vital_status", "survival_time"]
cat_cols = ["race", "ethnicity", "figo_chr"]


# Train models
for dset_idx in range(3):
    # Load and filter survival data
    survival_df = prep.load_survival_df(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/survival_data.tsv", event_code)
    filtered_survival_df = (
        prep.decode_figo_stage(survival_df[["sample_name"] + dep_cols + covariate_cols].dropna(), to="c")
            .query("vital_status == 1")
            .drop(["vital_status"], axis=1)
            .pipe(pd.get_dummies, columns=cat_cols)
            .reset_index(drop = True)
    )
    filtered_survival_df.columns = filtered_survival_df.columns.str.replace(' ', '_')

    # Load normalized matrisome count data
    norm_matrisome_counts_df = pd.read_csv(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/norm_matrisome_counts.tsv", sep='\t')
    norm_filtered_matrisome_counts_t_df = prep.transpose_df(
        norm_matrisome_counts_df[["geneID"] + list(filtered_survival_df.sample_name)], "geneID", "sample_name"
    )

    # Combine survival data and normalized count data
    joined_df = (
        pd.merge(filtered_survival_df, norm_filtered_matrisome_counts_t_df, on="sample_name")
            .set_index("sample_name")
    )

    rand.seed(seed)
    x_df, y_df = prep.shuffle_data(joined_df, rand)

    # Build models
    ev_gbr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_gbr_h_params_explained_variance.tsv", sep="\t")
    ev_gbrs = [
        GradientBoostingRegressor(
            **dict(zip(ev_gbr_h_param_df.columns[:-1], ev_gbr_h_param_df.iloc[i, :-1])), loss="ls", random_state=rand
        ) for i in range(ev_gbr_h_param_df.shape[0])
    ]

    mae_gbr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_gbr_h_params_neg_mean_absolute_error.tsv", sep="\t")
    mae_gbrs = [
        GradientBoostingRegressor(
            **dict(zip(mae_gbr_h_param_df.columns[:-1], mae_gbr_h_param_df.iloc[i, :-1])), loss="lad", random_state=rand
        ) for i in range(mae_gbr_h_param_df.shape[0])
    ]

    ev_rfr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_rfr_h_params_explained_variance.tsv", sep="\t")
    ev_rfrs = [
        RandomForestRegressor(
            **dict(zip(ev_rfr_h_param_df.columns[:-1], ev_rfr_h_param_df.iloc[i, :-1])), random_state=rand
        ) for i in range(ev_rfr_h_param_df.shape[0])
    ]

    mae_rfr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_rfr_h_params_neg_mean_absolute_error.tsv", sep="\t")
    mae_rfrs = [
        RandomForestRegressor(
            **dict(zip(mae_rfr_h_param_df.columns[:-1], mae_rfr_h_param_df.iloc[i, :-1])), random_state=rand
        ) for i in range(mae_rfr_h_param_df.shape[0])
    ]

    # GBR (EV)
    ev_gbr_mean_perm_res, ev_gbr_ref_scores, ev_gbr_perm_res_dfs = collect_feature_perm_results(
        ev_gbrs, x_df, y_df, rand, norm_filtered_matrisome_counts_t_df.columns[1:], "explained_variance"
    )
    ev_gbr_merge_df = merge_perm_results(ev_gbr_perm_res_dfs)
    ev_gbr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_ev_gbr_results.tsv", sep="\t", index=False)

    # GBR (MAE)
    mae_gbr_mean_perm_res, mae_gbr_ref_scores, mae_gbr_perm_res_dfs = collect_feature_perm_results(
        mae_gbrs, x_df, y_df, rand, norm_filtered_matrisome_counts_t_df.columns[1:], "neg_mean_absolute_error"
    )
    mae_gbr_merge_df = merge_perm_results(mae_gbr_perm_res_dfs)
    mae_gbr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_mae_gbr_results.tsv", sep="\t", index=False)

    # RFR (EV)
    ev_rfr_mean_perm_res, ev_rfr_ref_scores, ev_rfr_perm_res_dfs = collect_feature_perm_results(
        ev_rfrs, x_df, y_df, rand, norm_filtered_matrisome_counts_t_df.columns[1:], "explained_variance"
    )
    ev_rfr_merge_df = merge_perm_results(ev_rfr_perm_res_dfs)
    ev_rfr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_ev_rfr_results.tsv", sep="\t", index=False)

    # RFR (MAE)
    mae_rfr_mean_perm_res, mae_rfr_ref_scores, mae_rfr_perm_res_dfs = collect_feature_perm_results(
        mae_rfrs, x_df, y_df, rand, norm_filtered_matrisome_counts_t_df.columns[1:], "neg_mean_absolute_error"
    )
    mae_rfr_merge_df = merge_perm_results(mae_rfr_perm_res_dfs)
    mae_rfr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_mae_rfr_results.tsv", sep="\t", index=False)

    print(f"Completed dataset: {unified_dsets[dset_idx]}")
