import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, KFold
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.metrics import f1_score
from sklearn.inspection import permutation_importance

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt
import utils.feature_selection as feat_sel


# Helper functions
def make_lr(h_params, matrisome_genes):
    if pd.isna(h_params["class_weight"]):
        h_params["class_weight"] = None
    model = make_pipeline(
        ColumnTransformer([
            ("standard", StandardScaler(), ["age_at_diagnosis"] + list(matrisome_genes))
        ], remainder="passthrough"),
        LogisticRegression(
            C=h_params["C"],
            class_weight=h_params["class_weight"],
            solver=h_params["solver"],
            penalty=h_params["penalty"],
            random_state=h_params["random_state"],
            n_jobs=-1
        )
    )
    return model


def collect_feature_perm_results(models, x_df, y_df, r, gene_cols, score, verbose=True, to_array=True):
    all_mean_perm_results = []
    all_ref_scores = []
    all_perm_res_dfs = []
    
    for i, m in enumerate(models):
        if verbose:
            print(f"Running feature perm for model {i}")
        perm_results, ref_scores = opt.cv_permutation_importance(m, x_df, y_df, score, k=5, random_state=r, to_array=to_array)
        perm_importances = np.concatenate([r.importances for r in perm_results], axis=1)
        perm_importance_means = np.mean(perm_importances, axis=1)
        
        all_mean_perm_results.append(perm_importance_means)
        all_ref_scores.append(ref_scores)
        
        res_df = feat_sel.gather_perm_res(x_df, perm_importance_means, np.mean(ref_scores), gene_cols)
        res_df = res_df.rename(columns={"mean_imp": f"mean_imp_{i}", "score_pct_improvement": f"score_pct_improvement_{i}"})
        all_perm_res_dfs.append(res_df)
    
    return all_mean_perm_results, all_ref_scores, all_perm_res_dfs


def merge_perm_results(perm_res_dfs, importance_thresh=0):
    merge_df = perm_res_dfs[0]
    for i in range(1, len(perm_res_dfs)):
        merge_df = merge_df.merge(perm_res_dfs[i], on = "geneID", how = "inner")
    # merge_df = (
    #     merge_df.assign(consensus_imp_mean = merge_df.filter(regex="mean_imp").mean(axis=1))
    #         .assign(consensus_imp_std = merge_df.filter(regex="mean_imp").std(axis=1))
    # )
    # merge_df = merge_df.assign(consensus_imp_cv = merge_df.consensus_imp_std / merge_df.consensus_imp_mean)
    # merge_df["consensus_vote"] = (merge_df.set_index("geneID").filter(regex="mean_imp", axis=1) > importance_thresh).all(axis=1).values
    return merge_df


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
matrisome_list = f"{dirs.data_dir}/matrisome/matrisome_hs_masterlist.tsv"
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["figo_stage"]
cat_cols = ["race", "ethnicity"]


scoring_method = "f1_macro"


def main():
    # Train models
    for dset_idx in range(3):
        # Load and filter survival data
        survival_df = prep.load_survival_df(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/survival_data.tsv", event_code)
        filtered_survival_df = (
            prep.decode_figo_stage(survival_df[["sample_name"] + dep_cols + covariate_cols].dropna(), to="n")
                .pipe(pd.get_dummies, columns=cat_cols)
                .reset_index(drop = True)
                .pipe(prep.cols_to_front, ["sample_name", "figo_num"])
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

        # Prep for running models
        matrisome_genes = norm_filtered_matrisome_counts_t_df.columns[1:]

        # Build models
        l1_lr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_lr_h_params_l1_{scoring_method}.tsv", sep="\t")
        l1_lrs = []
        for i in range(l1_lr_h_param_df.shape[0]):
            l1_lr_h_params = {
                **dict(zip(l1_lr_h_param_df.columns[:-1], l1_lr_h_param_df.iloc[i, :-1])), "penalty": "l1", "random_state": rand
            }
            l1_lrs.append(make_lr(l1_lr_h_params, matrisome_genes))

        # l2_lr_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_lr_h_params_l2_{scoring_method}.tsv", sep="\t")
        # l2_lrs = []
        # for i in range(l2_lr_h_param_df.shape[0]):
        #     l2_lr_h_params = {
        #         **dict(zip(l2_lr_h_param_df.columns[:-1], l2_lr_h_param_df.iloc[i, :-1])), "penalty": "l2", "random_state": rand
        #     }
        #     l2_lrs.append(make_lr(l2_lr_h_params, matrisome_genes))

        gbc_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_gbc_h_params_{scoring_method}.tsv", sep="\t")
        gbcs = [
            GradientBoostingClassifier(
                **dict(zip(gbc_h_param_df.columns[:-1], gbc_h_param_df.iloc[i, :-1])), loss="deviance", random_state=rand
            ) for i in range(gbc_h_param_df.shape[0])
        ]

        # rfc_h_param_df = pd.read_csv(f"{unified_dsets[dset_idx]}_opt_rfc_h_params_{scoring_method}.tsv", sep="\t")
        # rfcs = [
        #     RandomForestClassifier(
        #         **dict(zip(rfc_h_param_df.columns[:-1], rfc_h_param_df.iloc[i, :-1])), random_state=rand, n_jobs=-1
        #     ) for i in range(rfc_h_param_df.shape[0])
        # ]

        # LR (L1)
        l1_lr_mean_perm_res, l1_lr_ref_scores, l1_lr_perm_res_dfs = collect_feature_perm_results(
            l1_lrs, x_df, y_df, rand, matrisome_genes, scoring_method, to_array=False
        )
        l1_lr_merge_df = merge_perm_results(l1_lr_perm_res_dfs)
        l1_lr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_l1_lr_results.tsv", sep="\t", index=False)
        l1_lr_mean_ref_scores = np.array(l1_lr_ref_scores).mean(axis=1)
        l1_lr_mean_ref_scores_df = pd.DataFrame({"model": range(len(l1_lr_mean_ref_scores)), "ref_score": l1_lr_mean_ref_scores})
        l1_lr_mean_ref_scores_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_l1_lr_ref_scores.tsv", sep="\t", index=False)

        # # LR (L2)
        # l2_lr_mean_perm_res, l2_lr_ref_scores, l2_lr_perm_res_dfs = collect_feature_perm_results(
        #     l2_lrs, x_df, y_df, rand, matrisome_genes, scoring_method, to_array=False
        # )
        # l2_lr_merge_df = merge_perm_results(l2_lr_perm_res_dfs)
        # l2_lr_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_l2_lr_results.tsv", sep="\t", index=False)
        # l2_lr_mean_ref_scores = np.array(l2_lr_ref_scores).mean(axis=1)
        # l2_lr_mean_ref_scores_df = pd.DataFrame({"model": range(len(l2_lr_mean_ref_scores)), "ref_score": l2_lr_mean_ref_scores})
        # l2_lr_mean_ref_scores_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_l2_lr_ref_scores.tsv", sep="\t", index=False)

        # GBC
        gbc_mean_perm_res, gbc_ref_scores, gbc_perm_res_dfs = collect_feature_perm_results(
            gbcs, x_df, y_df, rand, matrisome_genes, scoring_method, to_array=True
        )
        gbc_merge_df = merge_perm_results(gbc_perm_res_dfs)
        gbc_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_gbc_results.tsv", sep="\t", index=False)
        gbc_mean_ref_scores = np.array(gbc_ref_scores).mean(axis=1)
        gbc_mean_ref_scores_df = pd.DataFrame({"model": range(len(gbc_mean_ref_scores)), "ref_score": gbc_mean_ref_scores})
        gbc_mean_ref_scores_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_gbc_ref_scores.tsv", sep="\t", index=False)

        # # RFC
        # rfc_mean_perm_res, rfc_ref_scores, rfc_perm_res_dfs = collect_feature_perm_results(
        #     rfcs, x_df, y_df, rand, matrisome_genes, scoring_method, to_array=True
        # )
        # rfc_merge_df = merge_perm_results(rfc_perm_res_dfs)
        # rfc_merge_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_rfc_results.tsv", sep="\t", index=False)
        # rfc_mean_ref_scores = np.array(rfc_ref_scores).mean(axis=1)
        # rfc_mean_ref_scores_df = pd.DataFrame({"model": range(len(rfc_mean_ref_scores)), "ref_score": rfc_mean_ref_scores})
        # rfc_mean_ref_scores_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_rfc_ref_scores.tsv", sep="\t", index=False)


        print(f"Completed dataset: {unified_dsets[dset_idx]}")


if __name__ == "__main__":
    main()
