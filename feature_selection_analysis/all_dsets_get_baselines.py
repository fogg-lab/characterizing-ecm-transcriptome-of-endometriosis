import pandas as pd
import numpy as np
import os
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, explained_variance_score, f1_score
from skopt import gp_minimize

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt

# Define constants
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]


def main():
    seed = 123
    rand = np.random.RandomState()

    # Regression baselines
    event_code = {"Alive": 0, "Dead": 1}
    covariate_cols = ["figo_stage", "age_at_diagnosis", "race", "ethnicity"]
    dep_cols = ["vital_status", "survival_time"]
    cat_cols = ["race", "ethnicity", "figo_chr"]

    baseline_reg_df = pd.DataFrame({"baseline": ["L2", "L1", "R2", "explained_variance", "n"]})
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

        # Get baselines
        mean_baseline = mean_squared_error(y_df.values, np.repeat(np.mean(y_df.values.squeeze()), y_df.shape[0]))
        median_baseline = mean_absolute_error(y_df.values, np.repeat(np.median(y_df.values.squeeze()), y_df.shape[0]))
        r2_baseline = r2_score(y_df.values, np.repeat(np.mean(y_df.values.squeeze()), y_df.shape[0]))
        expl_var_baseline = explained_variance_score(y_df.values, np.repeat(np.mean(y_df.values.squeeze()), y_df.shape[0]))
        n = y_df.shape[0]

        print(f"******* Regression baselines for: {unified_dsets[dset_idx]} *******")
        print(f"L2 baseline: {mean_baseline}")
        print(f"L1 baseline: {median_baseline}")
        print(f"R2 baseline: {r2_baseline}")
        print(f"explained variance baseline: {expl_var_baseline}")
        print(f"Sample size: {n}")
        print()
        
        baseline_reg_df[f"{unified_dsets[dset_idx]}"] = np.array([mean_baseline, median_baseline, r2_baseline, expl_var_baseline, n])

    baseline_reg_t_df = prep.transpose_df(
        baseline_reg_df, "baseline", "dataset"
    )
    baseline_reg_t_df.to_csv(f"{dirs.analysis_dir}/meta/reg_baselines.tsv", sep="\t", index=False)

    # Classification baselines
    event_code = {"Alive": 0, "Dead": 1}
    covariate_cols = ["age_at_diagnosis", "race", "ethnicity"]
    dep_cols = ["figo_stage"]
    cat_cols = ["race", "ethnicity"]

    baseline_cls_df = pd.DataFrame({"baseline": ["f1_weighted_majority", "f1_weighted_MC", "n"]})
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

        
        label_value_counts_df = (
            pd.DataFrame(y_df.figo_num.value_counts()).reset_index()
                .rename(columns={"index": "label", "figo_num": "n"})
                .sort_values("n", ascending=False)
        )

        # Get baselines
        most_frequent_label = label_value_counts_df.label[0]
        most_frequent_baseline = f1_score(y_df.values.squeeze(), np.repeat(most_frequent_label, y_df.shape[0]), average="weighted")

        mc_baseline = opt.mc_classification_baseline(
            y=y_df.values.squeeze(),
            labels=label_value_counts_df.label.values,
            weights=label_value_counts_df.n.values / label_value_counts_df.n.values.sum(),
            metric=lambda y, yhat: f1_score(y, yhat, average="weighted"),
            n=1001
        )
        n = y_df.shape[0]

        print(f"******* Classification baselines for: {unified_dsets[dset_idx]} *******")
        print(f"F1 (weighted) majority guess baseline: {most_frequent_baseline}")
        print(f"F1 (weighted) Monte Carlo baseline: {mc_baseline.mean()}")
        print(f"Sample size: {n}")
        print()

        baseline_cls_df[f"{unified_dsets[dset_idx]}"] = np.array([most_frequent_baseline, mc_baseline.mean(), n])

    baseline_cls_t_df = prep.transpose_df(
        baseline_cls_df, "baseline", "dataset"
    )
    baseline_cls_t_df.to_csv(f"{dirs.analysis_dir}/meta/cls_baselines.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
