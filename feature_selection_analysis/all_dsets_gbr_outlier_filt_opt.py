import pandas as pd
import numpy as np
import os
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import mean_absolute_error
from skopt.space import Real, Integer, Categorical
from skopt import gp_minimize

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt


# Helper functions
def objective(h_params, X, y, loss_default, scoring_default, r, verbose=True):
    if verbose:
        print(h_params)
    model = GradientBoostingRegressor(
        # We use lad since it most closely matches up with hyper-parameter objective
        loss=loss_default,
        # Use this for the same reason
        learning_rate=h_params[0],
        n_estimators=h_params[1],
        max_depth=h_params[2],
        max_features=h_params[3],
        min_samples_split=h_params[4],
        min_samples_leaf=h_params[5],
        random_state=r
    )
    return -np.mean(cross_val_score(model, X, y, cv=KFold(n_splits=5), n_jobs=-1, scoring=scoring_default))


def run_optimization(x_df, y_df, space, loss_default, scoring_default, rand, n_initial, n_calls, callback_file):
    try:
        os.remove(callback_file)
    except OSError:
        pass

    res = gp_minimize(
        lambda h_ps: objective(h_ps, x_df, y_df.values.squeeze(), loss_default, scoring_default, rand),
        space,
        verbose=True,
        random_state=rand,
        n_initial_points=n_initial,
        n_calls = n_calls,
        n_jobs=-1,
        callback=lambda x: opt.save_callback(x, callback_file, n = 5, sep="\t")
    )


# Define constants
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["figo_stage", "age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["vital_status", "survival_time"]
cat_cols = ["race", "ethnicity", "figo_chr"]
space = [
    Real(1e-3, 1e-1, name="learning_rate"),
    Integer(int(1e2), int(1e3), name="n_estimators"),
    Integer(2, 5, name="max_depth"),
    Categorical(["auto", "sqrt", "log2"], name="max_features"),
    Integer(int(2), int(6), name="min_samples_split"),
    Integer(int(1), int(3), name="min_samples_leaf")
]
n_initial = 10 * (len(space))
n_calls = 50 * (len(space))


def main():
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
        filtered_joined_df = prep.filter_outliers_IQR(joined_df, "survival_time", coef=1.5)

        rand.seed(seed)
        x_df, y_df = prep.shuffle_data(filtered_joined_df, rand)

        # Optimize models
        run_optimization(
            x_df, y_df, space, "ls", "neg_mean_squared_error", rand, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_gbr_h_params_neg_mean_squared_error_removed.tsv"
        )

        print(f"Completed dataset: {unified_dsets[dset_idx]}")


if __name__ == "__main__":
    main()
