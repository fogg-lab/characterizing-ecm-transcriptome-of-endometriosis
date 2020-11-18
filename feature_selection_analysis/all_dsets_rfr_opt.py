import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, explained_variance_score
from skopt.space import Real, Integer, Categorical
from skopt import gp_minimize

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt


# Helper functions
def objective(h_params, X, y, scoring_default, r, verbose=True):
    if verbose:
        print(h_params)
    model = RandomForestRegressor(
        n_estimators=h_params[0],
        max_depth=h_params[1],
        max_features=h_params[2],
        min_samples_split=h_params[3],
        min_samples_leaf=h_params[4],
        bootstrap=h_params[5],
        n_jobs=-1,
        random_state=r
    )
    return -np.mean(cross_val_score(model, X, y, cv=KFold(n_splits=5), n_jobs=-1, scoring=scoring_default))


def run_optimization(x_df, y_df, space, scoring_default, rand, n_initial, n_calls, callback_file):
    try:
        os.remove(callback_file)
    except OSError:
        pass

    res = gp_minimize(
        lambda h_ps: objective(h_ps, x_df, y_df.values.squeeze(), scoring_default, rand),
        space,
        verbose=True,
        random_state=rand,
        n_initial_points=n_initial,
        n_calls = n_calls,
        n_jobs=-1,
        callback=lambda x: opt.save_callback(x, callback_file, n = 5, sep="\t")
    )


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["figo_stage", "age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["vital_status", "survival_time"]
cat_cols = ["race", "ethnicity", "figo_chr"]
space = [
    Integer(int(1e2), int(1e3), name="n_estimators"),
    Integer(int(10), int(100), name="max_depth"),
    Categorical(["auto", "sqrt", "log2"], name="max_features"),
    Integer(int(2), int(4), name="min_samples_split"),
    Integer(int(1), int(3), name="min_samples_leaf"),
    Categorical([True, False], name="bootstrap")
]
n_initial = 10 * len(space)
n_calls = 50 * len(space)


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

        rand.seed(seed)
        x_df, y_df = prep.shuffle_data(joined_df, rand)

        # Optimize models
        run_optimization(
            x_df, y_df, space, "explained_variance", rand, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_rfr_h_params_explained_variance.tsv"
        )

        run_optimization(
            x_df, y_df, space, "neg_mean_absolute_error", rand, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_rfr_h_params_neg_mean_absolute_error.tsv"
        )
        
        print(f"Completed dataset: {unified_dsets[dset_idx]}")


if __name__ == "__main__":
    main()
