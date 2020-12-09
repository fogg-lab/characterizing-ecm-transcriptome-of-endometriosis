import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, KFold
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.metrics import f1_score
from skopt.space import Real, Integer, Categorical
from skopt import gp_minimize

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.optimization as opt


# Helper functions
def objective(h_params, X, y, scoring_default, r, verbose=True):
    if verbose:
        print(h_params)
    model = RandomForestClassifier(
        n_estimators=h_params[0],
        max_depth=h_params[1],
        max_features=h_params[2],
        min_samples_split=h_params[3],
        min_samples_leaf=h_params[4],
        bootstrap=h_params[5],
        n_jobs=-1,
        random_state=r
    )
    return -np.mean(cross_val_score(
        model,
        X,
        y,
        cv=KFold(n_splits=5),
        n_jobs=-1,
        scoring=scoring_default
    ))


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


# Define constants
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
condition_map = {"healthy": 0, "tumor": 1}

seed = 123
rand = np.random.RandomState()

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

scoring_method = "f1_macro"


def main():
    for dset_idx in range(3):
        norm_matrisome_counts_df = pd.read_csv(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/norm_matrisome_counts.tsv", sep='\t')
        coldata_df = pd.read_csv(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/coldata.tsv", sep="\t")
        counts_t_df = prep.transpose_df(
            norm_matrisome_counts_df, "geneID", "sample_name"
        )
        joined_df = (
            pd.merge(coldata_df, counts_t_df, on="sample_name")
                .set_index("sample_name")
                .drop("data_source", axis=1)
                .assign(condition = lambda df: df.condition.apply(lambda x: condition_map[x]))
        )

        rand.seed(seed)
        x_df, y_df = prep.shuffle_data(joined_df, rand)

        run_optimization(
            x_df, y_df, space, scoring_method, rand, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_rfc_cancer_y_n_h_params_{scoring_method}.tsv"
        )

        print(f"Completed dataset: {unified_dsets[dset_idx]}")



if __name__ == "__main__":
    main()