import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LogisticRegression
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
def objective(h_params, X, y, penalty_default, scoring_default, r, c_transformer, verbose=True):
    if verbose:
        print(h_params)

    model = make_pipeline(
        c_transformer,
        LogisticRegression(
            C=h_params[0],
            l1_ratio=h_params[1],
            solver="saga",
            penalty=penalty_default,
            n_jobs=-1,
            random_state=r
        )
    )
    return -np.mean(cross_val_score(
        model,
        X,
        y,
        cv=KFold(n_splits=5),
        n_jobs=-1,
        scoring=scoring_default
    ))


def run_optimization(x_df, y_df, space, penalty_default, scoring_default, rand, matrisome_genes, n_initial, n_calls, callback_file):
    try:
        os.remove(callback_file)
    except OSError:
        pass
    c_transformer = ColumnTransformer([
        ("standard", StandardScaler(), matrisome_genes)
    ], remainder="passthrough")

    try:
        res = gp_minimize(
            lambda h_ps: objective(h_ps, x_df, y_df.values.squeeze(), penalty_default, scoring_default, rand, c_transformer),
            space,
            verbose=True,
            random_state=rand,
            n_initial_points=n_initial,
            n_calls = n_calls,
            n_jobs=-1,
            callback=lambda x: opt.save_callback(x, callback_file, n = 5, sep="\t")
        )
    # If get perfect scores, won't be able to learn & will throw error due to infs/NaNs
    except ValueError:
        pass


# Define constants
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
condition_map = {"healthy": 0, "tumor": 1}

seed = 123
rand = np.random.RandomState()

elasticnet_space = [
    Real(1e-1, 1e1, name="C"),
    Real(0, 1, name="l1_ratio")
]

n_initial = 10 * (len(elasticnet_space) + 1)
n_calls = 50 * (len(elasticnet_space) + 1)

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

        matrisome_genes = counts_t_df.columns[1:]

        rand.seed(seed)
        x_df, y_df = prep.shuffle_data(joined_df, rand)

        # Optimize models
        run_optimization(
            x_df, y_df, elasticnet_space, "elasticnet", scoring_method, rand, matrisome_genes, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_lr_cancer_y_n_h_params_elasticnet_{scoring_method}.tsv"
        )

        print(f"Completed dataset: {unified_dsets[dset_idx]}")



if __name__ == "__main__":
    main()