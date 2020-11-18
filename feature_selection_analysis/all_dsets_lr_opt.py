import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
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
            class_weight=h_params[1],
            solver=h_params[2],
            penalty=penalty_default,
            n_jobs=-1,
            random_state=r
        )
    )
    return -np.mean(cross_val_score(model, X, y, cv=KFold(n_splits=5), n_jobs=-1, scoring=scoring_default))


def run_optimization(x_df, y_df, space, penalty_default, scoring_default, rand, matrisome_genes, n_initial, n_calls, callback_file):
    try:
        os.remove(callback_file)
    except OSError:
        pass
    c_transformer = ColumnTransformer([
        ("standard", StandardScaler(), ["age_at_diagnosis"] + list(matrisome_genes))
    ], remainder="passthrough")

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


# Define constants
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["figo_stage"]
cat_cols = ["race", "ethnicity"]
l1_space = [
    Real(1e-1, 1e1, name="C"),
    Categorical(["balanced", None], name="class_weight"),
    Categorical(["liblinear", "saga"], name="solver")
]

l2_space = [
    Real(1e-1, 1e1, name="C"),
    Categorical(["balanced", None], name="class_weight"),
    Categorical(["newton-cg", "lbfgs", "sag", "saga"], name="solver")
]

no_penalty_space = [
    Real(1e-1, 1e1, name="C"),
    Categorical(["balanced", None], name="class_weight"),
    Categorical(["newton-cg", "lbfgs", "sag", "saga"], name="solver")
]

n_initial = 10 * len(l1_space)
n_calls = 50 * len(l1_space)


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

        # Optimize models
        run_optimization(
            x_df, y_df, l2_space, "l2", "f1_weighted", rand, matrisome_genes, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_lr_h_params_l2_f1_weighted.tsv"
        )

        run_optimization(
            x_df, y_df, l1_space, "l1", "f1_weighted", rand, matrisome_genes, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_lr_h_params_l1_f1_weighted.tsv"
        )

        run_optimization(
            x_df, y_df, no_penalty_space, "none", "f1_weighted", rand, matrisome_genes, n_initial, n_calls,
            f"{unified_dsets[dset_idx]}_opt_lr_h_params_none_f1_weighted.tsv"
        )

        print(f"Completed dataset: {unified_dsets[dset_idx]}")


if __name__ == "__main__":
    main()
