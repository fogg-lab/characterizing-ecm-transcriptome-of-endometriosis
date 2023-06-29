import os
import numpy as np
import pandas as pd

from typing import Tuple
from pathlib import Path
from scipy.optimize import OptimizeResult
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from skopt import gp_minimize
from skopt.space import Real


DATA_DIR = Path(__file__).parent.parent / "data"
GENE_SETS = ["all_genes", "all_matrisome", "core_matrisome"]
PHASES = ["all_phases", "early_secretory", "mid_secretory", "proliferative"]

class DataPaths:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)

    def __call__(self, phase: str, gene_set: str, fit: bool = True) -> Tuple[str, str]:
        subdir = "fit" if fit else "test"
        counts_filename = f"{phase}_{gene_set}_counts.tsv"
        coldata_filename = f"{phase}_coldata.tsv"
        directory = self.data_dir / subdir

        return str(directory / coldata_filename), str(directory / counts_filename)


def shuffle_data(df: pd.DataFrame, random_state: np.random.RandomState) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    shuffled_df = df.sample(frac=1, random_state=random_state)
    x_df = shuffled_df.iloc[:, 1:]
    y_df = shuffled_df.iloc[:, [0]]

    return x_df, y_df, shuffled_df


def save_callback(res: OptimizeResult, dest: str, n: int = 5, sep: str = "\t") -> None:
    # Hyper-parameters from most recent model
    new_h_param_df = pd.DataFrame(
        columns=res.space.dimension_names + ["loss_achieved"],
        data=[res.x_iters[-1] + [res.func_vals[-1]]]
    )

    try:
        # Existing optimal hyper-parameters
        h_param_df = pd.read_csv(dest, sep=sep)

        # Are all of the hyper-parameters for the newest model already here?
        # if redundant_h_params(new_h_param_df, h_param_df):
        #     print("*****redundant*****")
        # Not yet n models? Add the new one
        if h_param_df.shape[0] < n:
            h_param_df = pd.concat([h_param_df, new_h_param_df], axis=0)
        # New model better than the worst of previous n best? Replace it
        elif res.func_vals[-1] < h_param_df.loss_achieved.max():
            h_param_df.iloc[[int(h_param_df.loss_achieved.argmax())]] = new_h_param_df.values
    except FileNotFoundError:
        # First model? Save results
        h_param_df = new_h_param_df

    h_param_df.to_csv(dest, sep=sep, index=False)


def print_section_header(phase, gene_set, fit, coldata_path, counts_path):
    print('\n')
    print("-" * 80)
    print(phase, gene_set, "fit" if fit else "test")
    print(f"{coldata_path=}")
    print(f"{counts_path=}")
    print("-" * 80)


def join_and_process_dataframes(coldata_df, counts_df, condition_map):
    joined_df = (
        pd.merge(coldata_df, counts_df, left_on="sample_name", right_index=True).set_index("sample_name")
            .assign(condition=lambda df: df.condition.apply(lambda x: condition_map[x]))
    )
    return joined_df


def compute_cross_val_score(model, X, y, scoring_method):
    """Compute mean cross-validation score with KFold splits."""
    score = -np.mean(cross_val_score(
        model,
        X,
        y,
        cv=KFold(n_splits=5),
        n_jobs=-1,
        scoring=scoring_method
    ))
    return score


def create_logistic_regression_model(h_params, penalty_default, r, c_transformer):
    """Create a pipeline with a logistic regression model."""
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
    return model


def objective(h_params, X, y, penalty_default, scoring_default, r, c_transformer, verbose=True):
    if verbose:
        print(f"Hyperparameters: {h_params}")

    # Create model
    model = create_logistic_regression_model(h_params, penalty_default, r, c_transformer)

    # Compute cross-validation score
    score = compute_cross_val_score(model, X, y, scoring_default)

    return score


def run_optimization(x_df, y_df, space, penalty_default, scoring_default, rand, genes, n_initial, n_calls, callback_file):
    # Remove existing callback file if it exists
    if os.path.exists(callback_file):
        os.remove(callback_file)

    # Define column transformer
    c_transformer = ColumnTransformer([
        ("standard", StandardScaler(), genes)
    ], remainder="passthrough")

    try:
        res = gp_minimize(
            lambda h_ps: objective(h_ps, x_df, y_df.values.squeeze(), penalty_default, scoring_default, rand, c_transformer),
            space,
            verbose=True,
            random_state=rand,
            n_initial_points=n_initial,
            n_calls=n_calls,
            n_jobs=-1,
            callback=lambda x: save_callback(x, callback_file, n=5, sep="\t")
        )
    except ValueError as e:
        print(f"Optimization error: {e}")


def optimize_train_evaluate(x_df, y_df, elasticnet_space, penalty_default, scoring_method, rand,
                            genes, n_initial, n_calls, callback_file):
    # Optimize models
    run_optimization(x_df, y_df, elasticnet_space, penalty_default, scoring_method, rand,
                     genes, n_initial, n_calls, callback_file)

    # Load best hyperparameters
    best_params_dict = load_best_hyperparameters(callback_file)

    # Train the model
    print("Done optimizing... training...")
    best_model = train_model(x_df, y_df, genes, best_params_dict, rand)

    return best_model


def load_best_hyperparameters(callback_file):
    best_params_df = pd.read_csv(callback_file, sep="\t")
    best_params = best_params_df.sort_values("loss_achieved").iloc[0]
    best_params_dict = {"C": best_params["C"], "l1_ratio": best_params["l1_ratio"]}

    return best_params_dict


def train_model(x_df, y_df, matrisome_genes, best_params_dict, rand):
    c_transformer = ColumnTransformer([("standard", StandardScaler(), matrisome_genes)], remainder="passthrough")
    best_model = make_pipeline(
        c_transformer,
        LogisticRegression(
            C=best_params_dict["C"],
            l1_ratio=best_params_dict["l1_ratio"],
            solver="saga",
            penalty="elasticnet",
            n_jobs=-1,
            random_state=rand,
        ),
    )
    best_model.fit(x_df, y_df.values.squeeze())

    return best_model


def evaluate_model(best_model, x_df, y_df, shuffled_df, counts_path, phase, gene_set, coldata_path):
    # evaluate model from previous iteration
    predictions = best_model.predict(x_df)
    score = balanced_accuracy_score(y_df, predictions)

    results_df = pd.DataFrame({
        'sample_name': shuffled_df.index,
        'true_label': y_df.squeeze(),
        'predicted_label': predictions
    })

    # Save results
    results_df.to_csv(counts_path.replace("counts.tsv", "results.tsv"), sep="\t", index=False)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~TEST RESULTS~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"{phase=}, {gene_set=}, {coldata_path=}, {counts_path=}")
    print(f"Test score: {score}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

def main():
    data_paths = DataPaths(DATA_DIR)

    # Mapping conditions to numerical values
    condition_map = {"healthy": 0, "endometriosis": 1}

    # Setup random state and seed
    rand = np.random.RandomState()
    seed = 123

    # Setup the hyperparameter space for optimization
    elasticnet_space = [
        Real(1e-1, 1e1, name="C"),
        Real(0, 1, name="l1_ratio")
    ]

    # Scoring method
    scoring_method = "balanced_accuracy"

    # Initial points and calls for optimization
    n_initial = 10 * (len(elasticnet_space) + 1)
    n_calls = 50 * (len(elasticnet_space) + 1)

    model_trained = True
    model_evaluated = True

    for phase in PHASES:
        for gene_set in GENE_SETS:
            assert model_evaluated
            model_trained = False
            model_evaluated = False
            for fit in [True, False]:
                coldata_path, counts_path = data_paths(phase, gene_set, fit)

                print_section_header(phase, gene_set, fit, coldata_path, counts_path)

                counts_df = pd.read_csv(counts_path, sep="\t")
                counts_df = counts_df.groupby('symbol').mean().reset_index()
                counts_df.set_index('symbol', inplace=True)
                counts_df = counts_df.transpose()  # Now sample names are the index, and gene symbols are columns

                coldata_df = pd.read_csv(coldata_path, sep="\t", index_col=0)
                coldata_df = coldata_df.drop(["phase", "batch"], axis=1)
                coldata_df['condition'] = coldata_df['condition'].map(condition_map)

                joined_df = pd.merge(coldata_df, counts_df, left_index=True, right_index=True)

                genes = counts_df.columns

                rand.seed(seed)
                x_df, y_df, shuffled_df = shuffle_data(joined_df, rand)

                if fit:
                    callback_file = counts_path.replace("counts.tsv", "output.tsv")

                    best_model = optimize_train_evaluate(x_df, y_df, elasticnet_space, "elasticnet",
                                                         scoring_method, rand, genes, n_initial,
                                                         n_calls, callback_file)

                else:
                    assert model_trained
                    evaluate_model(best_model, x_df, y_df, shuffled_df, counts_path, phase, gene_set, coldata_path)

if __name__ == "__main__":
    main()
