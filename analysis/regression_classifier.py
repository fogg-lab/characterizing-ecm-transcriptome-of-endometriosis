"""
Train and evaluate a logistic regression classifier to predict endometriosis condition.

Basic flow of the script:
1. Data Processing: Prepare gene expression data for model training and evaluation.
2. Model Training: Uses Bayesian Optimization to optimize hyperparameters of the model.
3. Model Evaluation: Once trained, the model is evaluated on a test dataset.

Usage:

    python elasticnet_classification.py
"""

import os
from typing import Tuple
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeResult
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from skopt import gp_minimize
from skopt.space import Real

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output" / "condition_stratification"
GENE_SETS = ["all_genes", "all_matrisome", "core_matrisome"]
PHASES = ["all_phases", "early_secretory", "mid_secretory", "proliferative"]


class DataPaths:
    """
    Helper class for managing data file paths in a structured data directory.

    The data directory is assumed to be organized into 'fit' and 'test' subdirectories, with files
    named according to the format "{phase}_{gene_set}_counts.tsv" and "{phase}_coldata.tsv".
    """

    def __init__(self, data_dir):
        """
        Initialize the DataPaths object with the given data directory.

        Args:
            data_dir (str): Data directory.
        """
        self._data_dir = Path(data_dir)

    def __call__(self, phase: str, gene_set: str, fit: bool = True) -> Tuple[str, str]:
        """
        Returns the file paths for the coldata and counts files for the given phase
        and gene set in the appropriate (fit or test) subdirectory.

        Args:
            phase (str): The phase of data to be fetched. This forms part of the filenames.
            gene_set (str): The gene set of data to be fetched. This forms part of the filenames.
            fit (bool, optional): Whether to get training data paths. Otherwise, returns test paths.
                                  Defaults to True.

        Returns:
            Tuple[str, str]: File paths for coldata and counts, respectively.
        """
        subdir = "fit" if fit else "test"
        counts_filename = f"{phase}_{gene_set}_counts.tsv"
        coldata_filename = f"{phase}_coldata.tsv"
        directory = self._data_dir / subdir

        return str(directory / coldata_filename), str(directory / counts_filename)


def shuffle_data(
        df: pd.DataFrame, random_state: np.random.RandomState
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Shuffles the rows of a DataFrame and separates it into input (x_df) and target (y_df) DataFrames.

    Assumes that the first column of the DataFrame is the target variable and the remaining columns
    are input features. The DataFrame is shuffled using the given random state for reproducibility.

    Args:
        df (pd.DataFrame): The input DataFrame to be shuffled and separated into input and target.
        random_state (np.random.RandomState): Random state to ensure reproducibility of shuffling.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: A tuple containing:
            - x_df (pd.DataFrame): The input features DataFrame.
            - y_df (pd.DataFrame): The target variable DataFrame.
            - shuffled_df (pd.DataFrame): The shuffled DataFrame.
    """
    shuffled_df = df.sample(frac=1, random_state=random_state)
    x_df = shuffled_df.iloc[:, 1:]
    y_df = shuffled_df.iloc[:, [0]]

    return x_df, y_df, shuffled_df


def save_callback(res: OptimizeResult, dest: str, n: int = 5, sep: str = "\t") -> None:
    """
    Saves hyperparameters from the most recent model's optimization result, res, to a file.

    If there are already n models' hyperparameters saved and the new model's performance is 
    better than at least one of them, the worst performing model's hyperparameters are replaced
    by the new model's hyperparameters.

    If the file doesn't exist, it means it's the first model and its results are saved directly.

    Args:
        res (OptimizeResult): The optimization result object from scikit-optimize.
        dest (str): The destination path to save the model's hyperparameters.
        n (int): The number of best models' hyperparameters to save. Defaults to 5.
        sep (str): The separator to use in the output file. Defaults to "\t".
    """
    new_h_param_df = pd.DataFrame(
        columns=res.space.dimension_names + ["loss_achieved"],
        data=[res.x_iters[-1] + [res.func_vals[-1]]]
    )

    try:
        h_param_df = pd.read_csv(dest, sep=sep)

        if h_param_df.shape[0] < n:
            h_param_df = pd.concat([h_param_df, new_h_param_df], axis=0)
        elif res.func_vals[-1] < h_param_df.loss_achieved.max():
            h_param_df.iloc[[int(h_param_df.loss_achieved.argmax())]] = new_h_param_df.values
    except FileNotFoundError:
        h_param_df = new_h_param_df

    h_param_df.to_csv(dest, sep=sep, index=False)



def print_section_header(phase, gene_set, fit, coldata_path, counts_path):
    """ Print a section header for the given phase, gene set, and fit/test status. """
    print('\n')
    print("-" * 80)
    print(phase, gene_set, "fit" if fit else "test")
    print(f"{coldata_path=}")
    print(f"{counts_path=}")
    print("-" * 80)


def compute_cross_val_score(model, X, y, scoring_method):
    """
    Computes the mean cross-validation score for a given model using K-Fold splitting.

    Args:
        model (sklearn.pipeline.Pipeline): The model to evaluate.
        X (numpy.ndarray or pandas.DataFrame): Feature matrix for model.
        y (numpy.ndarray or pandas.Series): Target vector for model.
        scoring_method (str or callable): A str (see sklearn model evaluation documentation) 
                                          or a scorer callable object / function with signature
                                          scorer(estimator, X, y).

    Returns:
        float: The mean cross-validation score.
    """
    score = -np.mean(cross_val_score(
        model,
        X,
        y,
        cv=KFold(n_splits=5),
        n_jobs=-1,
        scoring=scoring_method
    ))
    return score


def create_logistic_regression_model(h_params, penalty_default, r, c_transformer) -> Pipeline:
    """
    Creates a pipeline with a logistic regression model and a given column transformer.

    Args:
        h_params (list): Hyperparameters for the logistic regression model. Expects 2 values, 
                         the regularization strength 'C' and 'l1_ratio' for the Elastic-Net mixing parameter.
        penalty_default (str): The norm used in the penalization ('l1', 'l2', 'elasticnet', 'none').
        r (np.random.RandomState): RandomState instance used by np.random as source of randomness.
        c_transformer (ColumnTransformer): The column transformer to use in the pipeline.

    Returns:
        sklearn.pipeline.Pipeline: The pipeline with the logistic regression model and column transformer.
    """
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
    """
    Objective function for model optimization, builds a logistic regression model with given hyperparameters,
    and returns its cross-validation score.

    Args:
        h_params (list): Hyperparameters for the logistic regression model. Expects 2 values, 
                         the regularization strength 'C' and 'l1_ratio' for the Elastic-Net mixing parameter.
        X (pd.DataFrame): Input features for the model.
        y (pd.Series or array-like): Target variable for the model.
        penalty_default (str): The norm used in the penalization ('l1', 'l2', 'elasticnet', 'none').
        scoring_default (str): A string indicating the scoring metric for the cross-validation.
        r (np.random.RandomState): RandomState instance used by np.random as source of randomness.
        c_transformer (ColumnTransformer): Column transformer to use in the pipeline.
        verbose (bool, optional): Whether to print the hyperparameters during the optimization process. 
                                  Defaults to True.

    Returns:
        float: The mean cross-validation score of the model.
    """
    if verbose:
        print(f"Hyperparameters: {h_params}")

    # Create model
    model = create_logistic_regression_model(h_params, penalty_default, r, c_transformer)

    # Compute cross-validation score
    score = compute_cross_val_score(model, X, y, scoring_default)

    return score


def run_optimization(x_df, y_df, space, penalty_default, scoring_default, rand, genes, n_initial, n_calls, callback_file):
    """
    Run the optimization process for a logistic regression model, utilizing a Gaussian process.

    Args:
        x_df (pd.DataFrame): Input features for the model.
        y_df (pd.Series or array-like): Target variable for the model.
        space (list): The hyperparameters space to explore for the optimization.
        penalty_default (str): The norm used in the penalization ('l1', 'l2', 'elasticnet', 'none').
        scoring_default (str): A string indicating the scoring metric for the cross-validation.
        rand (np.random.RandomState): RandomState instance used by np.random as source of randomness.
        genes (list): List of column names to be standardized in the input features.
        n_initial (int): Number of evaluations of the objective function with initial points
                         before starting the sampling of the optimizer.
        n_calls (int): Total number of calls to the objective function.
        callback_file (str): Path to the file where callback results will be stored.
    """

    # Remove existing callback file if it exists
    if os.path.exists(callback_file):
        os.remove(callback_file)

    # Define column transformer
    c_transformer = ColumnTransformer([
        ("standard", StandardScaler(), genes)
    ], remainder="passthrough")

    gp_minimize(
        lambda h_ps: objective(h_ps, x_df, y_df.values.squeeze(), penalty_default, scoring_default, rand, c_transformer),
        space,
        verbose=True,
        random_state=rand,
        n_initial_points=n_initial,
        n_calls=n_calls,
        n_jobs=-1,
        callback=lambda x: save_callback(x, callback_file, n=5, sep="\t")
    )


def optimize_and_train_model(x_df, y_df, elasticnet_space, penalty_default, scoring_method, rand,
                             genes, n_initial, n_calls, callback_file):
    """
    Run the optimization process on a logistic regression model using given data and hyperparameters space, 
    then train the model with the best found hyperparameters.

    Args:
        x_df (pd.DataFrame): Input features for the model.
        y_df (pd.DataFrame): Target variable for the model.
        elasticnet_space (list): The hyperparameters space to explore for the ElasticNet optimization.
        penalty_default (str): The norm used in the penalization ('l1', 'l2', 'elasticnet', 'none').
        scoring_method (str): A string indicating the scoring metric for the cross-validation.
        rand (np.random.RandomState): RandomState instance used by np.random as source of randomness.
        genes (list): List of column names to be standardized in the input features.
        n_initial (int): Number of evaluations of the objective function with initial points before starting the 
                         sampling of the optimizer.
        n_calls (int): Total number of calls to the objective function.
        callback_file (str): Path to the file where callback results will be stored.

    Returns:
        best_model (Pipeline): Trained model using the best found hyperparameters.
    """

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
    """
    Load the best hyperparameters from a given file.

    Args:
        callback_file (str): Path to the file containing hyperparameters.

    Returns:
        best_params_dict (dict): Dictionary containing the best hyperparameters.
    """

    best_params_df = pd.read_csv(callback_file, sep="\t")
    best_params = best_params_df.sort_values("loss_achieved").iloc[0]
    best_params_dict = {"C": best_params["C"], "l1_ratio": best_params["l1_ratio"]}

    return best_params_dict


def train_model(x_df, y_df, matrisome_genes, best_params_dict, rand):
    """
    Train a logistic regression model using the best hyperparameters.

    Args:
        x_df (pd.DataFrame): DataFrame containing the input features.
        y_df (pd.DataFrame): DataFrame containing the target variable.
        matrisome_genes (list): List of genes for the column transformer.
        best_params_dict (dict): Dictionary containing the best hyperparameters.
        rand (int): Random seed for the logistic regression model.

    Returns:
        best_model (sklearn.Pipeline): Trained logistic regression model.
    """
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
    """
    Evaluate the model on a test dataset and save the results.

    Args:
        best_model (sklearn.Pipeline): Trained logistic regression model.
        x_df (pd.DataFrame): DataFrame containing the test features.
        y_df (pd.DataFrame): DataFrame containing the test target variable.
        shuffled_df (pd.DataFrame): DataFrame with shuffled indices.
        counts_path (str): Path of the counts file.
        phase (str): Current phase of the data pipeline.
        gene_set (str): Name of the gene set.
        coldata_path (str): Path of the coldata file.
    """

    # Evaluate model from previous iteration
    predictions = best_model.predict(x_df)
    score = balanced_accuracy_score(y_df, predictions)

    # Prepare results DataFrame
    results_df = pd.DataFrame({
        'sample_name': shuffled_df.index,
        'true_label': y_df.squeeze(),
        'predicted_label': predictions
    })

    # Generate save path
    filename = Path(counts_path).name.replace("counts.tsv", "test_results.tsv")
    save_path = OUTPUT_DIR / filename

    # Save results
    results_df.to_csv(save_path, sep="\t", index=False)

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~TEST RESULTS~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"{phase=}, {gene_set=}, {coldata_path=}, {counts_path=}")
    print(f"Test score: {score}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


def preprocess_data(counts_path, coldata_path, condition_map):
    """
    Load, preprocess, and merge the data from counts_path and coldata_path.

    Args:
        counts_path (str): Path to the counts file.
        coldata_path (str): Path to the column data file.
        condition_map (dict): Mapping from condition names to numerical values.

    Returns:
        genes (list): List of genes in the counts data.
        joined_df (DataFrame): Merged dataframe containing both the counts and column data.
    """

    # Load counts data
    counts_df = pd.read_csv(counts_path, sep="\t")
    counts_df = counts_df.groupby('symbol').mean().reset_index()
    counts_df.set_index('symbol', inplace=True)
    counts_df = counts_df.transpose()  # Now sample names are the index, and gene symbols are columns

    # Load column data (sample conditions)
    coldata_df = pd.read_csv(coldata_path, sep="\t", index_col=0)
    coldata_df = coldata_df.drop(["phase", "batch"], axis=1)
    coldata_df['condition'] = coldata_df['condition'].map(condition_map)

    # Merge column and counts data
    joined_df = pd.merge(coldata_df, counts_df, left_index=True, right_index=True)

    # Extract gene names
    genes = counts_df.columns

    return genes, joined_df


def main():
    """
    Main execution function. It performs the following steps:
    - Set up the environment (directories, mappings, seed, optimization parameters)
    - Iteratively process each dataset (for each phase and gene set)
    - Load and preprocess the data
    - Optimize the model hyperparameters and train the model if 'fit' is True
    - Otherwise, evaluate the previously trained model on the test set
    """

    data_paths = DataPaths(DATA_DIR)

    # Create output dir if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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

    for phase in PHASES:
        for gene_set in GENE_SETS:
            for fit in [True, False]:
                coldata_path, counts_path = data_paths(phase, gene_set, fit)

                print_section_header(phase, gene_set, fit, coldata_path, counts_path)

                # Data loading and preprocessing
                genes, joined_df = preprocess_data(counts_path, coldata_path, condition_map)

                # Shuffle the data
                rand.seed(seed)
                x_df, y_df, shuffled_df = shuffle_data(joined_df, rand)

                # Fit or evaluate the model
                if fit:
                    callback_filename = Path(counts_path).name.replace("counts.tsv", "fit_results.tsv")
                    callback_file = OUTPUT_DIR / callback_filename

                    best_model = optimize_and_train_model(x_df, y_df, elasticnet_space, "elasticnet",
                                                          scoring_method, rand, genes, n_initial,
                                                          n_calls, callback_file)

                else:
                    evaluate_model(best_model, x_df, y_df, shuffled_df, counts_path, phase, gene_set, coldata_path)


if __name__ == "__main__":
    main()
