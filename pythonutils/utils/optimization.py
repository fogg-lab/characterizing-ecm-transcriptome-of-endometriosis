import pandas as pd
import numpy as np
from scipy.optimize import OptimizeResult
from sklearn.metrics import mean_absolute_error, explained_variance_score, f1_score
from sklearn.model_selection import cross_val_score, KFold
from sklearn.inspection import permutation_importance
from numpy.random import RandomState
from typing import Optional, Tuple, List


scoring_metrics = {
    "neg_mean_absolute_error": lambda y_true, y_pred: -mean_absolute_error(y_true, y_pred),
    "explained_variance": explained_variance_score,
    "f1_weighted": lambda y_true, y_pred: f1_score(y_true, y_pred, average="weighted")
}


def save_callback(res: OptimizeResult, dest: str, n: int = 5, sep: str = "\t") -> None:
    # Hyper-parameters from most recent model
    new_h_param_df = pd.DataFrame(
        columns=res.space.dimension_names + ["loss_achieved"],
        data=[res.x_iters[-1] + [res.func_vals[-1]]]
    )
    
    try:
        # Existing optimal hyper-parameters
        h_param_df = pd.read_csv(dest, sep=sep)
        
        # # Are all of the hyper-parameters for the newest model already here?
        # if redundant_h_params(new_h_param_df, h_param_df):
        #     print("*****redundant*****")
        # Not yet n models? Add the new one
        if h_param_df.shape[0] < n:
            h_param_df = pd.concat([h_param_df, new_h_param_df], axis=0)
        # New model better than the worst of previous n best? Replace it
        elif res.func_vals[-1] < h_param_df.loss_achieved.max():
            h_param_df.iloc[[h_param_df.loss_achieved.argmax()]] = new_h_param_df.values
    except FileNotFoundError:
        # First model? Save results
        h_param_df = new_h_param_df
    
    h_param_df.to_csv(dest, sep=sep, index=False)


def redundant_h_params(new_h_params_df: pd.DataFrame, old_h_params_df: pd.DataFrame) -> bool:
    # Are there repeats?
    # We ignore "loss_achieved" because this may be stochastic. Identical
    # hyper-parameters can yield different loss.
    res_df = pd.merge(
        new_h_params_df.drop("loss_achieved", axis=1),
        old_h_params_df.drop("loss_achieved", axis=1),
        on=list(new_h_params_df.columns[:-1])
    )
    return res_df.shape[0] > 0


def cv_permutation_importance(estimator: object, x_df: pd.DataFrame, y_df: pd.DataFrame, metric: str, k: int = 5, random_state: Optional[RandomState] = None, n_repeats: int = 11, to_array: bool = True) -> Tuple[List, np.ndarray]:
    metric_fun = scoring_metrics[metric]

    kf = KFold(n_splits=k)
    results = []
    ref_scores = []

    for train_idx, test_idx in kf.split(x_df):
        if to_array:
            y_train = y_df.iloc[train_idx].values.squeeze()
            y_test = y_df.iloc[test_idx].values.squeeze()
            x_train = x_df.iloc[train_idx].values
            x_test = x_df.iloc[test_idx].values
        else:
            y_train = y_df.iloc[train_idx].values.squeeze()
            y_test = y_df.iloc[test_idx].values.squeeze()
            x_train = x_df.iloc[train_idx]
            x_test = x_df.iloc[test_idx]

        # Train/retrain from scratch
        estimator.fit(x_train, y_train)
        result = permutation_importance(estimator, x_test, y_test, scoring=metric, n_jobs=-1, n_repeats=n_repeats, random_state=random_state)
        yhat = estimator.predict(x_test)
        ref_scores.append(metric_fun(y_test, yhat))
        results.append(result)
    return results, np.array(ref_scores)


def mc_classification_baseline(y: np.ndarray, labels: np.ndarray, weights: np.ndarray, metric: str, n: int = 11) -> np.ndarray:
    scores = []
    for _ in range(n):
        mc_yhat = np.random.choice(labels, size=y.shape[0], p = weights)
        scores.append(metric(y, mc_yhat))
    return np.array(scores)

