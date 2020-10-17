import pandas as pd
import numpy as np
from scipy.optimize import OptimizeResult
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import cross_val_score, KFold
from sklearn.inspection import permutation_importance
from numpy.random import RandomState
from typing import Optional, Tuple, List


scoring_metrics = {
    "neg_mean_absolute_error": lambda y_true, y_pred: -mean_absolute_error(y_true, y_pred)
}


def save_callback(res: OptimizeResult, dest: str, sep: str = "\t") -> None:
    # Does this model correspond with the best score seen so far?
    if res.func_vals[-1] == res.fun:
        h_param_df = pd.DataFrame({
            "param": res.space.dimension_names + ["loss_achieved"],
            "param_value": res.x + [res.fun]
        })
        h_param_df.to_csv(dest, sep=sep, index=False)


def cv_permutation_importance(estimator: object, x_df: pd.DataFrame, y_df: pd.DataFrame, metric: str, k: int = 5, random_state: Optional[RandomState] = None, n_repeats: int = 11) -> Tuple[List, np.ndarray]:
    metric_fun = scoring_metrics[metric]

    kf = KFold(n_splits=k)
    results = []
    ref_scores = []

    for train_idx, test_idx in kf.split(x_df):
        # Train/retrain from scratch
        estimator.fit(x_df.iloc[train_idx], y_df.iloc[train_idx])
        result = permutation_importance(estimator, x_df.iloc[test_idx], y_df.iloc[test_idx], scoring=metric, n_jobs=-1, n_repeats=n_repeats, random_state=random_state)
        yhat = estimator.predict(x_df.iloc[test_idx])
        ref_scores.append(metric_fun(y_df.iloc[test_idx].values, yhat))
        results.append(result)
    return results, np.array(ref_scores)
