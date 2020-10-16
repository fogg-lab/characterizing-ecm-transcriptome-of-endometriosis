import pandas as pd
from scipy.optimize import OptimizeResult


def save_callback(res: OptimizeResult, dest: str, sep: str = "\t") -> None:
    # Does this model correspond with the best score seen so far?
    if res.func_vals[-1] == res.fun:
        h_param_df = pd.DataFrame({
            "param": res.space.dimension_names + ["loss_achieved"],
            "param_value": res.x + [res.fun]
        })
        h_param_df.to_csv(dest, sep=sep, index=False)
