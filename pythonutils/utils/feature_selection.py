import pandas as pd
import numpy as np
from typing import List


def gather_perm_res(x_df: pd.DataFrame, mean_importances: np.ndarray, ref_score_mean: np.ndarray, gene_cols: List) -> pd.DataFrame:
    perm_res_df = pd.DataFrame({
        "geneID": x_df.columns,
        "mean_imp": mean_importances
    })
    perm_res_df = (
        perm_res_df[perm_res_df.geneID.isin(gene_cols)]
            .assign(score_pct_improvement = lambda x: x.mean_imp / np.abs(ref_score_mean) * 100)
            .reset_index(drop=True)
    )
    return perm_res_df
