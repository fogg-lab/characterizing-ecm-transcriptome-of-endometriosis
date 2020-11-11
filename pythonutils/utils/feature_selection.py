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

# Adapted from Mehmood, T. et al.: https://www.sciencedirect.com/science/article/pii/S0169743912001542
def VIP(plsr_model):
    T = plsr_model.x_scores_
    W = plsr_model.x_weights_
    Q = plsr_model.y_loadings_
    p, a = W.shape
    
    vip = np.zeros(p)
    # SSa for each A
    SSA = np.sum(T ** 2, axis=0) * np.sum(Q ** 2, axis=0)
    # Column-wise l2 norm of W
    W_norm = np.einsum("ij, ij -> j", W, W)
    
    vip = np.sqrt(p * np.sum(SSA * (W / W_norm) ** 2, axis=1) / np.sum(SSA, axis=0))
    return vip
