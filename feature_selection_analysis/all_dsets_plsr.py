import pandas as pd
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold
from sklearn.metrics import mean_absolute_error, explained_variance_score
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor

import utils.dev_config as dev_conf
import utils.preprocessing as prep
import utils.feature_selection as feat_sel


# Helper functions
def run_optimization(x_df, y_df, n_components_range, scoring, dset_idx):
    plsr_model = PLSRegression(scale=False)
    # Does worse with pre-processing, so just going to use vanilla
    plsr_pipeline = make_pipeline(plsr_model)
    ttr = TransformedTargetRegressor(regressor=plsr_pipeline)
    h_params = {"regressor__plsregression__n_components": n_components_range}
    cv_grid_search = GridSearchCV(ttr, h_params, scoring=scoring, cv=KFold(5), n_jobs=-1, verbose=1)
    cv_grid_search.fit(x_df, y_df)
    best_plsr = cv_grid_search.best_estimator_.regressor_["plsregression"]
    cv_score = cross_val_score(cv_grid_search.best_estimator_, x_df, y_df, cv=KFold(5), scoring=scoring, n_jobs=-1)

    # Save optimal h_params
    plsr_h_params_df = pd.DataFrame({"n_components": [best_plsr.n_components], "loss_achieved": [cv_score.mean()]})
    plsr_h_params_df.to_csv(f"{unified_dsets[dset_idx]}_opt_plsr_h_params_{scoring}.tsv", sep="\t", index=False)

    # Save VIP results
    vip = feat_sel.VIP(best_plsr)
    non_gene_vars = list(x_df.drop(list(norm_filtered_matrisome_counts_t_df.columns[1:]), axis=1).columns)
    plsr_res_df = (
        pd.DataFrame({"geneID": x_df.columns, "vip_scores": vip, "coeff": best_plsr.coef_.squeeze()})
            .pipe(lambda x: x[~x.geneID.isin(non_gene_vars)])
            .reset_index(drop=True)
    )
    plsr_res_df.to_csv(f"{dirs.analysis_dir}/{unified_dsets[dset_idx]}_plsr_{scoring}_results.tsv", sep="\t", index=False)


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
matrisome_list = f"{dirs.data_dir}/matrisome/matrisome_hs_masterlist.tsv"
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["figo_stage", "age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["vital_status", "survival_time"]
cat_cols = ["race", "ethnicity", "figo_chr"]


def main():
    # Train models
    for dset_idx in range(3):
        # Load and filter survival data
        survival_df = prep.load_survival_df(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/survival_data.tsv", event_code)
        filtered_survival_df = (
            prep.decode_figo_stage(survival_df[["sample_name"] + dep_cols + covariate_cols].dropna(), to="c")
                .query("vital_status == 1")
                .drop(["vital_status"], axis=1)
                .pipe(pd.get_dummies, columns=cat_cols)
                .reset_index(drop = True)
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

        # Optimize models
        run_optimization(x_df, y_df, range(2, 20), "neg_mean_absolute_error", dset_idx)
        run_optimization(x_df, y_df, range(2, 20), "explained_variance", dset_idx)


if __name__ == "__main__":
    main()
