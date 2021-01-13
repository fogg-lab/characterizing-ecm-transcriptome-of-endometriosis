import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_regression
from dask import compute, delayed
from multiprocessing import freeze_support

import utils.dev_config as dev_conf
import utils.preprocessing as prep


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
matrisome_list = f"{dirs.data_dir}/matrisome/matrisome_hs_masterlist.tsv"

seed = 123
rand = np.random.RandomState()

event_code = {"Alive": 0, "Dead": 1}
covariate_cols = ["age_at_diagnosis", "race", "ethnicity"]
dep_cols = ["vital_status", "survival_time"]
cat_cols = ["race", "ethnicity"]


def main():
    for dset_idx in range(3):
        # Load and filter survival data
        survival_df = prep.load_survival_df(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/survival_data.tsv", event_code)
        filtered_survival_df = (
            survival_df[["sample_name"] + dep_cols]
                .query("vital_status == 1")
                .dropna()
                .reset_index(drop=True)
        )

        # Load normalized matrisome count data
        norm_matrisome_counts_df = pd.read_csv(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/norm_matrisome_counts.tsv", sep='\t')
        norm_filtered_matrisome_counts_t_df = prep.transpose_df(
            norm_matrisome_counts_df[["geneID"] + list(filtered_survival_df.sample_name)], "geneID", "sample_name"
        )

        # Combine survival data and normalized count data
        joined_df = (
            pd.merge(filtered_survival_df, norm_filtered_matrisome_counts_t_df, on="sample_name")
                .drop("vital_status", axis=1)
                .set_index("sample_name")
        )

        # Examine mutual information
        X = joined_df.iloc[:, 1:].values
        y = joined_df.iloc[:, 0].values

        rand.seed(seed)
        sim_rounds = 101
        mi_delayed = [delayed(mutual_info_regression)(X, y, discrete_features=False, random_state=rand) for _ in range(sim_rounds)]
        res = compute(*mi_delayed, scheduler="processes")

        mi_df = pd.concat([
            pd.DataFrame({"geneID": joined_df.columns[1:]}),
            pd.DataFrame(np.column_stack(res), columns=[f"MI_est_{i + 1}" for i in range(sim_rounds)])
        ], axis=1)
        mi_df["MI_est_median"] = mi_df.iloc[:, 1:].median(axis=1)

        # Save results
        mi_df[["geneID", "MI_est_median"]].to_csv(f"{dirs.analysis_dir}/feature_selection/{unified_dsets[dset_idx]}_MI_survival_results.tsv", sep="\t", index=False)
        
        print(f"Completed dataset: {unified_dsets[dset_idx]}")


if __name__ == "__main__":
    freeze_support()
    main()
