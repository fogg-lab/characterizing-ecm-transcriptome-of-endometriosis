import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture

import utils.dev_config as dev_conf
import utils.preprocessing as prep


# Define constants and load data
dirs = dev_conf.get_dev_directories("../dev_paths.txt")
unified_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]
matrisome_list = f"{dirs.data_dir}/matrisome/matrisome_hs_masterlist.tsv"
seed = 123
rand = np.random.RandomState()
event_code = {"Alive": 0, "Dead": 1}
dep_cols = ["vital_status", "survival_time"]


def main():
    for dset_idx in range(3):
        survival_df = prep.load_survival_df(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/survival_data.tsv", event_code)
        norm_matrisome_counts_df = pd.read_csv(f"{dirs.data_dir}/{unified_dsets[dset_idx]}/norm_matrisome_counts.tsv", sep='\t')

        filtered_survival_df = (
            survival_df[["sample_name", "survival_time", "vital_status"]]
        )
        norm_filtered_matrisome_counts_t_df = prep.transpose_df(
            norm_matrisome_counts_df[["geneID"] + list(filtered_survival_df.sample_name)], "geneID", "sample_name"
        )
        joined_df = (
            pd.merge(filtered_survival_df, norm_filtered_matrisome_counts_t_df, on="sample_name")
                .set_index("sample_name")
        )

        gene_names = list(joined_df.columns[2:])
        n_genes = len(gene_names)
        cutoffs = np.zeros(n_genes)

        rand.seed(seed)
        for i in range(n_genes):
            gene_i = gene_names[i]
            gm = GaussianMixture(n_components=2, random_state=rand, n_init=5)
            X = np.array(joined_df[gene_i])[:, np.newaxis]
            gm.fit(X)
            Yhat = gm.predict_proba(X)
            # If didn't converge, just choose the median (since no evidence of bi-modal dist.)
            if gm.converged_:
                cutoff_loc = np.argmin(abs(Yhat[:, 0] - Yhat[:, 1]))
                cutoffs[i] = X[cutoff_loc].item()
            else:
                cutoffs[i] = np.median(X)

        cutoff_df = pd.DataFrame({"geneID": gene_names, "cutoff": cutoffs})
        cutoff_df.to_csv(f"{dirs.analysis_dir}/survival/{unified_dsets[dset_idx]}_expression_cutoffs.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()
