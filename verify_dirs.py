import os
import sys
from dataclasses import dataclass


all_present = True
cervical_dir = "unified_cervical_data"
matrisome_dir = "matrisome"
cancer_dsets = ["unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]


dev_paths_file = "dev_paths.txt"
counts_file = "counts.tsv"
coldata_file = "coldata.tsv"
matrisome_masterlist_file = "matrisome_hs_masterlist.tsv"


@dataclass
class SFM:
    red = "\u001b[31m"
    green = "\u001b[32m"
    reset = "\u001b[0m"
    success = f"{green}[SUCCESS]{reset}"
    failure = f"{red}[FAILURE]{reset}"
    all_succeeded = f"{green}[ALL SUCCEEDED]{reset}"
    failures_present = f"{red}[FAILURES PRESENT]{reset}"


# Check existence of dev_paths_file
if not os.path.exists(dev_paths_file):
    print(f"{SFM.failure} Could not find: {dev_paths_file}")
    print(f"You must create this file with the specifications listed in README.md.")
    sys.exit("Please create this file and try again.")
else:
    print(f"{SFM.success} Found: {dev_paths_file}")

with open("dev_paths.txt") as f:
    dirs = tuple(fn.strip() for fn in f.readlines())
    data_dir, analysis_dir, figures_dir = dirs

# Check existence of cancer data
for dset in cancer_dsets:
    if not os.path.exists(f"{data_dir}/{dset}/{coldata_file}"):
        print(f"{SFM.failure} Could not find: {data_dir}/{dset}/{coldata_file}")
        all_present = False
    else:
        print(f"{SFM.success} Found: {data_dir}/{dset}/{coldata_file}")

    if not os.path.exists(f"{data_dir}/{dset}/{counts_file}"):
        print(f"{SFM.failure} Could not find: {data_dir}/{dset}/{counts_file}")
        all_present = False
    else:
        print(f"{SFM.success} Found: {data_dir}/{dset}/{counts_file}")


# Check existence of matrisome data
if not os.path.exists(f"{data_dir}/{matrisome_dir}/{matrisome_masterlist_file}"):
    print(f"{SFM.failure} Could not find: {data_dir}/{matrisome_dir}/{matrisome_masterlist_file}")
    all_present = False
else:
    print(f"{SFM.success} Found: {data_dir}/{matrisome_dir}/{matrisome_masterlist_file}")

# Check existence of pathology data
if not os.path.exists(f"{data_dir}/THPA_v20_1_staining/pathology.tsv"):
    print(f"{SFM.failure} Could not find: {data_dir}/THPA_v20_1_staining/pathology.tsv")
    all_present = False
else:
    print(f"{SFM.success} Found: {data_dir}/THPA_v20_1_staining/pathology.tsv")



# Check existence of analysis directory
if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)
    print(f"{SFM.success} Created analysis directory: {analysis_dir}")
else:
    print(f"{SFM.success} Analysis directory already exists.")

# Check existence of analysis/deg directory
if not os.path.exists(f"{analysis_dir}/deg"):
    os.makedirs(f"{analysis_dir}/deg")
    print(f"{SFM.success} Created analysis/deg directory: {analysis_dir}/deg")
else:
    print(f"{SFM.success} Analysis deg directory already exists.")

# Check existence of analysis/model_opt directory
if not os.path.exists(f"{analysis_dir}/model_opt"):
    os.makedirs(f"{analysis_dir}/model_opt")
    print(f"{SFM.success} Created analysis/model_opt directory: {analysis_dir}/model_opt")
else:
    print(f"{SFM.success} Analysis model_opt directory already exists.")

# Check existence of analysis/feature_selection directory
if not os.path.exists(f"{analysis_dir}/feature_selection"):
    os.makedirs(f"{analysis_dir}/feature_selection")
    print(f"{SFM.success} Created analysis/feature_selection directory: {analysis_dir}/feature_selection")
else:
    print(f"{SFM.success} Analysis feature_selection directory already exists.")

# Check existence of analysis/meta directory
if not os.path.exists(f"{analysis_dir}/meta"):
    os.makedirs(f"{analysis_dir}/meta")
    print(f"{SFM.success} Created analysis/meta directory: {analysis_dir}/meta")
else:
    print(f"{SFM.success} Analysis meta directory already exists.")

# Check existence of analysis/network directory
if not os.path.exists(f"{analysis_dir}/network"):
    os.makedirs(f"{analysis_dir}/network")
    print(f"{SFM.success} Created analysis/network directory: {analysis_dir}/network")
else:
    print(f"{SFM.success} Analysis network directory already exists.")

# Check existence of analysis/enrichment directory
if not os.path.exists(f"{analysis_dir}/enrichment"):
    os.makedirs(f"{analysis_dir}/enrichment")
    print(f"{SFM.success} Created analysis/enrichment directory: {analysis_dir}/enrichment")
else:
    print(f"{SFM.success} Analysis enrichment directory already exists.")

# Check existence of analysis/gene_lists directory
if not os.path.exists(f"{analysis_dir}/gene_lists"):
    os.makedirs(f"{analysis_dir}/gene_lists")
    print(f"{SFM.success} Created analysis/gene_lists directory: {analysis_dir}/gene_lists")
else:
    print(f"{SFM.success} Analysis gene_lists directory already exists.")

# Check existence of analysis/survival directory
if not os.path.exists(f"{analysis_dir}/survival"):
    os.makedirs(f"{analysis_dir}/survival")
    print(f"{SFM.success} Created analysis/survival directory: {analysis_dir}/survival")
else:
    print(f"{SFM.success} Analysis survival directory already exists.")


# Check existence of figures directory
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)
    print(f"{SFM.success} Created figures directory: {figures_dir}")
else:
    print(f"{SFM.success} Figures directory already exists.")

# Check existence of figures/TCGA_overview directory
if not os.path.exists(f"{figures_dir}/TCGA_overview"):
    os.makedirs(f"{figures_dir}/TCGA_overview")
    print(f"{SFM.success} Created figures/TCGA_overview directory: {figures_dir}/TCGA_overview")
else:
    print(f"{SFM.success} Figures TCGA_overview directory already exists.")

# Check existence of figures/network directory
if not os.path.exists(f"{figures_dir}/network"):
    os.makedirs(f"{figures_dir}/network")
    print(f"{SFM.success} Created figures/network directory: {figures_dir}/network")
else:
    print(f"{SFM.success} Figures network directory already exists.")

# Check existence of figures/deg directory
if not os.path.exists(f"{figures_dir}/deg"):
    os.makedirs(f"{figures_dir}/deg")
    print(f"{SFM.success} Created figures/deg directory: {figures_dir}/deg")
else:
    print(f"{SFM.success} Figures deg directory already exists.")

# Check existence of figures/enrichment directory
if not os.path.exists(f"{figures_dir}/enrichment"):
    os.makedirs(f"{figures_dir}/enrichment")
    print(f"{SFM.success} Created figures/enrichment directory: {figures_dir}/enrichment")
else:
    print(f"{SFM.success} Figures enrichment directory already exists.")

# Check existence of figures/gene_lists directory
if not os.path.exists(f"{figures_dir}/gene_lists"):
    os.makedirs(f"{figures_dir}/gene_lists")
    print(f"{SFM.success} Created figures/gene_lists directory: {figures_dir}/gene_lists")
else:
    print(f"{SFM.success} Figures gene_lists directory already exists.")

# Check existence of figures/models directory
if not os.path.exists(f"{figures_dir}/models"):
    os.makedirs(f"{figures_dir}/models")
    print(f"{SFM.success} Created figures/models directory: {figures_dir}/models")
else:
    print(f"{SFM.success} Figures models directory already exists.")


# Check existence of figures/one_off directory
if not os.path.exists(f"{figures_dir}/one_off"):
    os.makedirs(f"{figures_dir}/one_off")
    print(f"{SFM.success} Created figures/one_off directory: {figures_dir}/one_off")
else:
    print(f"{SFM.success} Figures one_off directory already exists.")


if all_present:
    print(f"{SFM.all_succeeded} All necessary data/directories present.")
else:
    print(f"{SFM.failures_present} Check messages and make sure all needed directories/files are present.")
