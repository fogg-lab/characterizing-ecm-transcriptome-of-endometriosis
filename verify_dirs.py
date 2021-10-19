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


#### DATA ####

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

data_array = ["matrisome", "raw_unified_TCGA_GTEx", "saved_network_objects", "saved_RSE_objects", "tcga_biolinks_downloads", "TCGA_RNA_combined_matrix_count_data", "TCGA_RNA_matrix_count_data", "unified_cervical_data", "unified_uterine_data", "unified_uterine_endometrial_data"]

data_length = len(data_array)

for i in range(data_length):
    if os.path.exists(f"{data_dir}/{data_array[i]}"):
        print(f"{SFM.success} Data : {data_array[i]} exists")
    else: 
        print(f"{SFM.failure} {data_array[i]} Data : does not exist*")
        os.makedirs(f"{data_dir}/{data_array[i]}")
        print(f"{SFM.success} Created data/{data_array[i]} {data_dir}/{data_array[i]}")
        print(f"\tWARNING: directory created, but data may be missing.")


#### ANALYSIS ####

# Check existence of analysis directory
if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)
    print(f"{SFM.success} Created analysis directory: {analysis_dir}")
else:
    print(f"{SFM.success} Analysis directory already exists.")

analysis_array = ["deg", "model_opt", "feature_selection", "meta", "network", "enrichment", "gene_lists", "gene_lists_extra", "gene_lists_extra_network", "survival", "manual_notes", "one_off"]

analysis_length = len(analysis_array)

for i in range(analysis_length):
    if os.path.exists(f"{analysis_dir}/{analysis_array[i]}"):
        print(f"{SFM.success} Analysis :  {analysis_array[i]} exists")
    else: 
        print(f"{SFM.failure} Analysis : {analysis_array[i]} does not exist*")
        os.makedirs(f"{analysis_dir}/{analysis_array[i]}")
        print(f"{SFM.success} Created analysis/{analysis_array[i]} directory: {analysis_dir}/{analysis_array[i]}")
        print(f"\tWARNING: directory created, but data may be missing.")


#### Figures ####

# Check existence of figures directory
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)
    print(f"{SFM.success} Created figures directory: {figures_dir}")
else:
    print(f"{SFM.success} Figures directory already exists.")

figures_array = ["TCGA_overview", "network", "deg", "enrichment", "gene_lists", "models", "one_off", "paneled", "saved_obj"]

figures_length = len(figures_array)

for i in range(figures_length):
    if os.path.exists(f"{figures_dir}/{figures_array[i]}"):
        print(f"{SFM.success} Figures : {figures_array[i]} exists")
    else: 
        print(f"{SFM.failure} Figures : {figures_array[i]} does not exist*")
        os.makedirs(f"{figures_dir}/{figures_array[i]}")
        print(f"{SFM.success} Created figures/{figures_array[i]} directory: {figures_dir}/{figures_array[i]}")
        print(f"\tWARNING: directory created, but data may be missing.")


#### Final check ####
if all_present:
    print(f"{SFM.all_succeeded} All necessary data/directories present.")
else:
    print(f"{SFM.failures_present} Check messages and make sure all needed directories/files are present.")
