import os

all_present = True

with open("dev_paths.txt") as f:
    data_dir = f.readline().strip()

# Check existence of cervical cancer data
if not os.path.exists(f"{data_dir}/unified_cervical_data/coldata.tsv"):
    print(f"[ERROR] Could not find: {data_dir}/unified_cervical_data/coldata.tsv")
    all_present = False
if not os.path.exists(f"{data_dir}/unified_cervical_data/counts.tsv"):
    print(f"[ERROR] Could not find: {data_dir}/unified_cervical_data/counts.tsv")
    all_present = False

# Check existence of matrisome data
if not os.path.exists(f"{data_dir}/matrisome/matrisome_hs_masterlist.tsv"):
    print(f"[ERROR] Could not find: {data_dir}/matrisome/matrisome_hs_masterlist.tsv")
    all_present = False

if not os.path.exists(f"{data_dir}/analysis"):
    os.makedirs(f"{data_dir}/analysis")
    print(f"Created directory: {data_dir}/analysis")


if all_present:
    print(f"[SUCCESS] All necessary data/directories present.")
