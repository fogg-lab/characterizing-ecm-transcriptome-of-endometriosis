# Characterizing the extracellular matrix transcriptome of endometriosis

https://doi.org/10.1101/2023.03.08.531805

## Setup

### 1. Install dependencies (R and Python packages)

Please read through the R packages installation script (install_r_packages.r) before running it.
If you have R set up to install packages system-wide, you may need to do one of the following things:
- Run `Rscript install_r_packages.r` as superuser/admin
- Install these packages to a personal library
- Manually install the packages listed in the script

**Run the following commands in a terminal:**

```
git clone https://github.com/fogg-lab/endometriosis-microarray-analysis.git
cd endometriosis-microarray-analysis
pip install -r requirements.txt
Rscript install_r_packages.r
```

### 2. Prepare data for analysis

Run the Jupyter notebook, data_prep/prep.ipynb

## Unsupervised analysis (hierarchical clustering)

Run the Jupyter notebook, analysis/clustering.ipynb

## Condition stratification

Run the script:

```
cd analysis
python elasticnet_classification.py
```
