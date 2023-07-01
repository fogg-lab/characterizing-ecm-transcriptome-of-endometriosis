# Characterizing the extracellular matrix transcriptome of endometriosis

https://doi.org/10.1101/2023.03.08.531805

## Setup

### Prerequisites

- R (tested with R 4.2)
- Python 3.7+ (tested with Python 3.9)
- Jupyter Notebook

### 1. Install dependencies (R and Python packages)

Please read through the R packages installation script (install_r_packages.r) before running it. Some of the listed R packages may require additional system dependencies.

R packages installed by the script, `install_r_packages.r`:
- affy
- sva
- readr
- dplyr
- Biobase
- BiocGenerics
- BiocParallel
- genefilter
- hgu133plus2cdf
- jsonlite
- tibble
- limma
- yaml
- ggrepel

If you have R set up to install packages system-wide, you may need to do one of the following things:
- Run `Rscript install_r_packages.r` as superuser/administrator.
- Install these packages to a personal library.
- Manually install the packages listed in the script.

**Run the following commands in a terminal:**

```
git clone https://github.com/fogg-lab/characterizing-endometriosis-transcriptome.git
cd characterizing-endometriosis-transcriptome
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
