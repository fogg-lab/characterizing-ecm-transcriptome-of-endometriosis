# Characterizing the extracellular matrix transcriptome of endometriosis

https://link.springer.com/article/10.1007/s43032-023-01359-w

## Setup

### Prerequisites

- **Jupyter Notebook**
- **Python 3.7+**
- **R 4.2+**

### 1. Install dependencies (R and Python packages)

The following R packages are installed automatically by the script, `install_r_packages.r`:
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
- org.Hs.eg.db
- stringi
- tibble
- limma
- yaml
- ggrepel
- devtools
- IRkernel
- clusterProfiler

Some of the listed R packages may require additional system dependencies.

If you have R set up to install packages system-wide (rather than to a personal user library), you can either run the install script as admin/superuser, or manually install the packages listed above (note that IRkernel is installed via `devtools::install_github('IRkernel/IRkernel')`).

**Setup: Run the following commands at the command line:**

```zsh
git clone https://github.com/fogg-lab/characterizing-ecm-transcriptome-of-endometriosis.git
cd characterizing-ecm-transcriptome-of-endometriosis
pip install -r requirements.txt
Rscript install_r_packages.r
```

### 2. Prepare data for analysis

Run the Jupyter notebook, data_prep/prep.ipynb

## Unsupervised analysis (hierarchical clustering)

Run the Jupyter notebook, analysis/clustering.ipynb

## Condition stratification

Run the script:

```zsh
cd analysis
python regression classifier.py
```

## Compile condition stratification results and generate figures

Run the Jupyter notebook, analysis/get_classification_results.ipynb

## Enrichment analysis

Run the Jupyter notebook, analysis/enrichment_analysis.ipynb

## Differential expression analysis

Run the script, analysis/dgea.R

**Usage**

```zsh
Rscript dgea.R <counts_filepath> <coldata_filepath> <config_filepath> [<filter_filepath>] <output_dir>
```

**Example - Performing differential gene expression analysis with a filter list**

In this example, we will run the `dgea.R` script with the following parameters:

- `counts_filepath`: The file `all_phases_all_genes_counts.tsv` contains count data. 
- `coldata_filepath`: The file `all_phases_coldata.tsv` contains sample conditions, e.g. healthy/endometriosis.
- `config_filepath`: The YAML configuration file `dgea_config.yaml` is used.
- `filter_filepath` (optional argument): We are using the filter list `core_matrisome_genes.json`.
- `output_dir`: The results will be written to the `dgea_output` directory.

The command would be as follows:

```zsh
Rscript analysis/dgea.R data/all/all_phases_all_genes_counts.tsv data/all/all_phases_coldata.tsv analysis/dgea_config.yaml analysis/core_matrisome_genes.json dgea_output
```

**Command-line arguments (listed in positional order) for dgea.R**
- `-h` or `-help`: Print usage information and exit.
- `counts_filepath`: Path to the file containing count data.
- `coldata_filepath`: Path to the file containing column data.
- `config_filepath`: Path to the YAML file containing configuration settings.
- `filter_filepath`: (Optional) Path to the JSON file containing gene filter list.
- `output_dir`: Directory where the output file will be written.
