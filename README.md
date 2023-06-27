# endometriosis-microarray-analysis

Characterizing the extracellular matrix transcriptome of endometriosis

https://doi.org/10.1101/2023.03.08.531805

## Setup

```
git clone https://github.com/fogg-lab/endometriosis-microarray-analysis.git
cd endometriosis-microarray-analysis
pip install -r requirements.txt
```

## Data preparation

Run the Jupyter notebook, data_prep/prep.ipynb

## Unsupervised analysis (hierarchical clustering)

Run the Jupyter notebook, analysis/clustering.ipynb

## Condition stratification

```
cd analysis
python elasticnet_classification.py
```
