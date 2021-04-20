# Environment
## *conda* (*python* and base *R*)
Base *R*, *python*, and all *python* packages are installed via the `conda` utility. To set up the *conda* environment, run the following from the root directory
```
conda env create -f environment.yml
```
and follow instructions output in terminal.

### Troubleshooting
It is possible that a package will no longer be available from the original source. If *conda* environment creation fails when installing a specific package, it may be solvable by doing the following.

1. Find the line in the `environment.yml` file with the problem package, `scipy` for example.
2. The line will look something like this
```
  - scipy=1.6.0=py38hb2138dd_0
```
3. Change the line to look like this
```
  - scipy=1.6.0
```
4. Re-attempt *conda* environment creation.

## *R* packages
The base *R* installation is done through *conda* but most *R* packages are installed through *R*. The file `r_packages.csv` contains all installed packages & versions. These should be installable via `install.versions()` using the data saved in the `r_packages.csv` file.

# Project pathing
## Developer configuration file
Make sure you have a file, `dev_paths.txt`, in the root project directory laid out as follows (*do not include comments*):
```
root/path/for/storing/data      # location of data to be analyzed
root/path/for/storing/analysis  # location of analysis output (usually .csv files)
root/path/for/storing/figures   # location of figure and visualization output
<blank line>
```

## Directory verification script
Run
```
python verify_dirs.py
```
from the project root. This will verify that all needed files/directories exist. This script will also create all needed output directories for analysis and figures.

# Custom utilities
This project has custom packages (`rutils` and `pythonutils`) that must be installed in the *conda* environment. Installations will only occur within the *conda* environment and should not affect any other *R* or *python* installations on your system. 

## `rutils`
Change to the `rutils` directory and run
```
Rscript setup.R
```

All custom *R* utilities should now be installed in the *conda* environment.

## `pythonutils`
Change to the `pythonutils` directory and run
```
python setup.py
```

All custom *python* utilities should now be installed in the *conda* environment.
