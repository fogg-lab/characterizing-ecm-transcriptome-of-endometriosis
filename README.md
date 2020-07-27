# Project pathing
## Developer configuration file
Make sure you have a file, "dev_paths.txt", in the root project directory laid out as follows (*do not include comments*):
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
from the project root. This will verify that all needed files/directories exist.

# Custom `R` utilities
This project has a package for `R` utilities (`rutils`) that must be installed in the `R` environment.

## Install `rutils` using `devtools`
1. Launch an `R` shell session from the `rutils` directory.
2. Execute
    ```
    library(devtools)
    devtools::load_all()
    devtools::install()
    ```
3. The package should now be usable via `library(rutils)`

## Uninstall `rutils` using `devtools`
1. Launch an `R` shell session from the `rutils` directory.
2. Execute
    ```
    library(devtools)
    devtools::uninstall()
    ```
3. The package should now be removed from the `R` environment.
