# List of packages
packages <- c("affy",
              "sva",
              "readr",
              "dplyr",
              "Biobase",
              "BiocGenerics",
              "BiocParallel",
              "genefilter",
              "hgu133plus2cdf",
              "jsonlite",
              "tibble",
              "limma",
              "yaml",
              "ggrepel",
              "clusterProfiler",
              "stringi",
              "org.Hs.eg.db",
              "car")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# below is used to avoid a common multithreading bug when calling the rma function
BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)

for(package in packages) {
    if(!require(package, character.only = TRUE)) {
        if(package %in% rownames(installed.packages()) == FALSE) {
            print(paste("Installing", package))
            if(package %in% c("BiocManager", "readr", "dplyr", "jsonlite", "tibble", "yaml", "ggrepel", "stringi", "car")) {
                install.packages(package, dependencies = TRUE)
            } else {
                BiocManager::install(package)
            }
        }
    } else {
        print(paste(package, "already installed"))
    }
}

# Install IRkernel for Jupyter notebook
install.packages('devtools')
devtools::install_github('IRkernel/IRkernel')

#IRkernel::installspec()  # to register the kernel in the current R installation
# Comment-out above line as it may not be able to see jupyter installed via package managers
# Instead, we will add a cell to notebooks that need the R kernel with the following:
# !R -e "IRkernel::installspec()"
