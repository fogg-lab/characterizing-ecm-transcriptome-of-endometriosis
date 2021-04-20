library(tidyverse)

# Custom package
library(rutils)

installed.packages() %>%
    as_tibble() %>%
    rename_with(~ tolower(.)) %>%
    filter(!(priority %in% c("base"))) %>%
    select(package, version) %>%
    write_csv("r_packages.csv")
