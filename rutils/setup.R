library(devtools)

tryCatch(
    expr = {
        devtools::uninstall()
    },
    finally = {
        devtools::load_all()
        devtools::install()
    }
)
