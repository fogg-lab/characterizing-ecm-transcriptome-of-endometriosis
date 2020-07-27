

get_dev_directories <- function(dev_paths_file) {
    fn <- dev_paths_file
    con <- file(fn, open =  "r")
    dirs <- c()
    fn_lines <- readLines(con)
    for (i in seq_len(length(fn_lines))) {
        dirs <- c(dirs, fn_lines[i])
    }
    dev_dirs <- list(data_dir = dirs[1], analysis_dir = dirs[2], figures_dir = dirs[3])
    class(dev_dirs) <- "DevDirs"
    return(dev_dirs)
}