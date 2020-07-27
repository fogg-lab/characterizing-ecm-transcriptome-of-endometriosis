

get_dev_directories <- function(dev_paths_file) {
    fn <- dev_paths_file
    con <- file(fn, open =  "r")
    lines <- readLines(con)
    for (i in seq_len(length(lines))) {
        print(lines[i])
    }
}