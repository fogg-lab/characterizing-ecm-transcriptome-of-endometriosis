library(tidyverse)

ea_plot <- function(df, top = 40) {
    df <- df %>%
        dplyr::top_n(n = -top, wt = qval) %>%
        dplyr::arrange(-log10(qval)) %>%
        dplyr::mutate(type_fact = factor(type, levels = type))
    ggplot(df) +
        geom_bar(aes(x = type_fact, y = -log10(qval)), stat = "identity") +
        theme_classic() +
        theme(
            axis.title.y = element_blank()
        ) +
        labs(y = expression("-log"[10]*"(q-value)")) +
        coord_flip()
}
