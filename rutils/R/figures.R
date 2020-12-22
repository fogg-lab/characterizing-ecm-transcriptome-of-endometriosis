ea_plot <- function(df) {
    ggplot(df) +
        geom_bar(aes(x = reorder(type, (-log10(qval))), y = -log10(qval)), stat = "identity") +
        theme_classic() +
        theme(
            axis.title.y = element_blank()
        ) +
        labs(y = expression("-log"[10]*"(q-value)")) +
        coord_flip()
}

