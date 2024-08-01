split_violin_pbd <- function(brain_region){

    kznf_cortex <- kznf_exp_longer %>% filter(region==brain_region)
    kznf_cortex <- kznf_cortex %>% filter(value>0)
    kznf_cortex$species <-
        factor(kznf_cortex$species,
               levels = c("Human", "Chimpanzee", "Bonobo", "Macaque"))

    kznf_cortex_v <- ggplot(kznf_cortex, aes(x=species, y=value, fill=age)) +
        introdataviz::geom_split_violin(alpha=.8)+
        geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
        stat_compare_means(label = "p.signif", label.x=1.5, label.y=8.5, method="wilcox.test", size=8) +
        scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
        labs(fill="evolutionary age") +
        theme_bw() +
        rotate_x_text(angle = 20) +
        xlab("") +
        ylab("logExp.") +
        ggtitle(paste0("KRAB-ZNFs in ", brain_region)) +
        theme(
            plot.title = element_text(size = 20),      # title
            axis.title.x = element_text(size = 16),    # X axis
            axis.title.y = element_text(size = 16),    # Y axis
            axis.text.x = element_text(size = 20),     # X ticks
            axis.text.y = element_text(size = 20),     # Y ticks
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 14)
        ) +
        scale_y_continuous(limits = c(0,9))

    te_cortex <- te_exp_longer %>% filter(region==brain_region)
    te_cortex <- te_cortex %>% filter(value>0)
    te_cortex$species <-
        factor(te_cortex$species,
               levels = c("Human", "Chimpanzee", "Bonobo", "Macaque"))

    te_cortex_v <- ggplot(te_cortex, aes(x=species, y=value, fill=age)) +
        introdataviz::geom_split_violin(alpha=.8)+
        geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
        stat_compare_means(label = "p.signif", label.x=1.5, label.y=12.5, method="wilcox.test", size=8) +
        scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
        labs(fill="evolutionary age") +
        theme_bw() +
        rotate_x_text(angle = 20) +
        xlab("") +
        ylab("logExp.") +
        ggtitle(paste0("TEs in ", brain_region)) +
        theme(
            plot.title = element_text(size = 20),      # title
            axis.title.x = element_text(size = 16),    # X axis
            axis.title.y = element_text(size = 16),    # Y axis
            axis.text.x = element_text(size = 20),     # X ticks
            axis.text.y = element_text(size = 20),     # Y ticks
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 14)
        ) +
        scale_y_continuous(limits = c(0, 13))

    g_v <- ggarrange(kznf_cortex_v, te_cortex_v, ncol=1,
                     common.legend = TRUE, legend = "bottom")

    g_v
}


