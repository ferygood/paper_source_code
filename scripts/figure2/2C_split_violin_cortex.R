library(introdataviz)
library(dplyr)
library(ggplot2)

exp <- readRDS("~/Desktop/phd/pBrain/data/primateBrain_kznfTEexp.rds")

# split gene and TE dataframe
kznf_exp <- exp[,c(1:285, 1191, 1192)]
kznf_exp_longer <- kznf_exp %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = c(1:285), names_to = "KRAB-ZNFs",values_to = "value")

te_exp <- exp[,c(286:1192)]
te_exp_longer <- te_exp %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = c(1:905), names_to = "TEs", values_to = "value")

# add age inference
kznf_infer <- read.csv("~/Desktop/phd/pBrain/data/kznf_bucket.csv")
kznf_exp_longer <- kznf_exp_longer %>%
    left_join(kznf_infer[,c(2,6)], by=c("KRAB-ZNFs"="external_gene_name")) %>%
    mutate(value = log(value + 1))

te_infer <- read.csv("~/Desktop/phd/pBrain/data/Dfam_TE_simiiformes.csv")
te_exp_longer <- te_exp_longer %>%
    mutate(age = ifelse(TEs %in% te_infer$NM, "young", "old")) %>%
    mutate(value = log(value + 1))

# split violin plot
kznf_cortex <- kznf_exp_longer %>% filter(region=="Primary & Secondary Cortices")
kznf_cortex <- kznf_cortex %>% filter(value>0)
kznf_cortex$species <-
    factor(kznf_cortex$species,
           levels = c("Human", "Chimpanzee", "Bonobo", "Macaque"))

kznf_cortex_v <- ggplot(kznf_cortex, aes(x=species, y=value, fill=age)) +
    introdataviz::geom_split_violin(alpha=.8)+
    geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
    stat_compare_means(label = "p.signif", label.x=1.5, method="wilcox.test") +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    labs(fill="evolutionary age") +
    theme_bw() +
    rotate_x_text(angle = 20) +
    xlab("") +
    ylab("logExp.") +
    ggtitle("KRAB-ZNFs in Primary and Secondary Cortices")

te_cortex <- te_exp_longer %>% filter(region=="Primary & Secondary Cortices")
te_cortex <- te_cortex %>% filter(value>0)
te_cortex$species <-
    factor(te_cortex$species,
           levels = c("Human", "Chimpanzee", "Bonobo", "Macaque"))

te_cortex_v <- ggplot(te_cortex, aes(x=species, y=value, fill=age)) +
    introdataviz::geom_split_violin(alpha=.8)+
    geom_boxplot(width = .2, alpha = .6, show.legend = FALSE) +
    stat_compare_means(label = "p.signif", label.x=1.5, method="wilcox.test") +
    scale_fill_manual(values = c("#5A8FBB", "#E59E00")) +
    labs(fill="evolutionary age") +
    theme_bw() +
    rotate_x_text(angle = 20) +
    xlab("") +
    ylab("logExp.") +
    ggtitle("TEs in Primary and Secondary Cortices")

cortex <- ggarrange(kznf_cortex_v, te_cortex_v, ncol=1,
                    common.legend = TRUE, legend = "bottom")

ggsave(filename="figures/JPG/kznfs_TEs_violin_cortex.jpg", dpi=500, width=7, height=9)
ggsave(filename="figures/SVG/kznfs_TEs_violin_cortex.svg", dpi=500, width=7, height=9)
