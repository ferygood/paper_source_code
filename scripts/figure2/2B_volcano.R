library(twice)
library(ggpubr)
library(ggplot2)
library(tidyverse)

data("hmKZNFs337")
data("hg19rmsk_info")

# import brain data
read_rds_files("../../paper_source_code/data/results_rdata/")

HmPtC1_resKZNFs <- HmPtC1$DEobject$gene_res %>%
    data.frame() %>%
    filter(rownames(.) %in% hmKZNFs337$external_gene_name) %>%
    mutate(log2FoldChange = log2FoldChange * -1)
HmPtC1_resTE <- HmPtC1$DEobject$te_res %>%
    data.frame() %>%
    mutate(log2FoldChange = log2FoldChange * -1)


geneVol <- volcano_plot(HmPtC1_resKZNFs)
teVol <- volcano_plot(HmPtC1_resTE)
gvol <- ggarrange(geneVol, teVol, nrow=1)
ggsave(filename="figures/JPG/2B_volcano.jpg", dpi=500, width=10, height=5)
ggsave(filename="figures/SVG/2B_volcano.svg", dpi=500, width=10, height=5)
