# fig 3D
library(ComplexHeatmap)
library(TEKRABber)
library(twice)
library(tidyverse)
load("data/primateBrainData.RData")
data("hmKZNFs337")
data("hg19rmsk_info")

# prepare datasets including raw counts of KRAB-ZNFs and TEs
# convert to TPM
# genes
df_hm_gene <- hmGene[,c(-1)]
rownames(df_hm_gene) <- hmGene$geneID

# transposable elements
hsTEexp <- hmTE %>% select(-c(1,2,3))
rownames(hsTEexp) <- hmTE$name  #908 TEs

# genes convert to tpm
sample_counts <- colSums(df_hm_gene)

scaling_factor <- sample_counts / 1e6

df_hm_gene_tpm <- df_hm_gene
df_hm_gene_tpm <- t(t(df_hm_gene_tpm)/ scaling_factor + 1) * 1e6
df_hm_gene_tpm <- as.data.frame(df_hm_gene_tpm)

# tes convert to tpm
te_count <- colSums(hsTEexp)
te_scale <- te_count / 1e6
hsTE_tpm <- hsTEexp
hsTE_tpm <- t(t(hsTE_tpm)/ te_scale + 1) * 1e6
hsTE_tpm <- as.data.frame(hsTE_tpm)

hsKZNFexp <- df_hm_gene_tpm %>%
    mutate(geneName=rownames(.)) %>%
    inner_join(hmKZNFs337, join_by("geneName"=="ensembl_gene_id")) #337

rownames(hsKZNFexp) <- hsKZNFexp$external_gene_name

hsKZNFexp <- hsKZNFexp %>% select(-c(133, 134)) #keep only expression data

cluster_meta <- metadata %>%
    filter(Organism == "Homo sapiens") %>%
    inner_join(brain_meta, join_by("brain_region"=="region"))

# cluster_Corr: calculate correlation based on different brain cluster
cluster_Corr <- function(gene, te, cluster_num){

    cluster_id <- cluster_meta %>%
        filter(cluster == cluster_num) %>%
        select(1)

    cluster_gene <- gene %>% select(cluster_id$Run)
    cluster_te <- te %>% select(cluster_id$Run)

    st <- Sys.time()
    df_temp <- corrOrthologTE(
        geneInput = cluster_gene,
        teInput = cluster_te,
        numCore = 5
    )

    et <- Sys.time()
    print(et-st)

    df_temp <- df_temp %>%
        mutate(pair = paste0(teName, ":", geneName))

    df_temp
}

hsC1 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster1")
hsC2 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster2")
hsC3 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster3")
hsC4 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster4")
hsC5 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster5")
hsC6 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster6")
hsC7 <- cluster_Corr(hsKZNFexp, hsTE_tpm, "cluster7")

hsC1.sig <- hsC1 %>% filter(padj<0.01 & abs(coef) >= 0.4) #127797
hsC2.sig <- hsC2 %>% filter(padj<0.01 & abs(coef) >= 0.4) #49770
hsC3.sig <- hsC3 %>% filter(padj<0.01 & abs(coef) >= 0.4) #0
hsC4.sig <- hsC4 %>% filter(padj<0.01 & abs(coef) >= 0.4) #0
hsC5.sig <- hsC5 %>% filter(padj<0.01 & abs(coef) >= 0.4) #110
hsC6.sig <- hsC6 %>% filter(padj<0.01 & abs(coef) >= 0.4) #5
hsC7.sig <- hsC7 %>% filter(padj<0.01 & abs(coef) >= 0.4) #1

#write.csv(hsC1.sig, file="../tables/hsC1_corr_sig.csv", row.names=FALSE)
#write.csv(hsC2.sig, file="../tables/hsC2_corr_sig.csv", row.names=FALSE)

hsc1.sig <- read.csv("tables/hsC1_corr_sig.csv")
hsc2.sig <- read.csv("tables/hsC2_corr_sig.csv")

chipexo <- read.csv("tables/kznfs_TEs_ChIP_exo_modified.csv")
chipexo <- chipexo %>%
    mutate(pair = paste0(teName, ":", geneName))

hsc1.sig <- hsc1.sig %>%
    filter(pair %in% chipexo$pair) %>%
    mutate(group = ifelse(coef>0, "positive", "negative")) #total 869

hsc2.sig <- hsc2.sig %>%
    filter(pair %in% chipexo$pair) %>%
    mutate(group = ifelse(coef>0, "positive", "negative")) # total 399

hsc1.sig.p <- hsc1.sig %>% filter(group=="positive")
hsc1.sig.n <- hsc1.sig %>% filter(group=="negative")
hsc2.sig.p <- hsc2.sig %>% filter(group=="positive")
hsc2.sig.n <- hsc2.sig %>% filter(group=="negative")

c1c2 <- list(
    c1_p = hsc1.sig.p$pair,
    c1_n = hsc1.sig.n$pair,
    c2_p = hsc2.sig.p$pair,
    c2_n = hsc2.sig.n$pair
)

c1c2_m <- make_comb_mat(c1c2)

png("figures/JPG/3D_upset_plot_c1c2.png", width=5, height=3, units="in", res=400)
#svg("figures/SVG/3D_upset_plot_c1c2.svg", width=5, height=3)
u_c1c2 <- UpSet(c1c2_m, comb_order=order(-comb_size(c1c2_m)),
                top_annotation=upset_top_annotation(c1c2_m, add_numbers=TRUE))
print(u_c1c2)
dev.off()
