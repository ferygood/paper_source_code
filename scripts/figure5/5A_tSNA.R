# 5A
# tSNE variation in Figure 5.

library(twice)
library(dplyr)
library(ggpubr)
library(Rtsne)
data("hg38rmsk_info")
data("hmKZNFs337")
load("data/mayoTEKRABber_balance.RData")
mayo_meta <- read.csv("data/selectSample.csv")

cbeCtrlCorr <- mayoTEKRABber$cbeControlCorr
cbeADCorr <- mayoTEKRABber$cbeADCorr
cbeDE <- mayoTEKRABber$cbeDE

tcxCtrlCorr <- mayoTEKRABber$tcxControlCorr
tcxADCorr <- mayoTEKRABber$tcxADCorr
tcxDE <- mayoTEKRABber$tcxDE

# 1. combine gene and te table in cbe
# kznfs
cbeGene <- data.frame(cbeDE$normalized_gene_counts)
cbeKZNFs <- cbeGene %>% filter(rownames(cbeGene) %in% hmKZNFs337$ensembl_gene_id)
dfCbeKZNFs <- data.frame(t(cbeKZNFs))
rownames(dfCbeKZNFs) <- substr(rownames(dfCbeKZNFs), 2, 10000000)

# TEs
cbeTE <- data.frame(cbeDE$normalized_te_counts)
dfCbeTEs <- data.frame(t(cbeTE))
rownames(dfCbeTEs) <- substr(rownames(dfCbeTEs), 2, 10000000)

# combine column
dfCbe <- cbind(dfCbeKZNFs, dfCbeTEs)


# 2. combine gene and te table in tcx
# kznfs
tcxGene <- data.frame(tcxDE$normalized_gene_counts)
tcxKZNFs <- tcxGene %>% filter(rownames(tcxGene) %in% hmKZNFs337$ensembl_gene_id)
dfTcxKZNFs <- data.frame(t(tcxKZNFs))
rownames(dfTcxKZNFs) <- substr(rownames(dfTcxKZNFs), 2, 100000000)

# TEs
tcxTE <- data.frame(tcxDE$normalized_te_counts)
dfTcxTEs <- data.frame(t(tcxTE))
rownames(dfTcxTEs) <- substr(rownames(dfTcxTEs), 2, 100000000)

# conbine colum
dfTcx <- cbind(dfTcxKZNFs, dfTcxTEs)

# 3. combine cbe and tcx
df_mayo <- bind_rows(dfCbe, dfTcx)

# 4. add information of diagnosis
df_mayo$name <- rownames(df_mayo)
df_mayo_input <- df_mayo %>%
    left_join(mayo_meta, join_by(name==name))
    #select(c(1:1048))
rownames(df_mayo_input) <- rownames(df_mayo)

set.seed(42)

df_mayo_noNA <- df_mayo_input %>% select_if(~ !any(is.na(.)))

tsne <- Rtsne(
    as.matrix(df_mayo_noNA[,1:1048]),
    pca=TRUE,
    perplexity = 15,
    dims = 2
)

# create results dataframe
tsne_index <- data.frame(tsne$Y)
rownames(tsne_index) <- rownames(df_mayo_input)
tsne_index$diagnosis <- df_mayo_input$diagnosis
colnames(tsne_index)[c(1,2)] <- c("t-SNE1", "t-SNE2")

g1 <- ggplot(tsne_index, aes(x=`t-SNE1`, y=`t-SNE2`)) +
    geom_point(colour="black", shape=21, size=3,
               aes(fill=factor(diagnosis))) +
    scale_fill_manual(values=c("#f96161", "#d0b783", "#99CC99", "#CC99FF")) +
    labs(fill = "diagnosis") +
    #geom_text(label=rownames(tsne_index), size=1.5) +
    theme_bw()
ggsave(g1, filename="figures/JPG/5A_kznfsTEsTSNE.jpg", dpi=400,
       width=5, height=3, bg="white")
