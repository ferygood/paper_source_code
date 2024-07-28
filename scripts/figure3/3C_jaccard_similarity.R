# fig 3C
# Jaccard Index Statistical test comparing with ChIP-Exo

library(dplyr)
library(twice)
library(ggplot2)
library(ggpubr)

data("hmKZNFs337")
data("hg19rmsk_info")

# chipexo
chipexo <- read.csv("data/kznfs_TEs_ChIP_exo_modified.csv")
chipexo <- chipexo %>%
    mutate(pair=paste(teName, ":", geneName))

# cluster 1
hmc1 <- read.csv("data/hsC1_corr_sig.csv")
hmc1 <- hmc1 %>%
    mutate(pair=paste(teName, ":", geneName))

hmc2 <- read.csv("data/hsC2_corr_sig.csv")
hmc2 <- hmc2 %>%
    mutate(pair=paste(teName, ":", geneName))


## prepare a full combination list


tekrabznf_combination <- c()

for (i in hmKZNFs337$external_gene_name){
    for (j in hg19rmsk_info$gene_id){
        pair <- paste(j, ":", i)
        tekrabznf_combination <- c(tekrabznf_combination, pair)
    }
}


## write a jaccard similarity function

jaccard_similarity <- function(set1, set2) {

    inter_count <- length(intersect(set1, set2))
    union_count <- length(union(set1, set2))

    score <- inter_count / union_count

    score
}

num_selection <- 1000
num_items <- 127797
score_list_c1 <- c()

for (i in 1:num_selection) {
    # Randomly select 127797 items
    selection <- sample(tekrabznf_combination, num_items, replace = TRUE)
    # Store the selection
    jaccard_score <- jaccard_similarity(selection, chipexo$pair)
    score_list_c1 <- c(score_list_c1, jaccard_score)
}


df_c1_simulation <- data.frame(
    simulation = c(replicate(1000, "Primary & Secondary Cortices")),
    Jaccard_similarity = score_list_c1
)


## for cluster 2

num_selection <- 1000
num_items <- 49770
score_list_c2 <- c()

for (i in 1:num_selection) {
    # Randomly select 49110 items
    selection <- sample(tekrabznf_combination, num_items, replace = TRUE)
    # Store the selection
    jaccard_score <- jaccard_similarity(selection, chipexo$pair)
    score_list_c2 <- c(score_list_c2, jaccard_score)
}

df_c2_simulation <- data.frame(
    simulation = c(replicate(1000, "Limbic & Association Cortices")),
    Jaccard_similarity = score_list_c2
)

# combine cluster1 and cluster2
df_combine <- rbind(df_c1_simulation, df_c2_simulation)
df_combine$simulation <- factor(df_combine$simulation,
                                levels = c("Primary & Secondary Cortices",
                                           "Limbic & Association Cortices"))
write.csv(df_combine, file="data/jaccard_simu_table.csv", row.names = F)

g_statistic <- ggplot(df_combine, aes(x=simulation, y=Jaccard_similarity, fill=simulation)) +
    geom_boxplot(alpha=0.2) +
    geom_point(aes(x="Primary & Secondary Cortices", y=jaccard_similarity(hmc1$pair, chipexo$pair)),
               shape=21, fill="#e69138", colour="black", size=3) +
    geom_point(aes(x="Limbic & Association Cortices", y=jaccard_similarity(hmc2$pair, chipexo$pair)),
               shape=21, fill="#fecd3b", colour="black", size=3) +
    annotate(geom="text", x=1.3, y=0.0067, label="italic('p') < 0.001", parse=TRUE, size=5) +
    annotate(geom="text", x=2.3, y=0.0077, label="italic('p') < 0.001", parse=TRUE, size=5) +
    xlab("") +
    ylab("Jaccard Similarity") +
    theme_bw() +
    theme(
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.position = "none"
    ) +
    rotate_x_text(angle=20)

ggsave(g_statistic, file="figures/JPG/3C_simulation_overlap_randomkznf_chipexo_boxplot_c1c2.jpg", width=5, height=5)

