load("data/primateBrainData.RData")

brain_meta <- data.frame(
    region = sort(unique(metadata$brain_region))
)

brain_meta$cluster <- c("cluster1", "cluster1", rep("cluster2", 8),
                        "cluster3", "cluster1", rep("cluster3", 2),
                        rep("cluster4", 2), rep("cluster5", 6), "cluster1",
                        "cluster6", rep("cluster7", 3), rep("cluster1", 6)
)

brain_meta$cluster <- factor(brain_meta$cluster)
