# fig 3A and figure S5
# check the 1000 iterations of KRAB-ZNFs (with random gene sets) to TEs

library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(foreach)
library(doParallel)


get_count_parallel <- function(dir_path){
    # Create a vector of file names
    file_names <- paste0(dir_path, 1:1000, "_vs_TE_corr.csv")

    # Initialize parallel backend
    registerDoParallel(5)

    # Use foreach loop for parallel processing
    df_temp <- foreach(file_name = file_names, .combine = rbind) %dopar% {
        data <- read.csv(file_name)

        filtered_data <- data %>% filter(padj < 0.01)
        p_count <- filtered_data %>% filter(coef > 0) %>% nrow()
        n_count <- filtered_data %>% filter(coef < 0) %>% nrow()

        df_add <- data.frame(
            positive_count = p_count,
            negative_count = n_count
        )

        return(df_add)
    }

    # Stop parallel backend
    stopImplicitCluster()

    return(df_temp)
}

count_hm_c1 <- get_count_parallel("/group/ag_nowick/data/yao/results_c1/gene_")
count_hm_c2 <- get_count_parallel("/group/ag_nowick/data/yao/results_c2/gene_")

count_pt_c1 <- get_count_parallel("/group/ag_nowick/data/yao/results_c1_pt/gene_")
count_pt_c2 <- get_count_parallel("/group/ag_nowick/data/yao/results_c2_pt/gene_")

count_hm_cbe <- get_count_parallel("/group/ag_nowick/data/yao/results_cbe/gene_")
count_hm_tcx <- get_count_parallel("/group/ag_nowick/data/yao/results_tcx/gene_")

count_pp_c1 <- get_count_parallel("/group/ag_nowick/data/yao/results_c1_pp/gene_")
count_pp_c2 <- get_count_parallel("/group/ag_nowick/data/yao/results_c2_pp/gene_")

count_mm_c1 <- get_count_parallel("/group/ag_nowick/data/yao/results_c1_mm/gene_")
count_mm_c2 <- get_count_parallel("/group/ag_nowick/data/yao/results_c2_mm/gene_")

count_hm_c1_p <- count_hm_c1 %>% mutate(label="hm_cluster1")
count_hm_c2_p <- count_hm_c2 %>% mutate(label="hm_cluster2")
count_pt_c1_p <- count_pt_c1 %>% mutate(label="pt_cluster1")
count_pt_c2_p <- count_pt_c2 %>% mutate(label="pt_cluster2")
count_pp_c1_p <- count_pp_c1 %>% mutate(label="pp_cluster1")
count_pp_c2_p <- count_pp_c2 %>% mutate(label="pp_cluster2")
count_mm_c1_p <- count_mm_c1 %>% mutate(label="mm_cluster1")
count_mm_c2_p <- count_mm_c2 %>% mutate(label="mm_cluster2")
count_hm_cbe_p <- count_hm_cbe %>% mutate(label="cbe_control")
count_hm_tcx_p <- count_hm_tcx %>% mutate(label="tcx_control")

df_combine <- rbind(count_hm_c1_p, count_hm_c2_p, count_hm_cbe_p, count_hm_tcx_p,
                    count_mm_c1_p, count_mm_c2_p, count_pt_c1_p, count_pt_c2_p,
                    count_pp_c1_p, count_pp_c2_p)

write.csv(df_combine, file="tables/iterationAll_tables.csv", row.names=F)
df_combine <- read.csv("tables/iterationAll_tables.csv")

# visualize cluster 1 figure
df_cluster1_kznf <- data.frame(
    count = c(127797, 100987, 26810),
    group = c("all", "positive", "negative")
)

df_cluster1 <- df_combine %>%
    filter(label=="hm_cluster1") %>%
    select(c(1,2)) %>%
    mutate(all = positive_count + negative_count)

colnames(df_cluster1) <- c("positive", "negative", "all")
df_cluster1 <- pivot_longer(df_cluster1,
                            cols = c(all, positive, negative),
                            names_to = "group",
                            values_to = "count")

df_cluster1$label <- factor(df_cluster1$group,
                            levels=c("all", "positive", "negative"))

hm_cluster1 <- ggplot(data = df_cluster1, aes(x=group, y=log(count), fill=group)) +
    geom_boxplot() +
    geom_point(data = df_cluster1_kznf,
               aes(x=group, y=log(count)),
               color="#c27ba0",
               size=2) +
    scale_fill_manual(values = c("all" = "#999999", "positive" = "#cc3433", "negative" = "#336699")) +
    xlab("") +
    ylab("Normalized correlation counts") +
    ggtitle("Primary and Secondary Cortices (Cluster 1)") +
    theme_bw() +
    theme(plot.title=element_text(size=8),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8))

ggsave(hm_cluster1, file="figures/JPG/3A_correlation_check_human_cluster1.jpg",
       dpi=400, width=4, height=4)

# visualize cluster 2
df_cluster2_kznf <- data.frame(
    count = c(49770, 38295, 11475),
    group = c("all", "positive", "negative")
)

df_cluster2 <- df_combine %>%
    filter(label=="hm_cluster2") %>%
    select(c(1,2)) %>%
    mutate(all = positive_count + negative_count)

colnames(df_cluster2) <- c("positive", "negative", "all")
df_cluster2 <- pivot_longer(df_cluster2,
                            cols = c(all, positive, negative),
                            names_to = "group",
                            values_to = "count")

df_cluster2$group <- factor(df_cluster2$group,
                            levels=c("all", "positive", "negative"))

hm_cluster2 <- ggplot(data = df_cluster2,
                      aes(x=group, y=log(count), fill=group)) +
    geom_boxplot() +
    geom_point(data = df_cluster2_kznf,
               aes(x=group, y=log(count)),
               color="#c27ba0",
               size=2) +
    scale_fill_manual(values = c("all" = "#999999", "positive" = "#cc3433", "negative" = "#336699")) +
    xlab("") +
    ylab("Normalized correlation counts") +
    ggtitle("Limbic and Association Cortices(Cluster 2)") +
    theme_bw() +
    theme(plot.title=element_text(size=8),
          axis.title.x=element_text(size=8),
          axis.title.y=element_text(size=8))

ggsave(hm_cluster2, file="figures/JPG/S3A_correlation_check_human_cluster2.jpg",
       dpi=400, width=4, height=4)

# create comparison of human and NHPs
df_combine_c1 <- df_combine %>%
    filter(!label %in% c("cbe_control", "tcx_control")) %>%
    filter(label %in%
               c("hm_cluster1", "pt_cluster1", "pp_cluster1", "mm_cluster1")) %>%
    mutate(all = positive_count + negative_count)

colnames(df_combine_c1)[c(1,2)] <- c("positive", "negative")

df_c1_pivot <- df_combine_c1 %>%
    pivot_longer(cols = c(all, positive, negative),
                 names_to = "group",
                 values_to = "count")
df_c1_pivot$label <- factor(df_c1_pivot$label,
                            levels = c("hm_cluster1", "pt_cluster1",
                                       "pp_cluster1", "mm_cluster1"))
df_c1_pivot$group <- factor(df_c1_pivot$group, levels=c("all","positive","negative"))

# TE:KRAB-ZNF
df_c1_tekrabznf <- data.frame(
    label = c("hm_cluster1", "hm_cluster1", "hm_cluster1",
              "pt_cluster1", "pt_cluster1", "pt_cluster1",
              "pp_cluster1", "pp_cluster1", "pp_cluster1",
              "mm_cluster1", "mm_cluster1", "mm_cluster1"),
    group = c("all", "positive","negative",
              "all", "positive","negative",
              "all", "positive", "negative",
              "all", "positive", "negative"),
    count = c(127797, 100987, 26810,
              5793, 5741, 52,
              46759, 25293, 21466,
              2273, 2052, 221)
)

nhp_c1_g <- ggplot(df_c1_pivot, aes(x=label, y=log(count))) +
    geom_boxplot() +
    geom_point(data = df_c1_tekrabznf,
               aes(x=label, y=log(count)),
               color="#c27ba0",
               size=2) +
    facet_wrap(~group) +
    ylab("normalized correlation counts") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(nhp_c1_g, file="figures/JPG/S5A_nhp_c1_boxplot.jpg",
       dpi=400, width=5, height=4)

#compare cluster 2 with human and NHPs
df_combine_c2 <- df_combine %>%
    filter(!label %in% c("cbe_control", "tcx_control")) %>%
    filter(label %in%
               c("hm_cluster2", "pt_cluster2", "pp_cluster2", "mm_cluster2")) %>%
    mutate(all = positive_count + negative_count)

colnames(df_combine_c2)[c(1,2)] <- c("positive", "negative")

df_c2_pivot <- df_combine_c2 %>%
    pivot_longer(cols = c(all, positive, negative),
                 names_to = "group",
                 values_to = "count")
df_c2_pivot$label <- factor(df_c2_pivot$label,
                            levels = c("hm_cluster2", "pt_cluster2",
                                       "pp_cluster2", "mm_cluster2"))
df_c2_pivot$group <- factor(df_c2_pivot$group, levels=c("all","positive","negative"))

# TE:KRAB-ZNF
df_c2_tekrabznf <- data.frame(
    label = c("hm_cluster2", "hm_cluster2", "hm_cluster2",
              "pt_cluster2", "pt_cluster2", "pt_cluster2",
              "pp_cluster2", "pp_cluster2", "pp_cluster2",
              "mm_cluster2", "mm_cluster2", "mm_cluster2"),
    group = c("all", "positive","negative",
              "all", "positive", "negative",
              "all", "positive", "negative",
              "all", "positive", "negative"),
    count = c(49770, 38295, 11475,
              6185, 6018, 167,
              20070, 19254, 816,
              6253, 5305, 948)
)

nhp_c2_g <- ggplot(df_c2_pivot, aes(x=label, y=log(count))) +
    geom_boxplot() +
    geom_point(data = df_c2_tekrabznf,
               aes(x=label, y=log(count)),
               color="#c27ba0",
               size=2) +
    facet_wrap(~group) +
    ylab("normalized correlation counts") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(nhp_c2_g, file="figures/JPG/S5B_nhp_c2_boxplot.jpg",
       dpi=400, width=5, height=4)


