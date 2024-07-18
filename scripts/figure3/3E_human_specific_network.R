# fig 3E
# visualize human network in primary and secondary cortices & limbic
# and association cortices

library(dplyr)
library(RCy3)
library(netZooR)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


hm_c1_sig <- read.csv("tables/hmc1_sig_forNetwork.csv")

c1.node <- data.frame(
    id=c(unique(hm_c1_sig$geneName), unique(hm_c1_sig$teName))
)

# add age information
c1.link <- hm_c1_sig[,c(1,2,5,6)]
colnames(c1.link) <- c("source", "target", "coefficient", "age")

# use condor to calculate the bipartite modularity
node_pair <- hm_c1_sig[,c(2,1,3)]
node_pair$coef <- abs(node_pair$coef)
c1_condor_obj <- createCondorObject(node_pair)
c1_condor_obj <- condorCluster(c1_condor_obj, cs.method = "LEC", project=F)

kznf_com <- c1_condor_obj$blue.memb
colnames(kznf_com)[1] <- "id"
te_com <- c1_condor_obj$red.memb
colnames(te_com)[1] <- "id"
df_com <- rbind(kznf_com, te_com)

c1.node <- c1.node %>%
    left_join(df_com, join_by(id==id))
colnames(c1.node)[2] <- "group"
c1.node$group <- as.character(c1.node$group)

createNetworkFromDataFrames(c1.node, c1.link)

dfc1_condor <- c1.link %>%
    left_join(c1.node, join_by("source"=="id")) %>%
    left_join(c1.node, join_by("target"=="id"))
# divide condor 1
dfc1_condor_g1 <- dfc1_condor %>% filter(group.x=="1" & group.y=="1") #94
c1.g1.node <- data.frame(
    id=c(unique(dfc1_condor_g1$source), unique(dfc1_condor_g1$target))
)

c1.g1.link <- dfc1_condor_g1[,c(1:4)]
createNetworkFromDataFrames(c1.g1.node, c1.g1.link)

# divide condor 3 and 5
create_condor_network <- function(df, n_condor){
    df_subcondor <- df %>% filter(group.x==n_condor & group.y==n_condor) #94
    node <- data.frame(
        id=c(unique(df_subcondor$source), unique(df_subcondor$target))
    )

    link <- df_subcondor[,c(1:4)]
    createNetworkFromDataFrames(node, link)
}

create_condor_network(dfc1_condor, 1)
create_condor_network(dfc1_condor, 2)
create_condor_network(dfc1_condor, 3)
create_condor_network(dfc1_condor, 4)
create_condor_network(dfc1_condor, 5)

# we can do sth like this to calculate the degree
dfc1_condor %>% group_by(source, coefficient, age) %>%
    summarise(count=n()) %>%
    arrange(desc(count)) %>%
    filter(coefficient=="negative" & age=="young")

dfc1_condor_preprocess <- dfc1_condor %>%
    mutate(link_class = paste0(coefficient, "_", age)) %>%
    group_by(group.x, group.y, link_class) %>%
    summarise(count = n())

calculate_bipartite_connectivity <- function(adjacency_matrix, threshold_percent) {
    source_connectivity <- rowSums(adjacency_matrix)
    target_connectivity <- colSums(adjacency_matrix)

    # Create dataframe with node type (source or target), connectivity, and node name
    source_df <- data.frame(node_type = "KRAB-ZNF",
                            connectivity = source_connectivity,
                            node_name = rownames(adjacency_matrix))
    target_df <- data.frame(node_type = "TE",
                            connectivity = target_connectivity,
                            node_name = colnames(adjacency_matrix))

    # Combine dataframes
    connectivity_df <- rbind(source_df, target_df)

    # Calculate hubs based on threshold
    connectivity_df$is_hub <- FALSE
    hub_threshold <- quantile(connectivity_df$connectivity, 1 - threshold_percent/100)
    connectivity_df$is_hub <- ifelse(connectivity_df$connectivity >= hub_threshold, TRUE, FALSE)

    # Calculate gatekeepers
    connectivity_df$is_gatekeeper <- FALSE
    for (i in 1:nrow(connectivity_df)) {
        if (connectivity_df$node_type[i] == "KRAB-ZNF") {
            source_node <- connectivity_df$node_name[i]
            if (sum(adjacency_matrix[source_node,]) == 1) {
                connectivity_df$is_gatekeeper[i] <- TRUE
            }
        } else if (connectivity_df$node_type[i] == "TE") {
            target_node <- connectivity_df$node_name[i]
            if (sum(adjacency_matrix[,target_node]) == 1) {
                connectivity_df$is_gatekeeper[i] <- TRUE
            }
        }
    }

    rownames(connectivity_df) <- 1:nrow(connectivity_df)

    connectivity_df
}

# cluster 1
adjmatrix_c1 <- table(dfc1_condor[,c(1,2)])
c1_connect <- calculate_bipartite_connectivity(adjmatrix_c1, 5)

c1_connect_count <- c1_connect %>%
    group_by(node_type, connectivity) %>%
    summarise(count=n())

# because it is bipartite network, we need to normalized by KRAB-ZNF count and TE count
c1_k_num <- unique(hm_c1_sig$geneName) #101
c1_t_num <- unique(hm_c1_sig$teName) #253

c1_connect_count_normalize <- c1_connect_count %>%
    mutate(norm_degree = ifelse(node_type=="KRAB-ZNF", connectivity*1000/253, connectivity*1000/101))

# Perform multivariate regression
multivar_model <- lm(cbind(norm_degree, count) ~ node_type, data = c1_connect_count_normalize)
summary(multivar_model) #0.0002916 significant

g_c1_connect <- ggplot(c1_connect_count_normalize, aes(x=count, y=norm_degree, fill=node_type)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x = "Number of nodes", y="Normalized Degree") +
    scale_fill_manual(values = c("KRAB-ZNF"="#C994C7", "TE"="#7BCCC4")) +
    theme_bw()

ggsave(g_c1_connect, file="figures/JPG/3F_node_connectivity_c1.jpg",
       dpi=400,width=4,height=2.5)
