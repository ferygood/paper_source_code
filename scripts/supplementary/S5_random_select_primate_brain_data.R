#S5
#S4
# Create network using highly confident TE:KRAB-ZNF

In this script, we draw the cluster 1 and cluster 2 network in human. We also highlight the young pair of TE:KRAB-ZNF.

![](/figures/hmc1_869_age_network.png){width="649"}

![](/figures/networkImageC1/node_connectivity_barplot.jpg){width="597"}

Multivariate regression test (p\<0.001, node_typeTE) by combing the normalized degree and number of nodes.

![](/figures/networkImageC2/node_connectivity_barplot.jpg){width="602"}

```{r}
knitr::opts_chunk$set(echo=FALSE, message=FALSE, eval=FALSE)
```

```{r eval=FALSE, message=FALSE}
library(dplyr)
library(RCy3)
library(netZooR)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
```

Load data

In cluster 1

```{r eval=FALSE}
hm_c1_sig <- read.csv("../tables/hmc1_sig_forNetwork.csv")

c1.node <- data.frame(
    id=c(unique(hm_c1_sig$geneName), unique(hm_c1_sig$teName))
)

# add age information
c1.link <- hm_c1_sig[,c(1,2,5,6)]
colnames(c1.link) <- c("source", "target", "coefficient", "age")
```

```{r eval=FALSE}
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
```

```{r}
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
```

```{r}
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
```

```{r}
# we can do sth like this to calculate the degree
dfc1_condor %>% group_by(source, coefficient, age) %>%
    summarise(count=n()) %>%
    arrange(desc(count)) %>%
    filter(coefficient=="negative" & age=="young")

dfc1_condor_preprocess <- dfc1_condor %>%
    mutate(link_class = paste0(coefficient, "_", age)) %>%
    group_by(group.x, group.y, link_class) %>%
    summarise(count = n())



```

Calculate the properties of network

Create the adjacency matrix

```{r}
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
```

```{r}
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

ggsave(g_c1_connect, file="../figures/networkImageC1/node_connectivity_barplot.jpg",
       dpi=400,width=4,height=2.5)


```

## cluster 2

```{r}
hm_c2_sig <- read.csv("../tables/hmc2_sig_forNetwork.csv")

c2.node <- data.frame(
    id=c(unique(hm_c2_sig$geneName), unique(hm_c2_sig$teName))
)

# add age information
c2.link <- hm_c2_sig[,c(1,2,5,6)]
colnames(c2.link) <- c("source", "target", "coefficient", "age")

# use condor to calculate the bipartite modularity
node_pair <- hm_c2_sig[,c(2,1,3)]
node_pair$coef <- abs(node_pair$coef)
c2_condor_obj <- createCondorObject(node_pair)
c2_condor_obj <- condorCluster(c2_condor_obj, cs.method = "LEC", project=F)

kznf_com <- c2_condor_obj$blue.memb
colnames(kznf_com)[1] <- "id"
te_com <- c2_condor_obj$red.memb
colnames(te_com)[1] <- "id"
df_com <- rbind(kznf_com, te_com)

c2.node <- c2.node %>%
    left_join(df_com, join_by(id==id))
colnames(c2.node)[2] <- "group"
c2.node$group <- as.character(c2.node$group)

createNetworkFromDataFrames(c2.node, c2.link)

```

```{r}
dfc2_condor <- c2.link %>%
    left_join(c2.node, join_by("source"=="id")) %>%
    left_join(c2.node, join_by("target"=="id"))

# example using create_condor_network function to extract subnetwork
create_condor_network(dfc2_condor, 2)
create_condor_network(dfc2_condor, 3)
create_condor_network(dfc2_condor, 5)
create_condor_network(dfc2_condor, 12)
```

```{r}
dfc2_condor_preprocess <- dfc2_condor %>%
    mutate(link_class = paste0(coefficient, "_", age)) %>%
    group_by(group.x, group.y, link_class) %>%
    summarise(count = n())

# cluster 2
adjmatrix_c2 <- table(dfc2_condor[,c(1,2)])
c2_connect <- calculate_bipartite_connectivity(adjmatrix_c2, 5)

c2_connect_count <- c2_connect %>%
    group_by(node_type, connectivity) %>%
    summarise(count=n())

# because it is bipartite network, we need to normalized by KRAB-ZNF count and TE count
c2_k_num <- unique(hm_c2_sig$geneName) #66
c2_t_num <- unique(hm_c2_sig$teName) #181

c2_connect_count_normalize <- c2_connect_count %>%
    mutate(norm_degree = ifelse(node_type=="KRAB-ZNF", connectivity*1000/66, connectivity*1000/181))

# Perform multivariate regression
multivar_model <- lm(cbind(norm_degree, count) ~ node_type, data = c2_connect_count_normalize)
summary(multivar_model) #0.00461 significant

g_c2_connect <- ggplot(c2_connect_count_normalize, aes(x=count, y=norm_degree, fill=node_type)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x = "Number of nodes", y="Normalized Degree") +
    scale_fill_manual(values = c("KRAB-ZNF"="#C994C7", "TE"="#7BCCC4")) +
    theme_bw()

ggsave(g_c2_connect, file="../figures/networkImageC2/node_connectivity_barplot.jpg",
       dpi=400,width=4,height=2.5)


```

```{r}
# create a condor rds file with the results
library(twice)
data("hg19rmsk_info")
te_annot <- hg19rmsk_info
te_annot$family_id[te_annot$class_id=="SVA"] <- "SVA"

c1_condor <- dfc1_condor %>%
    left_join(te_annot[,c(1,2)], join_by(target==gene_id)) %>%
    mutate(link = paste0(coefficient, "-", age))

c2_condor <- dfc2_condor %>%
    left_join(te_annot[,c(1,2)], join_by(target==gene_id)) %>%
    mutate(link = paste0(coefficient, "-", age))

# save as rds
c1c2_condor_table <- list(
    "c1_condor" = c1_condor,
    "c2_condor" = c2_condor
)

saveRDS(c1c2_condor_table, file="../data/c1c2_condor_table.rds")
condorc1c2 <- readRDS("../data/c1c2_condor_table.rds")
```

## create heatmap

```{r}
c1_condor_group <- condorc1c2$c1_condor %>%
    dplyr::select(c(7,8)) %>%
    group_by(family_id, link) %>%
    summarise(count=n())

c2_condor_group <- condorc1c2$c2_condor %>%
    dplyr::select(c(7,8)) %>%
    group_by(family_id, link) %>%
    summarise(count=n())
```

```{r}
# Assuming you have already loaded the ComplexHeatmap package
modified_c1 <- c1_condor_group %>%
    pivot_wider(names_from = link, values_from = count, values_fill=NA)

matrix_c1 <- as.matrix(modified_c1[,-1])
rownames(matrix_c1) <- modified_c1$family_id
colnames(matrix_c1) <- c("N-O", "N-Y", "P-Y", "P-O")

col_fun = colorRamp2(c(0, 3, 3.5, 6), c("#F7EACE", "#F3CDBC", "#AC8583", "#D599CB"))

gc1 <- Heatmap(log(matrix_c1[,c(4,3,1,2)]), name="log correlations",
               width=ncol(matrix_c1)*unit(6, "mm"),
               height=nrow(matrix_c1)*unit(6, "mm"),
               cluster_columns = F,
               cluster_rows = F,
               rect_gp = gpar(col="black", lwd=0.8),
               na_col="white",
               col = col_fun)

tidyHeatmap::save_pdf(gc1, "../figures/networkImageC1/TEfamilyHeatmap.pdf",
                      width=5, height=5, units=c("in"))

# cluster 2
modified_c2 <- c2_condor_group %>%
    pivot_wider(names_from = link, values_from = count, values_fill=NA)

matrix_c2 <- as.matrix(modified_c2[,-1])
rownames(matrix_c2) <- modified_c2$family_id
colnames(matrix_c2) <- c("N-O", "N-Y", "P-Y", "P-O")

gc2 <- Heatmap(log(matrix_c2[,c(4,3,1,2)]), name="log correlations",
               width=ncol(matrix_c2)*unit(6, "mm"),
               height=nrow(matrix_c2)*unit(6, "mm"),
               cluster_columns = F,
               cluster_rows = F,
               rect_gp = gpar(col="black", lwd=0.8),
               na_col="white",
               col = col_fun)

tidyHeatmap::save_pdf(gc2, "../figures/networkImageC2/TEfamilyHeatmap.pdf",
                      width=5, height=5, units=c("in"))
```

## Use Chi-squared test to see their distribution

We use cluster 1 as an example

```{r}
matrix_c1[is.na(matrix_c1)] <- 0

chisq_test <- function(row_index, matrix_data) {
    row_data <- matrix_data[row_index, ]
    other_data <- matrix_data[-row_index, ]

    # Sum the other rows to get expected counts
    expected_counts <- colSums(other_data)

    # Combine row_data and expected_counts into a matrix
    observed <- rbind(row_data, expected_counts)

    # Perform the Chi-squared test
    test_result <- chisq.test(observed)

    return(test_result$p.value)
}

# Apply the test for each row
chisq_p_values <- sapply(1:nrow(matrix_c1), chisq_test, matrix_data = matrix_c1)

# Combine row names with p-values
chisq_result <- data.frame(
    row_names = rownames(matrix_c1),
    p_value = chisq_p_values
)

# Print the result
print(chisq_result)
```

Create barplot

```{r}
c1_condor_group <- c1_condor_group %>%
    mutate(link = case_when(
        link == "negative-old" ~ "N-O",
        link == "negative-young" ~ "N-Y",
        link == "positive-old" ~ "P-O",
        link == "positive-young" ~ "P-Y",
        TRUE ~ link  # This keeps the original value if no conditions are met
    ))

c1_condor_group$link <- factor(c1_condor_group$link,
                               levels=c("P-O", "P-Y", "N-O", "N-Y"))

gbarc1 <- ggplot(c1_condor_group, aes(x=link, y=count)) +
    geom_bar(stat = "identity") +
    facet_wrap(~family_id, scales = "free_y") +
    theme_bw() +
    ggtitle("cluster 1") +
    theme(axis.text.x = element_text(angle=70, vjust=0.5))

ggsave(gbarc1, file="../figures/networkImageC1/teFamilyFreqC1.svg")
ggsave(gbarc1, file="../figures/networkImageC1/teFamilyFreqC1.png")

c2_condor_group <- c2_condor_group %>%
    mutate(link = case_when(
        link == "negative-old" ~ "N-O",
        link == "negative-young" ~ "N-Y",
        link == "positive-old" ~ "P-O",
        link == "positive-young" ~ "P-Y",
        TRUE ~ link  # This keeps the original value if no conditions are met
    ))

c2_condor_group$link <- factor(c2_condor_group$link,
                               levels=c("P-O", "P-Y", "N-O", "N-Y"))

gbarc2 <- ggplot(c2_condor_group, aes(x=link, y=count)) +
    geom_bar(stat = "identity") +
    facet_wrap(~family_id, scales = "free_y") +
    theme_bw() +
    ggtitle("cluster 2") +
    theme(axis.text.x = element_text(angle=70, vjust=0.5))

ggsave(gbarc2, file="../figures/networkImageC2/teFamilyFreqC2.svg")
ggsave(gbarc2, file="../figures/networkImageC2/teFamilyFreqC2.png")

```

```{r}
# do statistical test
# cluster 1
modified_c1 <- c1_condor_group %>%
    pivot_wider(names_from = link, values_from = count, values_fill=0)

# Initialize an empty dataframe to store the p-values
p_values_dfc1 <- data.frame(family_id = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each row
for (i in 1:nrow(modified_c1)) {
    # Extract the family_id and counts for the current row
    row_family_id <- modified_c1$family_id[i]
    row_counts <- as.numeric(modified_c1[i, -1])

    # Subset the dataframe excluding the current row
    data_to_compare <- modified_c1[-i, ]

    # Combine counts of all other rows
    other_rows_counts <- colSums(data_to_compare[, -1], na.rm = TRUE)

    # Perform chi-square test
    chi_square_result <- chisq.test(rbind(row_counts, other_rows_counts))

    # Extract the p-value
    p_value <- chi_square_result$p.value

    # Add the family_id and p-value to the dataframe
    p_values_dfc1 <- rbind(p_values_dfc1, data.frame(family_id = row_family_id, p_value = p_value, stringsAsFactors = FALSE))
}

# Print the dataframe with p-values
print(p_values_dfc1)

p_values_dfc1 <- p_values_dfc1 %>%
    mutate(sig=ifelse(p_value < 0.001, "significant", "not significant"))

write.csv(p_values_dfc1, file = "../figures/networkImageC1/chisq_c1_teFamily.csv", row.names = F)

```

```{r}
# do statistical test
# cluster 2
modified_c2 <- c2_condor_group %>%
    pivot_wider(names_from = link, values_from = count, values_fill=0)

# Initialize an empty dataframe to store the p-values
p_values_dfc2 <- data.frame(family_id = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each row
for (i in 1:nrow(modified_c2)) {
    # Extract the family_id and counts for the current row
    row_family_id <- modified_c2$family_id[i]
    row_counts <- as.numeric(modified_c2[i, -1])

    # Subset the dataframe excluding the current row
    data_to_compare <- modified_c2[-i, ]

    # Combine counts of all other rows
    other_rows_counts <- colSums(data_to_compare[, -1], na.rm = TRUE)

    # Perform chi-square test
    chi_square_result <- chisq.test(rbind(row_counts, other_rows_counts))

    # Extract the p-value
    p_value <- chi_square_result$p.value

    # Add the family_id and p-value to the dataframe
    p_values_dfc2 <- rbind(p_values_dfc2, data.frame(family_id = row_family_id, p_value = p_value, stringsAsFactors = FALSE))
}

# Print the dataframe with p-values
print(p_values_dfc2)

p_values_dfc2 <- p_values_dfc2 %>%
    mutate(sig=ifelse(p_value < 0.001, "significant", "not significant"))

write.csv(p_values_dfc2, file = "../figures/networkImageC2/chisq_c2_teFamily.csv", row.names = F)
```

Let see how about KRAB-ZNFs

```{r}
znfc1 <- condorc1c2$c1_condor %>%
    group_by(source, link) %>%
    summarise(count=n())

source_count <- table(znfc1$source)
single_count <- names(source_count[source_count==2])


znfc1_filter <- znfc1 %>% filter(source %in% single_count)
```
