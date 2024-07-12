library(Rtsne)
library(tidyverse)
library(twice)
library(ggpubr)
theme_set(theme_bw(18))

data("hmKZNFs337")

# set the path to your expression raw count data
read_rds_files <- function(file_path){

    file_list <- list.files(path = folder_path, pattern="\\.rds$", full.names = TRUE)
    for (file_path in file_list){
        file_name <- basename(file_path)
        var_name <- gsub("\\.rds$", "", file_name)
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
    }
}

folder_path <- "../../paper_source_code/data/results_rdata/"
read_rds_files(folder_path)

# create a metadata for labeling sample
# 1. for human
create_human_dataset <- function(DEobject, region){

    df_kznfs <- DEobject$geneCorrInputRef %>%
        dplyr::filter(rownames(.) %in% hmKZNFs337$external_gene_name)

    df_Hm <- data.frame(cbind(
        t(df_kznfs),
        t(DEobject$teCorrInputRef)
    ))

    df_Hm <- df_Hm %>%
        mutate(species="Human", region=region)
    df_Hm
}

H1 <- create_human_dataset(HmPtC1$DEobject, "Primary & Secondary Cortices")
H2 <- create_human_dataset(HmPtC2$DEobject, "Limbic & Association Cortices")
H3 <- create_human_dataset(HmPtC3$DEobject, "Archicortex")
H4 <- create_human_dataset(HmPtC4$DEobject, "Thalamus & Hypothalamus")
H5 <- create_human_dataset(HmPtC5$DEobject, "Cerebellar White Matter")
H6 <- create_human_dataset(HmPtC6$DEobject, "Cerebellar Grey Matter")
H7 <- create_human_dataset(HmPtC7$DEobject, "Striatum")
df_Hm <- bind_rows(H1, H2, H3, H4, H5, H6, H7)
df_Hm <- df_Hm[,c(1:1167)]

# 2. for chimpanzee
create_chimp_dataset <- function(DEobject, region){

    df_kznfs <- DEobject$geneCorrInputCompare %>%
        dplyr::filter(rownames(.) %in% hmKZNFs337$external_gene_name)

    df_chimp <- data.frame(cbind(
        t(df_kznfs),
        t(DEobject$teCorrInputCompare)
    ))

    df_chimp <- df_chimp %>%
        mutate(species="Chimpanzee", region=region)
    df_chimp
}

C1 <- create_chimp_dataset(HmPtC1$DEobject, "Primary & Secondary Cortices")
C2 <- create_chimp_dataset(HmPtC2$DEobject, "Limbic & Association Cortices")
C3 <- create_chimp_dataset(HmPtC3$DEobject, "Archicortex")
C4 <- create_chimp_dataset(HmPtC4$DEobject, "Thalamus & Hypothalamus")
C5 <- create_chimp_dataset(HmPtC5$DEobject, "Cerebellar White Matter")
C6 <- create_chimp_dataset(HmPtC6$DEobject, "Cerebellar Grey Matter")
C7 <- create_chimp_dataset(HmPtC7$DEobject, "Striatum")
df_Chimp <- bind_rows(C1, C2, C3, C4, C5, C6, C7)
df_Chimp <- df_Chimp[,c(1:1167)]

# 3. for bonobo
create_bonobo_dataset <- function(DEobject, region){

    df_kznfs <- DEobject$geneCorrInputCompare %>%
        dplyr::filter(rownames(.) %in% hmKZNFs337$external_gene_name)

    df_bonobo <- data.frame(cbind(
        t(df_kznfs),
        t(DEobject$teCorrInputCompare)
    ))

    df_bonobo <- df_bonobo %>%
        mutate(species="Bonobo", region=region)
    df_bonobo
}

B1 <- create_bonobo_dataset(HmPpC1$DEobject, "Primary & Secondary Cortices")
B2 <- create_bonobo_dataset(HmPpC2$DEobject, "Limbic & Association Cortices")
B3 <- create_bonobo_dataset(HmPpC3$DEobject, "Archicortex")
B4 <- create_bonobo_dataset(HmPpC4$DEobject, "Thalamus & Hypothalamus")
B5 <- create_bonobo_dataset(HmPpC5$DEobject, "Cerebellar White Matter")
B6 <- create_bonobo_dataset(HmPpC6$DEobject, "Cerebellar Grey Matter")
B7 <- create_bonobo_dataset(HmPpC7$DEobject, "Striatum")
df_Bonobo <- bind_rows(B1, B2, B3, B4, B5, B6, B7)

# 4. for macaque
create_mm_dataset <- function(DEobject, region){

    df_kznfs <- DEobject$geneCorrInputCompare %>%
        dplyr::filter(rownames(.) %in% hmKZNFs337$external_gene_name)

    df_mm <- data.frame(cbind(
        t(df_kznfs),
        t(DEobject$teCorrInputCompare)
    ))

    df_mm <- df_mm %>%
        mutate(species="Macaque", region=region)
    df_mm
}

M1 <- create_mm_dataset(HmMmC1$HmMmC1_DE, "Primary & Secondary Cortices")
M2 <- create_mm_dataset(HmMmC2$HmMmC2_DE, "Limbic & Association Cortices")
M3 <- create_mm_dataset(HmMmC3$HmMmC3_DE, "Archicortex")
M4 <- create_mm_dataset(HmMmC4$HmMmC4_DE, "Thalamus & Hypothalamus")
M5 <- create_mm_dataset(HmMmC5$HmMmC5_DE, "Cerebellar White Matter")
M6 <- create_mm_dataset(HmMmC6$HmMmC6_DE, "Cerebellar Grey Matter")
M7 <- create_mm_dataset(HmMmC7$HmMmC7_DE, "Striatum")
df_Mm <- bind_rows(M1, M2, M3, M4, M5, M6, M7)

# combine all expression data
df_input <- bind_rows(df_Hm, df_Chimp, df_Bonobo, df_Mm)
df_input <- df_input %>% replace(is.na(.), 0)
df_input_reorder <- df_input[,c(1:1165, 1168:1196, 1166, 1167)] # put species and region at the end

# running tSNE
set.seed(42)
input_matrix <- df_input_reorder %>%
    select(-species, -region) %>%
    as.matrix()

tsne <- Rtsne(
    input_matrix,
    pca=TRUE,
    perplexity = floor((nrow(df_input_reorder) -1)/ 3),
    dims = 2
)

tsne_index <- as.data.frame(tsne$Y)
rownames(tsne_index) <- rownames(df_input_reorder)
tsne_index$brain_region <- df_input_reorder$region
tsne_index$species <- df_input_reorder$species
colnames(tsne_index)[c(1,2)] <- c("t-SNE1", "t-SNE2")

# convert column type into factor
tsne_index$species <- factor(
    tsne_index$species,
    levels = c("Human", "Chimpanzee", "Bonobo", "Macaque"))

tsne_index$brain_region <- factor(tsne_index$brain_region,
                                  levels=c("Primary & Secondary Cortices",
                                           "Limbic & Association Cortices",
                                           "Archicortex",
                                           "Thalamus & Hypothalamus",
                                           "Cerebellar White Matter",
                                           "Cerebellar Grey Matter",
                                           "Striatum"))

# draw tSNE plot and save ggplot
g1 <- ggplot(tsne_index, aes(x=`t-SNE1`, y=`t-SNE2`)) +
    geom_point(colour="black", shape=21, size=3,
               aes(fill=factor(species))) +
    #scale_fill_manual(values=c("#F0997B", "#6F92E9", "#AD73F6", "#8BE57D")) +
    scale_fill_manual(values=c("red", "blue", "purple", "darkgreen")) +
    labs(fill = "")

ggsave(g1, filename="figures/SVG/2A_PBD_tSNE_species.svg",
       dpi=400, width=7, height=5, bg="white")
ggsave(g1, filename="figures/JPG/2A_PBD_tSNE_species.jpg",
       dpi=400, width=7, height=5, bg="white")


g2 <- ggplot(tsne_index, aes(x=`t-SNE1`, y=`t-SNE2`)) +
    geom_point(colour="black", shape=21, size=3,
               aes(fill=factor(brain_region))) +
    scale_fill_manual(values=c("#ff9a85", "#e3c678", "#b07ff3", "#afdb72", "#656662", "#5adfd7", "#3d70ca")) +
    labs(fill = "")

ggsave(g2, filename="figures/SVG/2A_PBD_tSNE_region.svg",
       dpi=400, width=7, height=5, bg="white")
ggsave(g2, filename="figures/JPG/2A_PBD_tSNE_region.jpg",
       dpi=400, width=7, height=5, bg="white")

