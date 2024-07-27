library(dplyr)
library(twice)
library(TEKRABber)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)


data("hmKZNFs337") #337 KRAB-ZNFs in human
data("hg19rmsk_info")

# 1. create meta data
df_meta <- metadata %>%
    left_join(brain_meta, join_by(brain_region == region))


# 2. do TPM conversion in human and NHPs
# 2-1. TPM in human
# genes convert to tpm
hmGene_temp <- hmGene[,c(-1)]
sample_counts <- colSums(hmGene_temp)

scaling_factor <- sample_counts / 1e6

hmGene_tpm <- hmGene_temp
hmGene_tpm <- t(t(hmGene_tpm)/ scaling_factor + 1) * 1e6
hmGene_tpm <- as.data.frame(hmGene_tpm)
rownames(hmGene_tpm) <- hmGene$geneID

# sebset only KRAB-ZNFs
kznfs_tpm <- hmGene_tpm %>%
    filter(rownames(.) %in% hmKZNFs337$ensembl_gene_id) %>%
    mutate(name = rownames(.)) %>%
    left_join(hmKZNFs337, join_by("name"=="ensembl_gene_id"))
rownames(kznfs_tpm) <- kznfs_tpm$external_gene_name
kznfs_tpm <- kznfs_tpm[,c(1:132)]


# tes convert to tpm
hmTE_temp <- hmTE[,-c(1,2,3)]
te_count <- colSums(hmTE_temp)
te_scale <- te_count / 1e6
hmTE_tpm <- hmTE_temp
hmTE_tpm <- t(t(hmTE_tpm)/ te_scale + 1) * 1e6
hmTE_tpm <- as.data.frame(hmTE_tpm)
rownames(hmTE_tpm) <- hmTE$name

convertIDtoName <- function(df, species){
    df_gene <- df[,c(-1)]
    ensIDlist <- df$geneID
    geneName <- gprofiler2::gconvert(
        query = ensIDlist,
        organism = species,
        target="ENSG",
        mthreshold = Inf,
        filter_na = TRUE
    )

    geneName <- geneName %>% select(c(target, name))

    df_kznf <- df %>%
        inner_join(geneName, join_by(geneID==target)) %>%
        filter(name %in% hmKZNFs337$external_gene_name)

    rownames(df_kznf) <- df_kznf$name
    df_kznf <- df_kznf[,-c(1, ncol(df_kznf))]
    df_kznf
}

# 2-2 TPM in chimp
ptGene_temp <- ptGene[,c(-1)]
sample_counts <- colSums(ptGene_temp)

scaling_factor <- sample_counts / 1e6

ptGene_tpm <- ptGene_temp
ptGene_tpm <- t(t(ptGene_tpm)/ scaling_factor + 1) * 1e6
ptGene_tpm <- as.data.frame(ptGene_tpm)
ptGene_tpm$geneID <- ptGene$geneID
ptGene_tpm <- ptGene_tpm[,c(ncol(ptGene_tpm), 1:ncol(ptGene_tpm)-1)]

# sebset only KRAB-ZNFs and convert ID to name
pt_kznfs_tpm <- convertIDtoName(ptGene_tpm, "ptroglodytes")

# tes convert to tpm
ptTE_temp <- ptTE[,-c(1,2,3)]
te_count <- colSums(ptTE_temp)
te_scale <- te_count / 1e6
ptTE_tpm <- ptTE_temp
ptTE_tpm <- t(t(ptTE_tpm)/ te_scale + 1) * 1e6
ptTE_tpm <- as.data.frame(ptTE_tpm)
rownames(ptTE_tpm) <- ptTE$name

# 2-3 TPM in bonobo
ppGene_temp <- ppGene[,c(-1)]
sample_counts <- colSums(ppGene_temp)

scaling_factor <- sample_counts / 1e6

ppGene_tpm <- ppGene_temp
ppGene_tpm <- t(t(ppGene_tpm)/ scaling_factor + 1) * 1e6
ppGene_tpm <- as.data.frame(ppGene_tpm)
ppGene_tpm$geneID <- ppGene$geneID
ppGene_tpm <- ppGene_tpm[,c(ncol(ppGene_tpm), 1:ncol(ppGene_tpm)-1)]

# sebset only KRAB-ZNFs and convert ID to name
pp_kznfs_tpm <- convertIDtoName(ppGene_tpm, "ppaniscus")

# tes convert to tpm
ppTE_tmp <- ppTE[!duplicated(ppTE$name), ]
ppTE_temp <- ppTE_tmp[,-c(1,2,3)]
te_count <- colSums(ppTE_temp)
te_scale <- te_count / 1e6
ppTE_tpm <- ppTE_temp
ppTE_tpm <- t(t(ppTE_tpm)/ te_scale + 1) * 1e6
ppTE_tpm <- as.data.frame(ppTE_tpm)
rownames(ppTE_tpm) <- ppTE_tmp$name

# 2-4 TPM in macaque
mmGene_temp <- mmGene[,c(-1)]
sample_counts <- colSums(mmGene_temp)

scaling_factor <- sample_counts / 1e6

mmGene_tpm <- mmGene_temp
mmGene_tpm <- t(t(mmGene_tpm)/ scaling_factor + 1) * 1e6
mmGene_tpm <- as.data.frame(mmGene_tpm)
mmGene_tpm$geneID <- mmGene$geneID
mmGene_tpm <- mmGene_tpm[,c(ncol(mmGene_tpm), 1:ncol(mmGene_tpm)-1)]

# sebset only KRAB-ZNFs and convert ID to name
mm_kznfs_tpm <- convertIDtoName(mmGene_tpm, "mmulatta")

# tes convert to tpm
mmTE_tmp <- mmTE[!duplicated(mmTE$name), ]
mmTE_temp <- mmTE_tmp[,-c(1,2,3)]
te_count <- colSums(mmTE_temp)
te_scale <- te_count / 1e6
mmTE_tpm <- mmTE_temp
mmTE_tpm <- t(t(mmTE_tpm)/ te_scale + 1) * 1e6
mmTE_tpm <- as.data.frame(mmTE_tpm)
rownames(mmTE_tpm) <- mmTE_tmp$name

# overlap kznfs
h_kznf <- rownames(kznfs_tpm) #337
pt_kznf <- rownames(pt_kznfs_tpm) #292
pp_kznf <- rownames(pp_kznfs_tpm) #216
mm_kznf <- rownames(mm_kznfs_tpm) #245

# overlap TEs
overlapped_kznfs <- Reduce(intersect, list(h_kznf, pt_kznf, pp_kznf, mm_kznf)) #178 overlapped KRAB-ZNFs
overlapped_TEs <- Reduce(intersect, list(rownames(hmTE_tpm), rownames(ptTE_tpm),
                                         rownames(ppTE_tpm), rownames(mmTE_tpm))) #836 overlapped TEs

# we first get the cluster1 and cluster2 ID in each species
hmc1_id <- df_meta %>% filter(cluster=="cluster1" & Organism=="Homo sapiens")
hmc2_id <- df_meta %>% filter(cluster=="cluster2" & Organism=="Homo sapiens")

ptc1_id <- df_meta %>% filter(cluster=="cluster1" & Organism=="Pan troglodytes")
ptc2_id <- df_meta %>% filter(cluster=="cluster2" & Organism=="Pan troglodytes")

ppc1_id <- df_meta %>% filter(cluster=="cluster1" & Organism=="Pan paniscus")
ppc2_id <- df_meta %>% filter(cluster=="cluster2" & Organism=="Pan paniscus")

mmc1_id <- df_meta %>% filter(cluster=="cluster1" & Organism=="Macaca mulatta")
mmc2_id <- df_meta %>% filter(cluster=="cluster2" & Organism=="Macaca mulatta")

human_cor <- function(dfgene, dfte, idlist){

    corr_list <- list()

    indv_list <- c("ha", "hb", "hc", "hd")

    for (id in indv_list) {
        id_select <- idlist %>% filter(!individual %in% id)
        gene <- dfgene[overlapped_kznfs, id_select$Run]
        te <- dfte[overlapped_TEs, id_select$Run]

        corr <- corrOrthologTE(
            geneInput = gene,
            teInput = te,
            numCore = 3
        )

        corr_list <- append(corr_list, list(corr))

    }

    corr_list

}

hmc1_corr_leaveOne <- human_cor(kznfs_tpm, hmTE_tpm, hmc1_id)
hmc2_corr_leaveOne <- human_cor(kznfs_tpm, hmTE_tpm, hmc2_id)

# human c1
c1_a <- hmc1_corr_leaveOne[[1]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c1_b <- hmc1_corr_leaveOne[[2]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c1_c <- hmc1_corr_leaveOne[[3]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c1_d <- hmc1_corr_leaveOne[[4]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))

c1_overlap <- Reduce(intersect, list(c1_a$pair, c1_b$pair, c1_c$pair, c1_d$pair)) #31891

c2_a <- hmc2_corr_leaveOne[[1]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c2_b <- hmc2_corr_leaveOne[[2]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c2_c <- hmc2_corr_leaveOne[[3]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))
c2_d <- hmc2_corr_leaveOne[[4]] %>% filter(padj<0.01 & abs(coef)>0.4) %>% mutate(pair=paste0(geneName, ":", teName))

c2_overlap <- Reduce(intersect, list(c2_a$pair, c2_b$pair, c2_c$pair, c2_d$pair)) #4444


corr_set <- function(df, filter_list){

    all <- df %>% filter(pair %in% filter_list)
    pos <- all %>% filter(coef > 0)
    neg <- all %>% filter(coef < 0)

    result <- list(
        all = all$pair,
        pos = pos$pair,
        neg = neg$pair
    )

    result

}



ptc1_corr <- corrOrthologTE(
    geneInput = pt_kznfs_tpm[overlapped_kznfs, ptc1_id$Run],
    teInput = ptTE_tpm[overlapped_TEs, ptc1_id$Run],
    numCore = 3
)

ptc2_corr <- corrOrthologTE(
    geneInput = pt_kznfs_tpm[overlapped_kznfs, ptc2_id$Run],
    teInput = ptTE_tpm[overlapped_TEs, ptc2_id$Run],
    numCore = 3
)

ppc1_corr <- corrOrthologTE(
    geneInput = pp_kznfs_tpm[overlapped_kznfs, ppc1_id$Run],
    teInput = ppTE_tpm[overlapped_TEs, ppc1_id$Run],
    numCore = 3
)

ppc2_corr <- corrOrthologTE(
    geneInput = pp_kznfs_tpm[overlapped_kznfs, ppc2_id$Run],
    teInput = ppTE_tpm[overlapped_TEs, ppc2_id$Run],
    numCore = 3
)

mmc1_corr <- corrOrthologTE(
    geneInput = mm_kznfs_tpm[overlapped_kznfs, mmc1_id$Run],
    teInput = mmTE_tpm[overlapped_TEs, mmc1_id$Run],
    numCore =3
)

mmc2_corr <- corrOrthologTE(
    geneInput = mm_kznfs_tpm[overlapped_kznfs, mmc2_id$Run],
    teInput = mmTE_tpm[overlapped_TEs, mmc2_id$Run],
    numCore =3
)

# This function is for filling the excel table:
print_sig <- function(df){
    pos <- df %>% filter(padj<0.01 & coef>0.4) %>% nrow()
    neg <- df %>% filter(padj<0.01 & coef< -0.4) %>% nrow()

    print(paste0("positive: ", pos))
    print(paste0("negative: ", neg))
}

print_unique_len <- function(df){
    df_sig <- df %>% filter(padj<0.01 & abs(coef)>0.4)
    num_g <- length(unique(df_sig$geneName))
    num_t <- length(unique(df_sig$teName))

    print(paste0("krab-znf unique count: ", num_g))
    print(paste0("te unique count: ", num_t))

}

corr_obj <- function(df){

    df_sig <- df %>%
        filter(abs(coef)>0.4 & padj<0.01) %>%
        mutate(pair=paste0(teName, ":", geneName))

    all <- df_sig$pair
    pos <- df_sig %>% filter(coef>0)
    neg <- df_sig %>% filter(coef<0)

    result <- list(
        pos=pos$pair,
        neg=neg$pair
    )

    result
}

# cluster 1
hm_cluster1 <- c1_a %>%
    filter(pair %in% c1_overlap) # save human cluster1 using this
hm_cluster1 <- corr_obj(hm_cluster1)

chimp_cluster1 <- corr_obj(ptc1_corr)
bonobo_cluster1 <- corr_obj(ppc1_corr)
maca_cluster1 <- corr_obj(mmc1_corr)

pbd_c1 <- list(
    hs_p = hm_cluster1$pos,
    hs_n = hm_cluster1$neg,
    pt_p = chimp_cluster1$pos,
    pt_n = chimp_cluster1$neg,
    pp_p = bonobo_cluster1$pos,
    pp_n = bonobo_cluster1$neg,
    mm_p = maca_cluster1$pos,
    mm_n = maca_cluster1$neg)

m_pbd_c1 <- make_comb_mat(pbd_c1)

png("../figures/upset_plot_c1_pbd.png", width=6, height=5, units="in", res=400)
u_pbd_c1 <- UpSet(m_pbd_c1, comb_order=order(-comb_size(m_pbd_c1)),
                  top_annotation=upset_top_annotation(m_pbd_c1, add_numbers=TRUE))
print(u_pbd_c1)
dev.off()

# cluster 2
hm_cluster2 <- c2_a %>%
    filter(pair %in% c2_overlap)
hm_cluster2 <- corr_obj(hm_cluster2)

chimp_cluster2 <- corr_obj(ptc2_corr)
bonobo_cluster2 <- corr_obj(ppc2_corr)
maca_cluster2 <- corr_obj(mmc2_corr)

pbd_c2 <- list(
    hs_p = hm_cluster2$pos,
    hs_n = hm_cluster2$neg,
    pt_p = chimp_cluster2$pos,
    pt_n = chimp_cluster2$neg,
    pp_p = bonobo_cluster2$pos,
    pp_n = bonobo_cluster2$neg,
    mm_p = maca_cluster2$pos,
    mm_n = maca_cluster2$neg)

m_pbd_c2 <- make_comb_mat(pbd_c2)

png("../figures/upset_plot_c2_pbd.png", width=6, height=5, units="in", res=400)
u_pbd_c2 <- UpSet(m_pbd_c2, comb_order=order(-comb_size(m_pbd_c2)),
                  top_annotation=upset_top_annotation(m_pbd_c2, add_numbers=TRUE))
print(u_pbd_c2)
dev.off()

human_cor_all4 <- function(dfgene, dfte, idlist){

    gene <- dfgene[overlapped_kznfs, idlist$Run]
    te <- dfte[overlapped_TEs, idlist$Run]

    corr <- corrOrthologTE(
        geneInput = gene,
        teInput = te,
        numCore = 3
    )

    corr
}

hmc1_corr_all4 <- human_cor_all4(kznfs_tpm, hmTE_tpm, hmc1_id)
hmc2_corr_all4 <- human_cor_all4(kznfs_tpm, hmTE_tpm, hmc2_id)
print_sig(hmc1_corr_all4)
print_unique_len(hmc1_corr_all4)
print_sig(hmc2_corr_all4)
print_unique_len(hmc2_corr_all4)

# human
hmc1_corr <- c1_a %>%
    filter(pair %in% c1_overlap) %>%
    mutate(pair = paste0(teName, ":", geneName))

hmc2_corr <- c2_a %>%
    filter(pair %in% c2_overlap) %>%
    mutate(pair = paste0(teName, ":", geneName))

# chimp
ptc1_corr <- ptc1_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

ptc2_corr <- ptc2_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

# bonobo
ppc1_corr <- ppc1_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

ppc2_corr <- ppc2_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

# macaque
mmc1_corr <- mmc1_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

mmc2_corr <- mmc2_corr %>%
    filter(padj<0.01 & abs(coef)>0.4) %>%
    mutate(pair = paste0(teName, ":", geneName))

# filter
hmc1_corr_all4 <- hmc1_corr_all4 %>%
    mutate(pair=paste0(teName, ":", geneName)) %>%
    filter(pair %in% hmc1_corr$pair)

hmc2_corr_all4 <- hmc2_corr_all4 %>%
    mutate(pair=paste0(teName, ":", geneName)) %>%
    filter(pair %in% hmc2_corr$pair)

pbd_obj <- list(
    "hmc1_corr"=hmc1_corr_all4,
    "hmc2_corr"=hmc2_corr_all4,
    "ptc1_corr"=ptc1_corr,
    "ptc2_corr"=ptc2_corr,
    "ppc1_corr"=ppc1_corr,
    "ppc2_corr"=ppc2_corr,
    "mmc1_corr"=mmc1_corr,
    "mmc2_corr"=mmc2_corr
)

saveRDS(pbd_obj, file="../data/pbd_obj.rds")

pbd_obj <- readRDS("data/pbd_obj.rds")

# cluster 1
y_kznf <- kznf_infer %>%filter(age=="young")

hs_p <- hmc1_corr_all4 %>% filter(coef>0) %>% mutate(species="Hs")
pp_n <- ppc1_corr %>% filter(coef<0) %>% mutate(species="Pp")

intersect_hs_pp_276 <- intersect(hs_p$pair, pp_n$pair) #276
df_276 <- rbind(hs_p[hs_p$pair %in% intersect_hs_pp_276,],
                pp_n[pp_n$pair %in% intersect_hs_pp_276,])

hs_p_276 <- hs_p %>%
    filter(pair %in% intersect_hs_pp_276) %>%
    mutate(teAge = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(kznfAge = ifelse(geneName %in% y_kznf$external_gene_name, "young", "old")) %>%
    mutate(link = ifelse(teAge=="old" & kznfAge=="old", "old", "young")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

# create a 276 network
library(RCy3)
library(netZooR)

c1.276.kznf <- hs_p_276 %>% select(c(1,9)) %>% unique() #15
colnames(c1.276.kznf) <- c("id", "age")
c1.276.te <- hs_p_276 %>% select(c(2,8)) %>% unique() #180
colnames(c1.276.te) <- c("id", "age")

c1.276.node <- rbind(c1.276.kznf, c1.276.te)
c1.276.link <- hs_p_276[,c(1,2,3,10)]
colnames(c1.276.link) <- c("source", "target", "coefficient", "age")

createNetworkFromDataFrames(c1.276.node, c1.276.link)

library(ggplot2)
hs_pp_276_teFamily <- data.frame(t(table(hs_p_276[,c(10, 11)])))
colnames(hs_pp_276_teFamily) <- c("TE_family", "link", "Count")

g_hs_pp_276 <- ggplot(hs_pp_276_teFamily, aes(x=TE_family, y=Count, fill=link)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("old" = "#318ce7", "young" = "#f6b26b")) +
    coord_flip() +
    theme_bw()

ggsave(g_hs_pp_276, file="../figures/pbd_network_comparison/hs_pp_276_TEfamily.jpg", width=4, height=5)
ggsave(g_hs_pp_276, file="../figures/pbd_network_comparison/hs_pp_276_TEfamily.svg", width=4, height=5)

pbd_obj <- readRDS("data/pbd_obj.rds")
y_kznf <- kznf_infer %>% filter(age=="young")

hsc1 <- pbd_obj$hmc1_corr %>%
    filter(!pair %in% pbd_obj$ptc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$ppc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$mmc1_corr$pair) %>%
    mutate(te_age = ifelse(teName %in% te_infer$teName, "y", "o")) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer$external_gene_name, "y", "o")) %>%
    mutate(link_age = ifelse(te_age=="o" & kznf_age=="o", "o", "y")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

hsc1_specific_p <- hsc1 %>% filter(coef > 0)
hsc1_specific_n <- hsc1 %>% filter(coef < 0)

df_hs_p <- table(hsc1_specific_p$family_id, hsc1_specific_p$link_age) %>%
    data.frame() %>%
    mutate(coef = "positive")

df_hs_n <- table(hsc1_specific_n$family_id, hsc1_specific_n$link_age) %>%
    data.frame() %>%
    mutate(coef = "negative")

df_hs_merge <- rbind(df_hs_p, df_hs_n)
colnames(df_hs_merge) <- c("teFamily", "age", "count", "coef")

# normalize using the number of TE family
df_te <- hg19rmsk_info %>%
    filter(gene_id %in% overlapped_TEs)

df_teFamily <- table(df_te$family_id) %>%
    data.frame()

# merge with the raw counts and family
df_hs_merge_norm <- df_hs_merge %>%
    left_join(df_teFamily, join_by(teFamily==Var1)) %>%
    mutate(norm_count = log(count / Freq * 100)) %>%
    arrange(desc(norm_count))


g_h_specific <- ggplot(df_hs_merge_norm, aes(x=coef, y=teFamily)) +
    geom_point(aes(color=norm_count, size=norm_count, alpha=0.85)) +
    scale_colour_gradient(low="#c3d0d1", high="#6b8b8e") +
    scale_alpha(guide = 'none') +
    labs(x="", y="") +
    scale_x_discrete(labels=c("y-negative", "y-positive")) +
    ggtitle("Human specific correlation counts in TE family") +
    theme_bw()

ggsave(g_h_specific, file="../figures/pbd_hs_specific_link_TEfamilyCount.jpg", dpi=400, width=4)

ppc1 <- pbd_obj$ppc1_corr %>%
    filter(!pair %in% pbd_obj$hmc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$ptc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$mmc1_corr$pair) %>%
    mutate(te_age = ifelse(teName %in% te_infer$teName, "y", "o")) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer$external_gene_name, "y", "o")) %>%
    mutate(link_age = ifelse(te_age=="o" & kznf_age=="o", "o", "y")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

ppc1_specific_p <- ppc1 %>% filter(coef > 0) #3386
ppc1_specific_n <- ppc1 %>% filter(coef < 0) #3297

# separate positive and negative link
df_pp_p <- table(ppc1_specific_p$family_id, ppc1_specific_p$link_age) %>%
    data.frame() %>%
    mutate(coef = "positive")

df_pp_n <- table(ppc1_specific_n$family_id, ppc1_specific_n$link_age) %>%
    data.frame() %>%
    mutate(coef = "negative")

df_pp_merge <- rbind(df_pp_p, df_pp_n)
colnames(df_pp_merge) <- c("teFamily", "age", "count", "coef")


# merge with the raw counts and family
df_pp_merge_norm <- df_pp_merge %>%
    left_join(df_teFamily, join_by(teFamily==Var1)) %>%
    mutate(norm_count = log(count / Freq * 100)) %>%
    arrange(desc(norm_count))

g_p_specific <- ggplot(df_pp_merge_norm, aes(x=coef, y=teFamily)) +
    geom_point(aes(color=norm_count, size=norm_count, alpha=0.85)) +
    scale_colour_gradient(low="#c3d0d1", high="#6b8b8e") +
    scale_alpha(guide = 'none') +
    labs(x="", y="") +
    scale_x_discrete(labels=c("y-negative", "y-positive")) +
    ggtitle("Bonobo specific correlation counts in TE family") +
    theme_bw()

ggsave(g_p_specific, file="../figures/pbd_pp_specific_link_TEfamilyCount.jpg", dpi=400, width=4)

From the results, we could also find that bonobo has similar pattern as human, for example, all the correlation are young links and they are contributed from the young KRAB-ZNFs with old TEs.

Next, we go for the chimpanzee and macaque


ptc1 <- pbd_obj$ptc1_corr %>%
    filter(!pair %in% pbd_obj$hmc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$ppc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$mmc1_corr$pair) %>%
    mutate(te_age = ifelse(teName %in% te_infer$teName, "y", "o")) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer$external_gene_name, "y", "o")) %>%
    mutate(link_age = ifelse(te_age=="o" & kznf_age=="o", "o", "y")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

ptc1_specific_p <- ptc1 %>% filter(coef > 0) #20
ptc1_specific_n <- ptc1 %>% filter(coef < 0) #18

# separate positive and negative link
df_pt_p <- table(ptc1_specific_p$family_id, ptc1_specific_p$link_age) %>%
    data.frame() %>%
    mutate(coef = "positive")

df_pt_n <- table(ptc1_specific_n$family_id, ptc1_specific_n$link_age) %>%
    data.frame() %>%
    mutate(coef = "negative")

df_pt_merge <- rbind(df_pt_p, df_pt_n)
colnames(df_pt_merge) <- c("teFamily", "age", "count", "coef")


# merge with the raw counts and family
df_pt_merge_norm <- df_pt_merge %>%
    left_join(df_teFamily, join_by(teFamily==Var1)) %>%
    mutate(norm_count = log(count / Freq * 100)) %>%
    arrange(desc(norm_count))

# visualization
g_pt_specific <- ggplot(df_pt_merge_norm, aes(x=coef, y=teFamily)) +
    geom_point(aes(color=norm_count, size=norm_count, alpha=0.85)) +
    scale_colour_gradient(low="#c3d0d1", high="#6b8b8e") +
    scale_alpha(guide = 'none') +
    labs(x="", y="") +
    scale_x_discrete(labels=c("y-negative", "y-positive")) +
    ggtitle("Chimpanzee-specific correlation counts in TE family") +
    theme_bw()

ggsave(g_pt_specific, file="../figures/pbd_pt_specific_link_TEfamilyCount.jpg", dpi=400, width=4)

mmc1 <- pbd_obj$mmc1_corr %>%
    filter(!pair %in% pbd_obj$hmc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$ppc1_corr$pair) %>%
    filter(!pair %in% pbd_obj$ptc1_corr$pair) %>%
    mutate(te_age = ifelse(teName %in% te_infer$teName, "y", "o")) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer$external_gene_name, "y", "o")) %>%
    mutate(link_age = ifelse(te_age=="o" & kznf_age=="o", "o", "y")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

mmc1_specific_p <- mmc1 %>% filter(coef > 0) #3
mmc1_specific_n <- mmc1 %>% filter(coef < 0) #1

# separate positive and negative link
df_mm_p <- table(mmc1_specific_p$family_id, mmc1_specific_p$link_age) %>%
    data.frame() %>%
    mutate(coef = "positive")

df_mm_n <- table(mmc1_specific_n$family_id, mmc1_specific_n$link_age) %>%
    data.frame() %>%
    mutate(coef = "negative")

df_mm_merge <- rbind(df_mm_p, df_mm_n)
colnames(df_mm_merge) <- c("teFamily", "age", "count", "coef")


# merge with the raw counts and family
df_mm_merge_norm <- df_mm_merge %>%
    left_join(df_teFamily, join_by(teFamily==Var1)) %>%
    mutate(norm_count = log(count / Freq * 100)) %>%
    arrange(desc(norm_count))

# visualization
g_mm_specific <- ggplot(df_mm_merge_norm, aes(x=coef, y=teFamily)) +
    geom_point(aes(color=norm_count, size=norm_count, alpha=0.85)) +
    scale_colour_gradient(low="#c3d0d1", high="#6b8b8e") +
    scale_alpha(guide = 'none') +
    labs(x="", y="") +
    scale_x_discrete(labels=c("y-negative", "y-positive")) +
    ggtitle("macaque-specific correlation counts in TE family") +
    theme_bw()

ggsave(g_mm_specific, file="../figures/pbd_mm_specific_link_TEfamilyCount.jpg", dpi=400, width=4)

hsc1_n <- pbd_obj$hmc1_corr %>% filter(coef<0)
ppc1_p <- pbd_obj$ppc1_corr %>% filter(coef>0)

df_hs_pp_2 <- hsc1_n %>% filter(pair %in% ppc1_p$pair) #17

df_hs_pp_2 <- df_hs_pp_2 %>%
    mutate(te_age = ifelse(teName %in% te_infer$teName, "y", "o")) %>%
    mutate(kznf_age = ifelse(geneName %in% kznf_infer$external_gene_name, "y", "o")) %>%
    mutate(link_age = ifelse(te_age=="o" & kznf_age=="o", "o", "y")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

c1.17.kznf <- df_hs_pp_2 %>% select(c(1,8)) %>% unique() #4
colnames(c1.17.kznf) <- c("id", "age")
c1.17.te <- df_hs_pp_2 %>% select(c(2,7)) %>% unique() #13
colnames(c1.17.te) <- c("id", "age")

c1.17.node <- rbind(c1.17.kznf, c1.17.te)


c1.17.link <- df_hs_pp_2[,c(1,2,3,9)]
colnames(c1.17.link) <- c("source", "target", "coefficient", "age")


hs_pp_17_teFamily <- data.frame(t(table(df_hs_pp_2[,c(10, 9)])))
colnames(hs_pp_17_teFamily) <- c("link", "TE_family", "Count")

g_hs_pp_17 <- ggplot(hs_pp_17_teFamily, aes(x=TE_family, y=Count, fill=link)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("o" = "#318ce7", "y" = "#f6b26b")) +
    coord_flip() +
    theme_bw()

ggsave(g_hs_pp_17, file="../figures/pbd_network_comparison/hs_pp_17_TEfamily.jpg", width=4, height=5)
ggsave(g_hs_pp_17, file="../figures/pbd_network_comparison/hs_pp_17_TEfamily.svg", width=4, height=5)

# in cluster 1
hc1 <- pbd_obj$hmc1_corr
ppc1 <- pbd_obj$ppc1_corr

ppc1_neg <- ppc1 %>% filter(coef<0)

corr_276_list <- hc1 %>%
    filter(pair %in% ppc1_neg$pair) %>%
    filter(padj<0.01 & coef>0)

corr_276_list <- corr_276_list %>%
    mutate(teAge = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(kznfAge = ifelse(geneName %in% y_kznf$external_gene_name, "young", "old")) %>%
    mutate(link = ifelse(teAge=="old" & kznfAge=="old", "old", "young")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

# load chimpanzee and macaque, we want the unfilter data
load("~/github/primate_network/data/c1c2_NHPs.RData")
ptc1 <- c1c2_NHPs$pt_c1
ptc1_filter <- ptc1 %>%
    filter(pair %in% corr_276_list$pair) #76 positive and 200 negative, but not significant

mmc1 <- c1c2_NHPs$mm_c1
mmc1_filter <- mmc1 %>%
    filter(pair %in% corr_276_list$pair) #144 positive and 132 negative, but not significant


