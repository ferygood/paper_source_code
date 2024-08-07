library(ggplot2)
library(twice)

data(hg19rmsk_info)
pbd_obj <- readRDS("data/pbd_obj.rds")
y_kznf <- kznf_infer %>% filter(age=="young")
overlaps <- readRDS("data/pbd_kznfs_tes_overlap.rds")

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
    filter(gene_id %in% overlaps$overlapped_TEs)

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

ggsave(g_h_specific, file="figures/JPG/4F_pbd_hs_specific_link_TEfamilyCount.jpg",
       dpi=400, width=4, height=6)

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

#ggsave(g_p_specific, file="figures/JPG/pbd_pp_specific_link_TEfamilyCount.jpg",
#       dpi=400, width=4, height=6)

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

#ggsave(g_pt_specific, file="figures/JPG/pbd_pt_specific_link_TEfamilyCount.jpg",
       #dpi=400, width=4, height=6)

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

#ggsave(g_mm_specific, file="figures/JPG/pbd_mm_specific_link_TEfamilyCount.jpg",
       #dpi=400, width=4, height=6)
