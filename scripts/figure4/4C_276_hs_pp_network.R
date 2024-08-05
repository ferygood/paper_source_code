library(RCy3)
library(netZooR)
library(twice)

pbd_obj <- readRDS("data/pbd_obj.rds")
data("hg19rmsk_info")

# cluster 1
y_kznf <- kznf_infer %>%filter(age=="young")

hs_p <- pbd_obj$hmc1_corr %>% filter(coef>0) %>% mutate(species="Hs")
pp_n <- pbd_obj$ppc1_corr %>% filter(coef<0) %>% mutate(species="Pp")

intersect_hs_pp_276 <- base::intersect(hs_p$pair, pp_n$pair) #276
df_276 <- rbind(hs_p[hs_p$pair %in% intersect_hs_pp_276,],
                pp_n[pp_n$pair %in% intersect_hs_pp_276,])

hs_p_276 <- hs_p %>%
    filter(pair %in% intersect_hs_pp_276) %>%
    mutate(teAge = ifelse(teName %in% te_infer$NM, "young", "old")) %>%
    mutate(kznfAge = ifelse(geneName %in% y_kznf$external_gene_name, "young", "old")) %>%
    mutate(link = ifelse(teAge=="old" & kznfAge=="old", "old", "young")) %>%
    left_join(hg19rmsk_info[,c(1,2)], join_by(teName==gene_id))

write.table(hs_p_276, file = "tables/hs_pp_276_c1.csv", sep=",", row.names = F)

c1.276.kznf <- hs_p_276 %>% select(c(1,9)) %>% unique() #15
colnames(c1.276.kznf) <- c("id", "age")
c1.276.te <- hs_p_276 %>% select(c(2,8)) %>% unique() #180
colnames(c1.276.te) <- c("id", "age")

c1.276.node <- rbind(c1.276.kznf, c1.276.te)
c1.276.link <- hs_p_276[,c(1,2,3,10)]
colnames(c1.276.link) <- c("source", "target", "coefficient", "age")

createNetworkFromDataFrames(c1.276.node, c1.276.link)
