library(twice)
library(dplyr)
library(ComplexHeatmap)

data("hg19rmsk_info")
load("data/mayoTEKRABber_balance.RData") # load mayo RNA seq result

c1_confirm <- read.csv("tables/hmc1_sig_forNetwork.csv")
c2_confirm <- read.csv("tables/hmc2_sig_forNetwork.csv")
confirm_mix <- c(c1_confirm$pair, c2_confirm$pair)
confirm_mix <- confirm_mix[!duplicated(confirm_mix)] # merge the confirm list, total 901

chipexo <- read.csv("tables/kznfs_TEs_ChIP_exo_modified.csv")
chipexo_process <- chipexo %>%
    mutate(pair=paste0(teName, ":", geneName))

# Find overlap with ChIP-exo
tcx_control <- mayoTEKRABber$tcxControlCorr %>%
    mutate(pair = paste0(teName, ":", geneName)) %>%
    filter(pair %in% chipexo_process$pair & padj<0.01 & abs(coef)>=0.4)  #59

tcx_ad <- mayoTEKRABber$tcxADCorr %>%
    mutate(pair = paste0(teName, ":", geneName)) %>%
    filter(pair %in% chipexo_process$pair & padj<0.01 & abs(coef)>=0.4)  #3

cbe_control <- mayoTEKRABber$cbeControlCorr %>%
    mutate(pair = paste0(teName, ":", geneName)) %>%
    filter(pair %in% chipexo_process$pair & padj<0.01 & abs(coef)>=0.4)  #91

cbe_ad <- mayoTEKRABber$cbeADCorr %>%
    mutate(pair = paste0(teName, ":", geneName)) %>%
    filter(pair %in% chipexo_process$pair & padj<0.01 & abs(coef)>=0.4)  #35

# upsetplot for the overlap of mayo data
lt_mayo <- list(
    tcx_control = tcx_control$pair,
    tcx_AD = tcx_ad$pair,
    cbe_control = cbe_control$pair,
    cbe_AD = cbe_ad$pair)

m_mayo <- make_comb_mat(lt_mayo)

png("figures/JPG/6A_upset_plot_mayo.jpg", width=5, height=3, units="in", res=400)
u_mayo <- UpSet(m_mayo, comb_order=order(-comb_size(m_mayo)), top_annotation=upset_top_annotation(m_mayo, add_numbers=TRUE))
print(u_mayo)
dev.off()
