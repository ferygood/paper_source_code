# fig 3B
# From previous results, we know that only cluster1 and cluster2 have correlations detected.
# In this section, I am going to analyze the 1000 iterations result:
# (1) check with pvalue (2) check the idea range of coefficient
# (3) overlapped the results with ChIP-exo. Then we will finally get how many
# correlations that we want to go further.

library(dplyr)
library(purrr)
library(ggplot2)

corr_c1 <- read.csv("tables/hsC1_corr_sig.csv")
hsC1_sig <- corr_c1 %>% filter(padj<0.01) #127797

corr_c2 <- read.csv("tables/hsC2_corr_sig.csv")
hsC2_sig <- corr_c2 %>% filter(padj<0.01) #49770

chip_exo <- read.csv("tables/kznfs_TEs_ChIP_exo_modified.csv")
chip_exo_pair <- chip_exo %>%
    mutate(pair = paste0(teName, ":", geneName))

overlapped <- intersect(hsC1_sig$pair, chip_exo_pair$pair) #869

# create a dataframe, col1=count, col2=coefficient
df_temp <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

df_pos <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

df_neg <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

for (i in seq(from=0, to=1, by=0.05)){

    filter_c1 <- hsC1_sig %>% filter(pair %in% overlapped)

    count <- filter_c1 %>% filter(abs(coef)>=i) %>% nrow()
    pos_count <- filter_c1 %>% filter(coef>=i) %>% nrow()
    neg_count <- filter_c1 %>% filter(coef<= i*(-1)) %>% nrow()
    coef <- i

    df_temp <- rbind(df_temp, c(count, coef, "all"))
    df_pos <- rbind(df_pos, c(pos_count, coef, "positive"))
    df_neg <- rbind(df_neg, c(neg_count, coef, "negative"))

}

colnames(df_temp) <- c("count", "coefficient", "group")
colnames(df_pos) <- c("count", "coefficient", "group")
colnames(df_neg) <- c("count", "coefficient", "group")

df_merge <- rbind(df_temp, df_pos, df_neg)
df_merge$count <- as.numeric(df_merge$count)
df_merge$coefficient <- as.numeric(df_merge$coefficient)
df_merge$group <- factor(df_merge$group, levels = c("all", "positive", "negative"))

# cluster 1, primary and secondary cortices
gc1 <- ggplot(df_merge, aes(x=coefficient, y=count, color=group)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks=seq(0, 1, by=0.1)) +
    scale_y_continuous(breaks=seq(0, 1000, by=50)) +
    theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=1))+
    ylab("overlapped correlations") +
    ggtitle("padj<0.01 (Primary and Secondary cortices)") +
    scale_color_manual(values=c("black", "#cc3433", "#336699")) +
    theme_bw() +
    theme(text = element_text(size=8))

ggsave(gc1, file="figures/JPG/3B_correlation_check_c1.jpg", dpi=400, width=4, height=4)

# cluster 2, limbic and association cortices
overlapped_c2 <- intersect(hsC2_sig$pair, chip_exo_pair$pair) #399

# create a dataframe, col1=count, col2=coefficient
df_temp <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

df_pos <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

df_neg <- data.frame(
    count=numeric(),
    coef=numeric(),
    group=character()
)

for (i in seq(from=0, to=1, by=0.05)){

    filter_c2 <- hsC2_sig %>% filter(pair %in% overlapped)

    count <- filter_c2 %>% filter(abs(coef)>=i) %>% nrow()
    pos_count <- filter_c2 %>% filter(coef>=i) %>% nrow()
    neg_count <- filter_c2 %>% filter(coef<= i*(-1)) %>% nrow()
    coef <- i

    df_temp <- rbind(df_temp, c(count, coef, "all"))
    df_pos <- rbind(df_pos, c(pos_count, coef, "positive"))
    df_neg <- rbind(df_neg, c(neg_count, coef, "negative"))

}

colnames(df_temp) <- c("count", "coefficient", "group")
colnames(df_pos) <- c("count", "coefficient", "group")
colnames(df_neg) <- c("count", "coefficient", "group")

df_merge <- rbind(df_temp, df_pos, df_neg)
df_merge$count <- as.numeric(df_merge$count)
df_merge$coefficient <- as.numeric(df_merge$coefficient)
df_merge$group <- factor(df_merge$group, levels = c("all", "positive", "negative"))

gc2 <- ggplot(df_merge, aes(x=coefficient, y=count, color=group)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks=seq(0, 1, by=0.1)) +
    scale_y_continuous(breaks=seq(0, 1000, by=50)) +
    theme(axis.text.x = element_text(angle=45, vjust=0.5, hjust=1))+
    ylab("overlapped correlations") +
    ggtitle("padj<0.01 (Limbic and Association cortices)") +
    scale_color_manual(values=c("black", "#cc3433", "#336699")) +
    theme_bw() +
    theme(text = element_text(size=8))

ggsave(gc2, file="figures/JPG/S3B_correlation_check_c2.jpg", dpi=400, width=4, height=4)
