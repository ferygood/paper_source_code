library(ggplot2)

hs_pp_276 <- read.csv("tables/hs_pp_276_c1.csv")

hs_pp_276_teFamily <- data.frame(t(table(hs_pp_276[,c(10, 11)])))
colnames(hs_pp_276_teFamily) <- c("TE_family", "link", "Count")

g_hs_pp_276 <- ggplot(hs_pp_276_teFamily, aes(x=TE_family, y=Count, fill=link)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("old" = "#318ce7", "young" = "#f6b26b")) +
    coord_flip() +
    theme_bw()

ggsave(g_hs_pp_276, file="figures/JPG/4E_hs_pp_276_TEfamily.jpg", width=4, height=5)
ggsave(g_hs_pp_276, file="figures/SVG/4E_hs_pp_276_TEfamily.svg", width=4, height=5)
