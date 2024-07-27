source("R/brain_meta.R")

read_rds_files <- function(folder_path){

    file_list <- list.files(path = folder_path, pattern="\\.rds$", full.names = TRUE)
    for (file_path in file_list){
        file_name <- basename(file_path)
        var_name <- gsub("\\.rds$", "", file_name)
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
    }
}


volcano_plot <- function(df){

    v <- EnhancedVolcano::EnhancedVolcano(
        df,
        lab = rownames(df),
        x = 'log2FoldChange',
        y = 'pvalue',
        labSize = 3.0,
        FCcutoff = 1.5,
        pCutoff = 0.05,
        legendPosition = 'none',
        legendLabSize = 10,
        title = "",
        subtitle = "",
        caption = ""
    )

    v

}


compare_violin <- function(kznfs, TEs, brain_region){
    c_kznfs <- kznfs %>%
        filter(region == brain_region)

    c_kznfs$species <- factor(c_kznfs$species,
                              levels=c("Human", "Chimpanzee", "Bonobo", "Macaque"))

    gv <- ggviolin(c_kznfs, x="age", y="value", fill="age",
                   palette=c("#5A8FBB", "#E59E00"),
                   add = "boxplot", add.params = list(fill="white")) +
        stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
        theme(legend.position = "bottom", axis.text.x=element_blank()) +
        facet_wrap(~species, ncol=4) +
        ylab("Log expression") +
        xlab("") +
        ggtitle(paste0("KRAB-ZNFs in ", brain_region))

    #ggsave(filename=gfile, gv, dpi=500, width=10, height=4)

    c_te <- TEs %>%
        filter(region == brain_region)

    c_te$species <- factor(c_te$species,
                           levels=c("Human", "Chimpanzee", "Bonobo", "Macaque"))
    tv <- ggviolin(c_te, x="age", y="value", fill="age",
                   palette=c("#5A8FBB", "#E59E00"),
                   add = "boxplot", add.params = list(fill="white")) +
        stat_compare_means(label="p.signif", label.x=1.5, method="wilcox.test") +
        theme(legend.position = "bottom", axis.text.x=element_blank()) +
        facet_wrap(~species, ncol=4) +
        ylab("Log expression") +
        xlab("") +
        ggtitle(paste0("TEs in ", brain_region))

    #ggsave(filename = tfile, tv, dpi=500, width=10, height=4)

    result <- list("gv"=gv, "tv"=tv)
    result
}
