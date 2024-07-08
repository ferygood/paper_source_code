# Paper title

This repo including the analysis code in the paper, \_\_\_\_\_. You can refer to the script based on the description below.

## preTEKRABber pipeline is provided

You can find the steps including: (1) using fastp to remove low quality reads and trimmed adapters (2) using STAR to align reads to reference genome (3) use TEtranscripts to quantify the expression of genes and transposable elements. [repo link](https://github.com/ferygood/preTEKRABber_pipe)

## in-house software used

To reduce the redundancy of codes, several functions is from R twice package. You can download it from github:

``` r
devtools::install_github("ferygood/twice")
```

## Datasets

Two independent RNA-seq dataset were used in this study. - Cross-species: [GSE127898](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127898) - Control-Alzheimer's disease: [syn5550404](https://www.synapse.org/#!Synapse:syn5550404)

## Scripts

1.  Figure 2A Primate Brain Data t-SNE plot [[link](google.com)]

2.  Violin Plot expression

3.  Percentage barplot

4.  Expression Heatmap

## Contact

Email: [yao-chung.chen\@fu-berlin.de](mailto:yao-chung.chen@fu-berlin.de){.email}
