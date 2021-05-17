# functional enrichment for GTEx paper

### Guangchuang Yu

#### 07/15, 2015

------------------------------------------------------------------------

We re-run the comparison with latest version of DAVID (v6.8) and
clusterProfiler (v4) in May 17, 2021. The result can be found in the
[new-result-2021](https://github.com/GuangchuangYu/enrichment4GTEx_clusterProfiler/tree/master/new-result-2021)
foulder.

------------------------------------------------------------------------

Refer to the [issue](https://github.com/GuangchuangYu/DOSE/issues/6),
this document will reproduce functional analysis using
[DOSE](http://www.bioconductor/packages/DOSE) and
[clusterProfiler](http://www.bioconductor.org/packages/clusterProfiler)
packages.

The paper and supplemental files are located at [paper folder](paper).

In page 63 of [supplemental file](paper/Mele.SM.pdf), the authors
mentioned that *genes with high contribution of individuals to splicing
variation* were used in GO enrichment analysis.

![](figures/Screenshot%202015-07-15%2020.03.05.png)

The **Individual and Tissue contribution to variation of splicing** was
stored in [supplemental table
17](paper/LargeSupplementTABLES_May1st_2015.xlsx), which was exported as
a [csv file](tableS17.csv).

``` r
require(magrittr)
require(DOSE)
```

    ## Loading required package: DOSE

    ## 

    ## DOSE v3.17.2  For help: https://guangchuangyu.github.io/software/DOSE
    ## 
    ## If you use DOSE in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609

``` r
require(RDAVIDWebService)
```

    ## Loading required package: RDAVIDWebService

    ## Loading required package: graph

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: GOstats

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: Category

    ## Loading required package: stats4

    ## Loading required package: AnnotationDbi

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

    ## 
    ## Attaching package: 'GOstats'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     makeGOGraph

    ## Loading required package: ggplot2

``` r
require(clusterProfiler)
```

    ## Loading required package: clusterProfiler

    ## clusterProfiler v3.99.2  For help: https://guangchuangyu.github.io/software/clusterProfiler
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

``` r
table17 <- read.csv("paper/tableS17.csv")

gene <- with(table17, gene_id[ Individuals > quantile(Individuals, 0.98)])
gene %<>% as.character %>% gsub("\\.\\d+", "", .)

head(gene)
```

    ## [1] "ENSG00000002586" "ENSG00000011007" "ENSG00000014164" "ENSG00000042753"
    ## [5] "ENSG00000051620" "ENSG00000055130"

``` r
length(gene)
```

    ## [1] 212

``` r
save(gene, file="cache/gene.rda")
```

The authors did not mention how they define ***high contribution***,
here I use the top 2% of genes and get 212 genes with 168 genes that
could be mapped in the DAVID database. Slightly large than the number
(139) reported in [supplemental file](paper/Mele.SM.pdf).

The GO enrichment result reported in [supplemental table
18](paper/LargeSupplementTABLES_May1st_2015.xlsx) is:

![](figures/Screenshot%202015-07-15%2019.49.01.png)

I pasted the genes into DAVID and got similar results.

![](figures/Screenshot%202015-07-15%2021.10.11.png)
![](figures/Screenshot%202015-07-15%2021.10.28.png)

Although the gene list I selected here is slighly different from the one
author selected (which we don’t know), it can reproduce the results
reported in the paper.

Here, I use these genes to perform enrichment analyses and compare
results from DAVID and clusterProfiler.

-   [GO (BP) enrichment analysis](GO_BP.md)
-   [KEGG enrichment analysis](KEGG.md)
-   [DO enrichment analysis](DO.md)
-   [Reactome enrichment analysis](Reactome.md)
