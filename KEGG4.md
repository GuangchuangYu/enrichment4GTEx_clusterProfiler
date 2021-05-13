## KEGG enrichment with DAVID

``` r
require(magrittr)
require(DOSE)
require(RDAVIDWebService)
require(clusterProfiler)
```

``` r
load("cache/gene.rda")
david_KEGG4 <- enrichDAVID(gene, idType="ENSEMBL_GENE_ID", annotation="KEGG_PATHWAY", david.user="gcyu@connect.hku.hk")
save(david_KEGG4, file = "david_KEGG4.rda")
```

``` r
summary(david_KEGG)[, -8]
```

    ## Warning in summary(david_KEGG): summary method to convert the object to
    ## data.frame is deprecated, please use as.data.frame instead.

    ##                ID Description GeneRatio BgRatio       pvalue     p.adjust
    ## hsa03010 hsa03010    Ribosome     23/83 87/5085 1.804017e-21 1.731856e-19
    ##                qvalue Count
    ## hsa03010 7.595859e-21    23

Ribosome is the only term reported by DAVID consistent with the result
reported in the supplemental file.

## KEGG enrichment with clusterProfiler

``` r
eg=bitr(gene, "ENSEMBL", "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]
clusterProfiler_KEGG4 <- enrichKEGG(eg)
save(clusterProfiler_KEGG4, file = "clusterProfiler_KEGG4.rda")
```

``` r
dim(summary(clusterProfiler_KEGG))
```

    ## Warning in summary(clusterProfiler_KEGG): summary method to convert the object
    ## to data.frame is deprecated, please use as.data.frame instead.

    ## [1] 9 9

``` r
head(summary(clusterProfiler_KEGG)[, -8])
```

    ## Warning in summary(clusterProfiler_KEGG): summary method to convert the object
    ## to data.frame is deprecated, please use as.data.frame instead.

    ##                ID                         Description GeneRatio  BgRatio
    ## hsa03010 hsa03010                            Ribosome    25/104 137/6895
    ## hsa04612 hsa04612 Antigen processing and presentation     8/104  77/6895
    ## hsa05330 hsa05330                 Allograft rejection     5/104  37/6895
    ## hsa05332 hsa05332           Graft-versus-host disease     5/104  41/6895
    ## hsa04940 hsa04940            Type I diabetes mellitus     5/104  43/6895
    ## hsa05320 hsa05320          Autoimmune thyroid disease     5/104  52/6895
    ##                pvalue     p.adjust       qvalue Count
    ## hsa03010 6.111415e-21 1.008384e-18 9.392280e-19    25
    ## hsa04612 1.819033e-05 1.500702e-03 1.397783e-03     8
    ## hsa05330 2.105386e-04 1.157962e-02 1.078549e-02     5
    ## hsa05332 3.450666e-04 1.423400e-02 1.325782e-02     5
    ## hsa04940 4.327776e-04 1.428166e-02 1.330222e-02     5
    ## hsa05320 1.049593e-03 2.886381e-02 2.688432e-02     5

``` r
dotplot(clusterProfiler_KEGG, showCategory=9)
```

![](KEGG4_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
dim(clusterProfiler_KEGG4)
```

    ## [1] 8 9

``` r
head(clusterProfiler_KEGG4[, -8])
```

    ##                ID                         Description GeneRatio  BgRatio
    ## hsa03010 hsa03010                            Ribosome    26/112 158/8106
    ## hsa05171 hsa05171      Coronavirus disease - COVID-19    26/112 232/8106
    ## hsa04612 hsa04612 Antigen processing and presentation     8/112  78/8106
    ## hsa05330 hsa05330                 Allograft rejection     5/112  38/8106
    ## hsa05332 hsa05332           Graft-versus-host disease     5/112  42/8106
    ## hsa04940 hsa04940            Type I diabetes mellitus     5/112  43/8106
    ##                pvalue     p.adjust       qvalue Count
    ## hsa03010 2.083064e-21 4.186959e-19 3.859151e-19    26
    ## hsa05171 4.295167e-17 4.316643e-15 3.978681e-15    26
    ## hsa04612 1.084559e-05 7.266545e-04 6.697627e-04     8
    ## hsa05330 1.606510e-04 8.072715e-03 7.440680e-03     5
    ## hsa05332 2.605642e-04 9.769639e-03 9.004747e-03     5
    ## hsa04940 2.916310e-04 9.769639e-03 9.004747e-03     5

``` r
dotplot(clusterProfiler_KEGG4, showCategory=9)
```

![](KEGG4_files/figure-markdown_github/unnamed-chunk-9-1.png)

## Compare KEGG enrichment result obtained from DAVID and clusterProfiler

``` r
merge_result(list(david=david_KEGG, clusterProfiler=clusterProfiler_KEGG, clusterProfiler4=clusterProfiler_KEGG4)) %>%
    dotplot(., showCategory=10)
```

![](KEGG4_files/figure-markdown_github/unnamed-chunk-10-1.png)

-   DAVID’s KEGG annotate 5085 genes in background while clusterProfiler
    use latest online version that annotate 6895 genes.
-   DAVID only annotate 83 genes, while clusterProfiler can annotate 104
    genes of the gene list with 207 genes in total.
-   DAVID enrich 1 KEGG terms, while clusterProfiler enrich 9 KEGG
    terms.
-   All enriched terms reported in DAVID were also reported by
    clusterProfiler.

## Session info

``` r
date()
```

    ## [1] "Thu May 13 20:36:11 2021"

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=Chinese (Simplified)_China.936 
    ## [2] LC_CTYPE=Chinese (Simplified)_China.936   
    ## [3] LC_MONETARY=Chinese (Simplified)_China.936
    ## [4] LC_NUMERIC=C                              
    ## [5] LC_TIME=Chinese (Simplified)_China.936    
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] clusterProfiler_3.19.1  RDAVIDWebService_1.28.0 ggplot2_3.3.2          
    ##  [4] GOstats_2.56.0          Category_2.56.0         Matrix_1.2-18          
    ##  [7] AnnotationDbi_1.52.0    IRanges_2.24.0          S4Vectors_0.28.1       
    ## [10] Biobase_2.50.0          graph_1.68.0            BiocGenerics_0.36.0    
    ## [13] DOSE_3.16.0             magrittr_1.5           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-149           bitops_1.0-6           ggtree_2.5.0.992      
    ##  [4] enrichplot_1.11.2.995  bit64_4.0.5            RColorBrewer_1.1-2    
    ##  [7] httr_1.4.2             Rgraphviz_2.34.0       tools_4.0.3           
    ## [10] R6_2.5.0               lazyeval_0.2.2         DBI_1.1.0             
    ## [13] colorspace_1.4-1       withr_2.3.0            tidyselect_1.1.0      
    ## [16] gridExtra_2.3          bit_4.0.4              compiler_4.0.3        
    ## [19] scatterpie_0.1.5       labeling_0.4.2         shadowtext_0.0.7      
    ## [22] scales_1.1.1           genefilter_1.72.0      RBGL_1.66.0           
    ## [25] stringr_1.4.0          digest_0.6.27          rmarkdown_2.5         
    ## [28] AnnotationForge_1.32.0 pkgconfig_2.0.3        htmltools_0.5.0       
    ## [31] rlang_0.4.8            RSQLite_2.2.1          generics_0.1.0        
    ## [34] farver_2.0.3           jsonlite_1.7.1         BiocParallel_1.24.0   
    ## [37] GOSemSim_2.17.1        dplyr_1.0.2            RCurl_1.98-1.2        
    ## [40] GO.db_3.12.1           patchwork_1.1.1        Rcpp_1.0.5            
    ## [43] munsell_0.5.0          ape_5.4-1              viridis_0.5.1         
    ## [46] lifecycle_0.2.0        stringi_1.5.3          yaml_2.2.1            
    ## [49] ggraph_2.0.3           MASS_7.3-53            plyr_1.8.6            
    ## [52] qvalue_2.22.0          grid_4.0.3             blob_1.2.1            
    ## [55] ggrepel_0.9.0          DO.db_2.9              crayon_1.3.4          
    ## [58] lattice_0.20-41        cowplot_1.1.0          graphlayouts_0.7.1    
    ## [61] splines_4.0.3          annotate_1.68.0        knitr_1.30            
    ## [64] pillar_1.4.6           fgsea_1.16.0           igraph_1.2.6          
    ## [67] reshape2_1.4.4         fastmatch_1.1-0        XML_3.99-0.5          
    ## [70] glue_1.4.2             evaluate_0.14          downloader_0.4        
    ## [73] BiocManager_1.30.10    data.table_1.13.2      treeio_1.14.0         
    ## [76] vctrs_0.3.4            tweenr_1.0.1           gtable_0.3.0          
    ## [79] purrr_0.3.4            polyclip_1.10-0        tidyr_1.1.2           
    ## [82] xfun_0.19              ggforce_0.3.2          xtable_1.8-4          
    ## [85] tidygraph_1.2.0        tidytree_0.3.4         survival_3.2-7        
    ## [88] viridisLite_0.3.0      tibble_3.0.4           rvcheck_0.1.8         
    ## [91] aplot_0.0.6            rJava_0.9-13           memoise_1.1.0         
    ## [94] ellipsis_0.3.1         ROCR_1.0-11            GSEABase_1.52.0
