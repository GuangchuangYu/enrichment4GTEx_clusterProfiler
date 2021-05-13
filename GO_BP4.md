## GO (BP) enrichment with DAVID

``` r
require(magrittr)
require(DOSE)
require(RDAVIDWebService)
require(clusterProfiler)
```

``` r
load("cache/gene.rda")
david_bp4 <- enrichDAVID(gene, idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_ALL", david.user="gcyu@connect.hku.hk")
save(david_bp4, file = "david_bp4.rda")
```

``` r
summary(david_bp)[, -8]
```

    ## Warning in summary(david_bp): summary method to convert the object to data.frame
    ## is deprecated, please use as.data.frame instead.

    ##                    ID                                      Description
    ## GO:0006414 GO:0006414                         translational elongation
    ## GO:0006412 GO:0006412                                      translation
    ## GO:0010467 GO:0010467                                  gene expression
    ## GO:0006396 GO:0006396                                   RNA processing
    ## GO:0044085 GO:0044085                    cellular component biogenesis
    ## GO:0022613 GO:0022613             ribonucleoprotein complex biogenesis
    ## GO:0044267 GO:0044267               cellular protein metabolic process
    ## GO:0042254 GO:0042254                              ribosome biogenesis
    ## GO:0016070 GO:0016070                            RNA metabolic process
    ## GO:0044260 GO:0044260         cellular macromolecule metabolic process
    ## GO:0034645 GO:0034645      cellular macromolecule biosynthetic process
    ## GO:0009059 GO:0009059               macromolecule biosynthetic process
    ## GO:0006397 GO:0006397                                  mRNA processing
    ## GO:0032268 GO:0032268 regulation of cellular protein metabolic process
    ##            GeneRatio    BgRatio       pvalue     p.adjust       qvalue Count
    ## GO:0006414    25/168  101/14116 4.833921e-25 6.298599e-22 9.718725e-23    25
    ## GO:0006412    30/168  331/14116 1.916789e-17 1.248788e-14 1.926877e-15    30
    ## GO:0010467    66/168 2999/14116 1.313901e-07 5.706548e-05 8.805443e-06    66
    ## GO:0006396    20/168  547/14116 2.561454e-05 8.309328e-03 1.287468e-03    20
    ## GO:0044085    28/168 1001/14116 4.549872e-05 1.178722e-02 1.829527e-03    28
    ## GO:0022613    11/168  180/14116 5.566722e-05 1.201662e-02 1.844751e-03    11
    ## GO:0044267    49/168 2355/14116 6.422825e-05 1.188482e-02 1.844751e-03    49
    ## GO:0042254     9/168  122/14116 9.956471e-05 1.608661e-02 2.434515e-03     9
    ## GO:0016070    26/168  938/14116 1.089796e-04 1.565485e-02 2.434515e-03    26
    ## GO:0044260    86/168 5214/14116 1.478656e-04 1.908386e-02 2.972876e-03    86
    ## GO:0034645    54/168 2812/14116 1.896884e-04 2.222098e-02 3.467032e-03    54
    ## GO:0009059    54/168 2832/14116 2.282688e-04 2.448429e-02 3.824503e-03    54
    ## GO:0006397    13/168  321/14116 4.392293e-04 4.307854e-02 6.792939e-03    13
    ## GO:0032268    16/168  474/14116 5.095183e-04 4.632629e-02 7.317143e-03    16

``` r
dotplot(david_bp, showCategory=12)
```

<img src="GO_BP4_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

## GO (BP) enrichment with clusterProfiler

``` r
dim(summary(clusterProfiler_bp))
```

    ## Warning in summary(clusterProfiler_bp): summary method to convert the object to
    ## data.frame is deprecated, please use as.data.frame instead.

    ## [1] 222   9

``` r
head(summary(clusterProfiler_bp)[, -8])
```

    ## Warning in summary(clusterProfiler_bp): summary method to convert the object to
    ## data.frame is deprecated, please use as.data.frame instead.

    ##                    ID
    ## GO:0006614 GO:0006614
    ## GO:0006613 GO:0006613
    ## GO:0045047 GO:0045047
    ## GO:0019080 GO:0019080
    ## GO:0072599 GO:0072599
    ## GO:0044033 GO:0044033
    ##                                                               Description
    ## GO:0006614    SRP-dependent cotranslational protein targeting to membrane
    ## GO:0006613                  cotranslational protein targeting to membrane
    ## GO:0045047                                        protein targeting to ER
    ## GO:0019080                                          viral gene expression
    ## GO:0072599 establishment of protein localization to endoplasmic reticulum
    ## GO:0044033                               multi-organism metabolic process
    ##            GeneRatio   BgRatio       pvalue     p.adjust       qvalue Count
    ## GO:0006614    24/194 108/18585 1.975510e-25 3.923344e-22 3.349772e-22    24
    ## GO:0006613    24/194 110/18585 3.183361e-25 3.923344e-22 3.349772e-22    24
    ## GO:0045047    24/194 112/18585 5.078427e-25 3.923344e-22 3.349772e-22    24
    ## GO:0019080    29/194 199/18585 5.520006e-25 3.923344e-22 3.349772e-22    29
    ## GO:0072599    24/194 116/18585 1.256069e-24 7.142007e-22 6.097883e-22    24
    ## GO:0044033    29/194 209/18585 2.336069e-24 8.598730e-22 7.341641e-22    29

``` r
dotplot(clusterProfiler_bp, showCategory=20)
```

![](GO_BP4_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
eg=bitr(gene, "ENSEMBL", "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]
clusterProfiler_bp4 <- enrichGO(eg, ont="BP", OrgDb = org.Hs.eg.db)
```

``` r
dim(clusterProfiler_bp4)
```

    ## [1] 85  9

``` r
head(clusterProfiler_bp4[, -8])
```

    ##                    ID
    ## GO:0006614 GO:0006614
    ## GO:0006413 GO:0006413
    ## GO:0019080 GO:0019080
    ## GO:0006613 GO:0006613
    ## GO:0019083 GO:0019083
    ## GO:0000184 GO:0000184
    ##                                                                    Description
    ## GO:0006614         SRP-dependent cotranslational protein targeting to membrane
    ## GO:0006413                                            translational initiation
    ## GO:0019080                                               viral gene expression
    ## GO:0006613                       cotranslational protein targeting to membrane
    ## GO:0019083                                                 viral transcription
    ## GO:0000184 nuclear-transcribed mRNA catabolic process, nonsense-mediated decay
    ##            GeneRatio   BgRatio       pvalue     p.adjust       qvalue Count
    ## GO:0006614    23/193 105/18866 1.883501e-24 2.680929e-21 2.436185e-21    23
    ## GO:0006413    28/193 192/18866 2.155120e-24 2.680929e-21 2.436185e-21    28
    ## GO:0019080    28/193 195/18866 3.352869e-24 2.680929e-21 2.436185e-21    28
    ## GO:0006613    23/193 109/18866 4.765420e-24 2.680929e-21 2.436185e-21    23
    ## GO:0019083    27/193 178/18866 5.091016e-24 2.680929e-21 2.436185e-21    27
    ## GO:0000184    23/193 120/18866 5.052559e-23 1.900484e-20 1.726988e-20    23

``` r
dotplot(clusterProfiler_bp4, showCategory=20)
```

![](GO_BP4_files/figure-markdown_github/unnamed-chunk-9-1.png)

## Compare GO (BP) enrichment result obtained from DAVID and clusterProfiler

``` r
merge_result(list(david=david_bp, clusterProfiler=clusterProfiler_bp, clusterProfiler4=clusterProfiler_bp4)) %>%
    dotplot(., showCategory=10)
```

![](GO_BP4_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
summary(david_bp)[, "ID"] %in% summary(clusterProfiler_bp)[, "ID"]
```

    ## Warning in summary(david_bp): summary method to convert the object to data.frame
    ## is deprecated, please use as.data.frame instead.

    ## Warning in summary(clusterProfiler_bp): summary method to convert the object to
    ## data.frame is deprecated, please use as.data.frame instead.

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

-   DAVID only annotate 168 genes, while clusterProfiler can annotate
    194 genes of the gene list with 207 genes in total.
-   DAVID enrich 14 BP terms, while clusterProfiler enrich 222 BP terms.
-   All enriched terms reported in DAVID were also reported by
    clusterProfiler.

## Session info

``` r
date()
```

    ## [1] "Thu May 13 21:21:54 2021"

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
