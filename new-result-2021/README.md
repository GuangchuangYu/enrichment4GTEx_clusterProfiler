# functional enrichment for GTEx paper

### Guangchuang Yu

#### 05/17, 2021


The original result can be found on <https://github.com/GuangchuangYu/enrichment4GTEx_clusterProfiler>. 

This folder contains new results using DAVID (v6.8) and clusterProfiler (v4) with [Gene Ontology](https://github.com/GuangchuangYu/enrichment4GTEx_clusterProfiler/blob/master/new-result-2021/GO-BP-2021.md) and [KEGG](https://github.com/GuangchuangYu/enrichment4GTEx_clusterProfiler/blob/master/new-result-2021/KEGG-2021.md) analyses.


The results are consistent with [previous results](https://github.com/GuangchuangYu/enrichment4GTEx_clusterProfiler). clusterProfiler uses update-to-date background annotations and is able to uncover more biological processes and pathways. 


Run the following command to generate the markdown files.

```r
render("KEGG-2021.Rmd", md_document(variant='gfm'))
render("GO-BP-2021.Rmd", md_document(variant='gfm'))
```
