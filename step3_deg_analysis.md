---
title: "Step 3: differential Expression analysis"
date: "2024-02-11"
author: Yongjin Park
bibliography: "MS_Ref.bib"
output:
  html_document:
    toc: true
    keep_md: true
    self_contained: true
---




```r
.data <- fileset.list("result/step1/final_matrix")
```


```r
.features <- readLines(.data$row)
.hash <- .features[str_detect(.features, "Hash")]
.hash.hdr <- "result/step1/hash"
.hash.data <- fileset.list(.hash.hdr)
.hash.info <- read.hash(.hash.data)
```

## 1. MS vs. HC at the major cell type level


```r
final.cell.type <- 
     fread("Tab/step2_cell_type.txt.gz") %>%
     left_join(.hash.info) %>%
     na.omit()
```



### a. Prepare data to just include cells annotated for their cell types and phenotypes


```r
.mkdir("result/step3/deg/")
.deg.data <- fileset.list("result/step3/deg/hc_ms")
if.needed(.deg.data, {
    .deg.data <-
        rcpp_mmutil_copy_selected_columns(.data$mtx,
                                          .data$row,
                                          .data$col,
                                          unique(final.cell.type$tag),
                                          "result/step3/deg/hc_ms")
})
```


```r
.file <- "result/step3/deg/hc_ms.rds"
if.needed(.file,
{
    ## cell types
    .celltype <-
        final.cell.type %>%
        select(tag, celltype) %>%
        as.data.frame()
    ## cell -> individual
    .cell2indv <- final.cell.type %>% 
        left_join(.hash.info) %>%
        select(tag, subject) %>%
        unique %>%
        as.data.frame()
    ## individual -> disease
    .indv2exp <- .cell2indv %>%
        select(subject) %>%
        mutate(disease = substr(`subject`, 1, 2)) %>%
        as.data.frame()

    .deg.stat <-
        make.cocoa(.deg.data, .celltype, .cell2indv, .indv2exp,
                   knn = 50, .rank = 50, .take.ln = TRUE,
                   impute.by.knn = TRUE, num.threads = 16)
    ## save for computed results
    saveRDS(.deg.stat, .file)
})
.deg.stat <- readRDS(.file)

.cts <- unique(final.cell.type$celltype)
.indvs <- unique(final.cell.type$subject)

.hc.ms.dt <-
    list(tot = sort.col(.deg.stat$sum, .cts, .indvs),
         cfa = sort.col(.deg.stat$resid.ln.mu, .cts, .indvs),
         cfa.sd = sort.col(.deg.stat$resid.ln.mu.sd, .cts, .indvs)) %>%
    combine.statistics() %>%
    na.omit() %>%
    as.data.table() %>% 
    (function(x) {
        x[, c("sample", "celltype") := tstrsplit(as.character(Var2), split="_")];
        x[, disease := substr(`sample`, 1, 2)];
        x[, gene := as.character(Var1)];
        x
    }) %>% 
    dplyr::select(-Var1, -Var2) %>% 
    as.data.table()

hc.ms.deg <-
    summarize.deg(.hc.ms.dt, tot.cutoff = 10) %>%
    parse.gene()
```


```r
plt <-
    .gg.plot(hc.ms.deg, aes(pv, fill=`celltype`)) +
    geom_histogram(bins=20, linewidth=.1, color="black") +
    xlab("p-value") +
    scale_fill_brewer(palette = "Paired")
print(plt)
```

![](Fig/STEP3/Fig_DEG_hist-1.png)<!-- -->


[PDF](Fig/STEP3//Fig_DEG_hist.pdf)

[**DOWNLOAD:** DEG MS vs HC](Tab/DEG_MS_vs_HC.txt.gz)

### b. Comparison with the bulk DEG results


```r
read.bulk <- function(.file) {
    fread(.file) %>%
        rename(hgnc_symbol = gene_name) %>%
        rename(bulk.t = t, bulk.pv = P.Value, bulk.qv = adj.P.Val) %>%
        dplyr::select(hgnc_symbol, starts_with("bulk"))
}

.file <- "data/DEG/20180513/deg.ms.hc.treg.mem.exvivo.sex_covar.ruv.txt"
bulk.mTreg.dt <- read.bulk(.file)

.file <- "data/DEG/20180513/deg.ms.hc.teff.mem.exvivo.sex_covar.ruv.txt"
bulk.mTconv.dt <- read.bulk(.file)
```




#### Memory Treg genes



##### Correlation 0.35 with p-value = 1.08e-59

![](Fig/STEP3/Fig_DEG_comparison_mTreg-1.png)<!-- -->

[PDF](Fig/STEP3//Fig_DEG_comparison_mTreg.pdf)

* Here, we computed correlation between the bulk and scRNA-seq results after selecting on the DEGs marginally significant in the bulk analysis ($|t| > 2$).

* Two numbers below the marked genes: first, we show p-value of the bulk analysis; second, we show Bonferroni-adjusted p-values of the scRNA-seq analysis.

#### Memory Tconv genes


```r
.sc <- hc.ms.deg[celltype == "mTconv"]

.bulk <- bulk.mTconv.dt

.out <- plot.bulk.sc(.sc,
                     .bulk,
                     .key.genes = c("PRDM1","SGK1"),
                     fwer.cutoff = .01,
                     qv.cutoff = .2,
                     n.top = 5)
```

#### Correlation 0.19 with p-value = 6.79e-33

![](Fig/STEP3/Fig_DEG_comparison_mTconv-1.png)<!-- -->

[PDF](Fig/STEP3//Fig_DEG_comparison_mTconv.pdf)

* Here, we computed correlation between the bulk and scRNA-seq results after selecting on the DEGs marginally significant in the bulk analysis ($|t| > 2$).

* Two numbers below the marked genes: first, we show p-value of the bulk analysis; second, we show Bonferroni-adjusted p-values of the scRNA-seq analysis.

### Found 349 unique genes strongly perturbed by MS with FWER 5%

* Up-regulated: 163

* Down-regulated:  277

* Total pairs of genes and cell types: 36,565


```r
count.deg <- function(.dt, fwer.cutoff = .01) {
    .dt[fwer < fwer.cutoff &
        sign(ADD) == sign(ADE) &
        sign(ADC) == sign(ADE),
        .(n = .N),
        by = .(celltype,
               direction = if_else(z > 0, "up", "down"))
        ] %>%
        mutate(direction = factor(direction, c("up", "down"))) %>%
        group_by(celltype) %>%
        arrange(desc(direction)) %>%
        mutate(nc = cumsum(n)) %>%
        ungroup
}
```


```r
.df <- count.deg(hc.ms.deg)

plt <-
    .gg.plot(.df, aes(x = celltype, y = n, group = direction, fill = direction)) +
    theme(legend.position = c(1,1)) +
    theme(legend.justification = c(1,1)) +
    xlab("cell type") + ylab("# DEG") +
    geom_bar(stat="identity") +
    geom_text(aes(y=nc, label=n)) +
    scale_fill_manual(values = c("#ef8a62","#67a9cf"))

print(plt)
```

![](Fig/STEP3/Fig_DEG_count-1.png)<!-- -->


[PDF](Fig/STEP3//Fig_DEG_count.pdf)

## 2. Show genes with FWER $< 10^{-2}$ overlapping with the bulk DEGs




```r
.bulk.genes <-
    rbind(bulk.mTreg.dt[bulk.qv < .2, .(hgnc_symbol)],
          bulk.mTconv.dt[bulk.qv < .2, .(hgnc_symbol)]) %>%
    unique() %>%
    unlist()

hc.ms.deg[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`gene`, split="_")]

.genes.show <-
    hc.ms.deg[fwer < 0.05 &
              hgnc_symbol %in% .bulk.genes &
              sign(ADD) == sign(ADE) &
              sign(ADC) == sign(ADE)] %>%
    select(gene) %>%
    unique %>%
    unlist %>%
    as.character
```





<table class=" lightable-paper lightable-striped" style="font-family: Helvetica; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> celltype </th>
   <th style="text-align:left;"> hgnc_symbol </th>
   <th style="text-align:left;"> z </th>
   <th style="text-align:left;"> pv </th>
   <th style="text-align:left;"> fwer </th>
   <th style="text-align:left;"> .link </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="47"> mTconv </td>
   <td style="text-align:left;"> CD27 </td>
   <td style="text-align:left;"> -43.51 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD27_CD27.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV38-2DV8 </td>
   <td style="text-align:left;"> -7.74 </td>
   <td style="text-align:left;"> 1.02e-14 </td>
   <td style="text-align:left;"> 3.72e-10 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211817_TRAV38-2DV8.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PNP </td>
   <td style="text-align:left;"> -4.89 </td>
   <td style="text-align:left;"> 1.03e-06 </td>
   <td style="text-align:left;"> 3.68e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198805_PNP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> VTA1 </td>
   <td style="text-align:left;"> -5.71 </td>
   <td style="text-align:left;"> 1.14e-08 </td>
   <td style="text-align:left;"> 4.12e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000009844_VTA1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> P2RY8 </td>
   <td style="text-align:left;"> -5.29 </td>
   <td style="text-align:left;"> 1.20e-07 </td>
   <td style="text-align:left;"> 4.34e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182162_P2RY8.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV24-1 </td>
   <td style="text-align:left;"> 10.02 </td>
   <td style="text-align:left;"> 1.28e-23 </td>
   <td style="text-align:left;"> 4.65e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211750_TRBV24-1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PASK </td>
   <td style="text-align:left;"> 4.83 </td>
   <td style="text-align:left;"> 1.33e-06 </td>
   <td style="text-align:left;"> 4.79e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000115687_PASK.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD226 </td>
   <td style="text-align:left;"> 12.99 </td>
   <td style="text-align:left;"> 1.37e-38 </td>
   <td style="text-align:left;"> 5.00e-34 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD226_CD226.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DENND2D </td>
   <td style="text-align:left;"> -10 </td>
   <td style="text-align:left;"> 1.46e-23 </td>
   <td style="text-align:left;"> 5.31e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162777_DENND2D.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> KLRG1 </td>
   <td style="text-align:left;"> -6.74 </td>
   <td style="text-align:left;"> 1.54e-11 </td>
   <td style="text-align:left;"> 5.59e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000139187_KLRG1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CITED2 </td>
   <td style="text-align:left;"> 7.68 </td>
   <td style="text-align:left;"> 1.57e-14 </td>
   <td style="text-align:left;"> 5.72e-10 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000164442_CITED2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP2 </td>
   <td style="text-align:left;"> -5.63 </td>
   <td style="text-align:left;"> 1.80e-08 </td>
   <td style="text-align:left;"> 6.50e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000106560_GIMAP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DYNLL1 </td>
   <td style="text-align:left;"> -6.01 </td>
   <td style="text-align:left;"> 1.80e-09 </td>
   <td style="text-align:left;"> 6.52e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000088986_DYNLL1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MRPS6 </td>
   <td style="text-align:left;"> -5.21 </td>
   <td style="text-align:left;"> 1.91e-07 </td>
   <td style="text-align:left;"> 6.89e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000243927_MRPS6.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV13-1 </td>
   <td style="text-align:left;"> -7.35 </td>
   <td style="text-align:left;"> 1.97e-13 </td>
   <td style="text-align:left;"> 7.16e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211788_TRAV13-1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TC2N </td>
   <td style="text-align:left;"> -5.99 </td>
   <td style="text-align:left;"> 2.13e-09 </td>
   <td style="text-align:left;"> 7.72e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000165929_TC2N.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA2 </td>
   <td style="text-align:left;"> 5.57 </td>
   <td style="text-align:left;"> 2.56e-08 </td>
   <td style="text-align:left;"> 9.25e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182718_ANXA2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DDIT4 </td>
   <td style="text-align:left;"> 10.61 </td>
   <td style="text-align:left;"> 2.64e-26 </td>
   <td style="text-align:left;"> 9.62e-22 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000168209_DDIT4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PRDM1 </td>
   <td style="text-align:left;"> 5.94 </td>
   <td style="text-align:left;"> 2.92e-09 </td>
   <td style="text-align:left;"> 1.06e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000057657_PRDM1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD40LG </td>
   <td style="text-align:left;"> -5.11 </td>
   <td style="text-align:left;"> 3.27e-07 </td>
   <td style="text-align:left;"> 1.18e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000102245_CD40LG.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PLP2 </td>
   <td style="text-align:left;"> 5.52 </td>
   <td style="text-align:left;"> 3.46e-08 </td>
   <td style="text-align:left;"> 1.25e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000102007_PLP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PHC3 </td>
   <td style="text-align:left;"> -5.91 </td>
   <td style="text-align:left;"> 3.49e-09 </td>
   <td style="text-align:left;"> 1.26e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000173889_PHC3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> HERPUD1 </td>
   <td style="text-align:left;"> 8.43 </td>
   <td style="text-align:left;"> 3.49e-17 </td>
   <td style="text-align:left;"> 1.27e-12 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000051108_HERPUD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAC </td>
   <td style="text-align:left;"> -6.62 </td>
   <td style="text-align:left;"> 3.53e-11 </td>
   <td style="text-align:left;"> 1.28e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000277734_TRAC.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ID3 </td>
   <td style="text-align:left;"> -5.05 </td>
   <td style="text-align:left;"> 4.51e-07 </td>
   <td style="text-align:left;"> 1.62e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000117318_ID3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DHX36 </td>
   <td style="text-align:left;"> -5.45 </td>
   <td style="text-align:left;"> 4.93e-08 </td>
   <td style="text-align:left;"> 1.78e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000174953_DHX36.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA1 </td>
   <td style="text-align:left;"> -7.82 </td>
   <td style="text-align:left;"> 5.40e-15 </td>
   <td style="text-align:left;"> 1.96e-10 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000135046_ANXA1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> XPC </td>
   <td style="text-align:left;"> -5.01 </td>
   <td style="text-align:left;"> 5.49e-07 </td>
   <td style="text-align:left;"> 1.97e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000154767_XPC.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ARPC5 </td>
   <td style="text-align:left;"> -8.64 </td>
   <td style="text-align:left;"> 5.58e-18 </td>
   <td style="text-align:left;"> 2.03e-13 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162704_ARPC5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM128 </td>
   <td style="text-align:left;"> -5 </td>
   <td style="text-align:left;"> 5.65e-07 </td>
   <td style="text-align:left;"> 2.03e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000132406_TMEM128.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> FTH1 </td>
   <td style="text-align:left;"> 10.1 </td>
   <td style="text-align:left;"> 5.74e-24 </td>
   <td style="text-align:left;"> 2.09e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000167996_FTH1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CHST12 </td>
   <td style="text-align:left;"> -5.42 </td>
   <td style="text-align:left;"> 6.02e-08 </td>
   <td style="text-align:left;"> 2.17e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000136213_CHST12.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV11-2 </td>
   <td style="text-align:left;"> -8.36 </td>
   <td style="text-align:left;"> 6.04e-17 </td>
   <td style="text-align:left;"> 2.20e-12 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000241657_TRBV11-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> RNF19A </td>
   <td style="text-align:left;"> 4.98 </td>
   <td style="text-align:left;"> 6.21e-07 </td>
   <td style="text-align:left;"> 2.23e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000034677_RNF19A.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ICAM2 </td>
   <td style="text-align:left;"> -8.08 </td>
   <td style="text-align:left;"> 6.41e-16 </td>
   <td style="text-align:left;"> 2.33e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000108622_ICAM2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PSMA5 </td>
   <td style="text-align:left;"> -4.97 </td>
   <td style="text-align:left;"> 6.55e-07 </td>
   <td style="text-align:left;"> 2.35e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000143106_PSMA5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV17 </td>
   <td style="text-align:left;"> -6.17 </td>
   <td style="text-align:left;"> 6.77e-10 </td>
   <td style="text-align:left;"> 2.45e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211797_TRAV17.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAT1 </td>
   <td style="text-align:left;"> -4.97 </td>
   <td style="text-align:left;"> 6.82e-07 </td>
   <td style="text-align:left;"> 2.45e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163519_TRAT1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> OSM </td>
   <td style="text-align:left;"> 7.49 </td>
   <td style="text-align:left;"> 6.91e-14 </td>
   <td style="text-align:left;"> 2.51e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000099985_OSM.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV29DV5 </td>
   <td style="text-align:left;"> -6.51 </td>
   <td style="text-align:left;"> 7.35e-11 </td>
   <td style="text-align:left;"> 2.66e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211810_TRAV29DV5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV21 </td>
   <td style="text-align:left;"> -6.84 </td>
   <td style="text-align:left;"> 7.70e-12 </td>
   <td style="text-align:left;"> 2.80e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211801_TRAV21.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PSMB8-AS1 </td>
   <td style="text-align:left;"> -6.15 </td>
   <td style="text-align:left;"> 7.90e-10 </td>
   <td style="text-align:left;"> 2.86e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000204261_PSMB8-AS1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PPWD1 </td>
   <td style="text-align:left;"> -5.76 </td>
   <td style="text-align:left;"> 8.45e-09 </td>
   <td style="text-align:left;"> 3.06e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000113593_PPWD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MRPL18 </td>
   <td style="text-align:left;"> -5.76 </td>
   <td style="text-align:left;"> 8.55e-09 </td>
   <td style="text-align:left;"> 3.09e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000112110_MRPL18.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP6 </td>
   <td style="text-align:left;"> -6.83 </td>
   <td style="text-align:left;"> 8.55e-12 </td>
   <td style="text-align:left;"> 3.10e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000133561_GIMAP6.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> XRCC5 </td>
   <td style="text-align:left;"> -4.9 </td>
   <td style="text-align:left;"> 9.43e-07 </td>
   <td style="text-align:left;"> 3.39e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000079246_XRCC5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV8-2 </td>
   <td style="text-align:left;"> -5.73 </td>
   <td style="text-align:left;"> 9.95e-09 </td>
   <td style="text-align:left;"> 3.60e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211786_TRAV8-2.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="69"> mTreg </td>
   <td style="text-align:left;"> CD226 </td>
   <td style="text-align:left;"> 53.1 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD226_CD226.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD27 </td>
   <td style="text-align:left;"> -106.42 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD27_CD27.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV30 </td>
   <td style="text-align:left;"> 5.72 </td>
   <td style="text-align:left;"> 1.05e-08 </td>
   <td style="text-align:left;"> 3.79e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000259092_TRAV30.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP1 </td>
   <td style="text-align:left;"> -7.42 </td>
   <td style="text-align:left;"> 1.14e-13 </td>
   <td style="text-align:left;"> 4.14e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000213203_GIMAP1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> RSBN1L </td>
   <td style="text-align:left;"> -6.44 </td>
   <td style="text-align:left;"> 1.23e-10 </td>
   <td style="text-align:left;"> 4.45e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000187257_RSBN1L.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PRDX3 </td>
   <td style="text-align:left;"> -8.28 </td>
   <td style="text-align:left;"> 1.23e-16 </td>
   <td style="text-align:left;"> 4.48e-12 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000165672_PRDX3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP6 </td>
   <td style="text-align:left;"> -6.42 </td>
   <td style="text-align:left;"> 1.32e-10 </td>
   <td style="text-align:left;"> 4.79e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000133561_GIMAP6.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DUSP4 </td>
   <td style="text-align:left;"> 10.01 </td>
   <td style="text-align:left;"> 1.33e-23 </td>
   <td style="text-align:left;"> 4.84e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000120875_DUSP4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> C14orf119 </td>
   <td style="text-align:left;"> -5.26 </td>
   <td style="text-align:left;"> 1.41e-07 </td>
   <td style="text-align:left;"> 5.08e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000179933_C14orf119.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ODC1 </td>
   <td style="text-align:left;"> 5.67 </td>
   <td style="text-align:left;"> 1.44e-08 </td>
   <td style="text-align:left;"> 5.22e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000115758_ODC1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TUBA1B </td>
   <td style="text-align:left;"> -7.39 </td>
   <td style="text-align:left;"> 1.44e-13 </td>
   <td style="text-align:left;"> 5.24e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000123416_TUBA1B.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV13-1 </td>
   <td style="text-align:left;"> -6.74 </td>
   <td style="text-align:left;"> 1.60e-11 </td>
   <td style="text-align:left;"> 5.79e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211788_TRAV13-1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GADD45A </td>
   <td style="text-align:left;"> 5.23 </td>
   <td style="text-align:left;"> 1.65e-07 </td>
   <td style="text-align:left;"> 5.95e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000116717_GADD45A.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PRDM1 </td>
   <td style="text-align:left;"> 7.06 </td>
   <td style="text-align:left;"> 1.71e-12 </td>
   <td style="text-align:left;"> 6.23e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000057657_PRDM1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PSMA5 </td>
   <td style="text-align:left;"> -5.64 </td>
   <td style="text-align:left;"> 1.74e-08 </td>
   <td style="text-align:left;"> 6.30e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000143106_PSMA5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TNFAIP8L2 </td>
   <td style="text-align:left;"> -7.67 </td>
   <td style="text-align:left;"> 1.74e-14 </td>
   <td style="text-align:left;"> 6.33e-10 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163154_TNFAIP8L2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP8 </td>
   <td style="text-align:left;"> -5.22 </td>
   <td style="text-align:left;"> 1.76e-07 </td>
   <td style="text-align:left;"> 6.33e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000171115_GIMAP8.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LINC00402 </td>
   <td style="text-align:left;"> -6.72 </td>
   <td style="text-align:left;"> 1.78e-11 </td>
   <td style="text-align:left;"> 6.44e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000235532_LINC00402.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV12-3 </td>
   <td style="text-align:left;"> -6.01 </td>
   <td style="text-align:left;"> 1.83e-09 </td>
   <td style="text-align:left;"> 6.63e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000274752_TRBV12-3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA2R </td>
   <td style="text-align:left;"> -5.63 </td>
   <td style="text-align:left;"> 1.84e-08 </td>
   <td style="text-align:left;"> 6.64e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000177721_ANXA2R.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SIT1 </td>
   <td style="text-align:left;"> -6.01 </td>
   <td style="text-align:left;"> 1.84e-09 </td>
   <td style="text-align:left;"> 6.67e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000137078_SIT1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ARPC5 </td>
   <td style="text-align:left;"> -8.75 </td>
   <td style="text-align:left;"> 2.06e-18 </td>
   <td style="text-align:left;"> 7.51e-14 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162704_ARPC5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SLC25A20 </td>
   <td style="text-align:left;"> -6.36 </td>
   <td style="text-align:left;"> 2.08e-10 </td>
   <td style="text-align:left;"> 7.54e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000178537_SLC25A20.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SATB1 </td>
   <td style="text-align:left;"> -5.6 </td>
   <td style="text-align:left;"> 2.09e-08 </td>
   <td style="text-align:left;"> 7.53e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182568_SATB1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV12-2 </td>
   <td style="text-align:left;"> 10.84 </td>
   <td style="text-align:left;"> 2.14e-27 </td>
   <td style="text-align:left;"> 7.80e-23 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211789_TRAV12-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD58 </td>
   <td style="text-align:left;"> 19.34 </td>
   <td style="text-align:left;"> 2.25e-83 </td>
   <td style="text-align:left;"> 8.24e-79 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD58_CD58.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MRPS6 </td>
   <td style="text-align:left;"> -6.67 </td>
   <td style="text-align:left;"> 2.54e-11 </td>
   <td style="text-align:left;"> 9.23e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000243927_MRPS6.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> BCL3 </td>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> 2.60e-12 </td>
   <td style="text-align:left;"> 9.46e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000069399_BCL3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LY9 </td>
   <td style="text-align:left;"> -5.94 </td>
   <td style="text-align:left;"> 2.86e-09 </td>
   <td style="text-align:left;"> 1.03e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000122224_LY9.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LIMS1 </td>
   <td style="text-align:left;"> -8.18 </td>
   <td style="text-align:left;"> 2.86e-16 </td>
   <td style="text-align:left;"> 1.04e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000169756_LIMS1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA1 </td>
   <td style="text-align:left;"> -5.13 </td>
   <td style="text-align:left;"> 2.96e-07 </td>
   <td style="text-align:left;"> 1.06e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000135046_ANXA1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> M6PR </td>
   <td style="text-align:left;"> -6.3 </td>
   <td style="text-align:left;"> 3.03e-10 </td>
   <td style="text-align:left;"> 1.10e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000003056_M6PR.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DDIT4 </td>
   <td style="text-align:left;"> 12.57 </td>
   <td style="text-align:left;"> 3.15e-36 </td>
   <td style="text-align:left;"> 1.15e-31 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000168209_DDIT4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> COPB2 </td>
   <td style="text-align:left;"> -5.53 </td>
   <td style="text-align:left;"> 3.17e-08 </td>
   <td style="text-align:left;"> 1.15e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000184432_COPB2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> IRF1 </td>
   <td style="text-align:left;"> -5.52 </td>
   <td style="text-align:left;"> 3.34e-08 </td>
   <td style="text-align:left;"> 1.21e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000125347_IRF1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> IMPDH2 </td>
   <td style="text-align:left;"> -5.52 </td>
   <td style="text-align:left;"> 3.36e-08 </td>
   <td style="text-align:left;"> 1.21e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000178035_IMPDH2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PTPRCAP </td>
   <td style="text-align:left;"> 6.63 </td>
   <td style="text-align:left;"> 3.46e-11 </td>
   <td style="text-align:left;"> 1.25e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000213402_PTPRCAP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV29DV5 </td>
   <td style="text-align:left;"> -9.69 </td>
   <td style="text-align:left;"> 3.47e-22 </td>
   <td style="text-align:left;"> 1.27e-17 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211810_TRAV29DV5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GRN </td>
   <td style="text-align:left;"> 6.27 </td>
   <td style="text-align:left;"> 3.59e-10 </td>
   <td style="text-align:left;"> 1.30e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000030582_GRN.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV24-1 </td>
   <td style="text-align:left;"> -12 </td>
   <td style="text-align:left;"> 3.77e-33 </td>
   <td style="text-align:left;"> 1.37e-28 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211750_TRBV24-1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> FTH1 </td>
   <td style="text-align:left;"> 18.22 </td>
   <td style="text-align:left;"> 3.86e-74 </td>
   <td style="text-align:left;"> 1.41e-69 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000167996_FTH1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SGK1 </td>
   <td style="text-align:left;"> 5.49 </td>
   <td style="text-align:left;"> 4.03e-08 </td>
   <td style="text-align:left;"> 1.45e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000118515_SGK1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DYNLL1 </td>
   <td style="text-align:left;"> -5.49 </td>
   <td style="text-align:left;"> 4.12e-08 </td>
   <td style="text-align:left;"> 1.49e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000088986_DYNLL1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA2 </td>
   <td style="text-align:left;"> 8.13 </td>
   <td style="text-align:left;"> 4.38e-16 </td>
   <td style="text-align:left;"> 1.59e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182718_ANXA2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> VBP1 </td>
   <td style="text-align:left;"> -5.47 </td>
   <td style="text-align:left;"> 4.56e-08 </td>
   <td style="text-align:left;"> 1.64e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000155959_VBP1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CTLA4 </td>
   <td style="text-align:left;"> 7.24 </td>
   <td style="text-align:left;"> 4.59e-13 </td>
   <td style="text-align:left;"> 1.67e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163599_CTLA4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PNP </td>
   <td style="text-align:left;"> -5.86 </td>
   <td style="text-align:left;"> 4.60e-09 </td>
   <td style="text-align:left;"> 1.66e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198805_PNP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ICAM2 </td>
   <td style="text-align:left;"> -10.11 </td>
   <td style="text-align:left;"> 4.77e-24 </td>
   <td style="text-align:left;"> 1.74e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000108622_ICAM2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> NAMPT </td>
   <td style="text-align:left;"> 5.84 </td>
   <td style="text-align:left;"> 5.35e-09 </td>
   <td style="text-align:left;"> 1.93e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000105835_NAMPT.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SRGN </td>
   <td style="text-align:left;"> 6.89 </td>
   <td style="text-align:left;"> 5.74e-12 </td>
   <td style="text-align:left;"> 2.08e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000122862_SRGN.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM128 </td>
   <td style="text-align:left;"> -4.99 </td>
   <td style="text-align:left;"> 5.94e-07 </td>
   <td style="text-align:left;"> 2.14e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000132406_TMEM128.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRGV10 </td>
   <td style="text-align:left;"> 5.42 </td>
   <td style="text-align:left;"> 6.00e-08 </td>
   <td style="text-align:left;"> 2.17e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211694_TRGV10.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM35B </td>
   <td style="text-align:left;"> -5.8 </td>
   <td style="text-align:left;"> 6.51e-09 </td>
   <td style="text-align:left;"> 2.36e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000243749_TMEM35B.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DENND2D </td>
   <td style="text-align:left;"> -8.62 </td>
   <td style="text-align:left;"> 6.65e-18 </td>
   <td style="text-align:left;"> 2.42e-13 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162777_DENND2D.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV38-2DV8 </td>
   <td style="text-align:left;"> -13.39 </td>
   <td style="text-align:left;"> 7.12e-41 </td>
   <td style="text-align:left;"> 2.60e-36 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211817_TRAV38-2DV8.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MT2A </td>
   <td style="text-align:left;"> 7.48 </td>
   <td style="text-align:left;"> 7.35e-14 </td>
   <td style="text-align:left;"> 2.67e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000125148_MT2A.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SDHA </td>
   <td style="text-align:left;"> -5.78 </td>
   <td style="text-align:left;"> 7.47e-09 </td>
   <td style="text-align:left;"> 2.70e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000073578_SDHA.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PSMB8-AS1 </td>
   <td style="text-align:left;"> -5.38 </td>
   <td style="text-align:left;"> 7.58e-08 </td>
   <td style="text-align:left;"> 2.73e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000204261_PSMB8-AS1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> APBB1IP </td>
   <td style="text-align:left;"> -5.37 </td>
   <td style="text-align:left;"> 7.69e-08 </td>
   <td style="text-align:left;"> 2.77e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000077420_APBB1IP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> NQO2 </td>
   <td style="text-align:left;"> -6.51 </td>
   <td style="text-align:left;"> 7.75e-11 </td>
   <td style="text-align:left;"> 2.81e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000124588_NQO2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV11-2 </td>
   <td style="text-align:left;"> 8.6 </td>
   <td style="text-align:left;"> 7.80e-18 </td>
   <td style="text-align:left;"> 2.84e-13 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000241657_TRBV11-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> HERPUD1 </td>
   <td style="text-align:left;"> 9.12 </td>
   <td style="text-align:left;"> 7.81e-20 </td>
   <td style="text-align:left;"> 2.85e-15 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000051108_HERPUD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> RGS14 </td>
   <td style="text-align:left;"> -5.77 </td>
   <td style="text-align:left;"> 7.96e-09 </td>
   <td style="text-align:left;"> 2.88e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000169220_RGS14.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM106B </td>
   <td style="text-align:left;"> -5.77 </td>
   <td style="text-align:left;"> 8.09e-09 </td>
   <td style="text-align:left;"> 2.93e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000106460_TMEM106B.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP2 </td>
   <td style="text-align:left;"> -7.16 </td>
   <td style="text-align:left;"> 8.18e-13 </td>
   <td style="text-align:left;"> 2.97e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000106560_GIMAP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TUBB </td>
   <td style="text-align:left;"> -8.05 </td>
   <td style="text-align:left;"> 8.31e-16 </td>
   <td style="text-align:left;"> 3.02e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000196230_TUBB.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PLP2 </td>
   <td style="text-align:left;"> 5.34 </td>
   <td style="text-align:left;"> 9.18e-08 </td>
   <td style="text-align:left;"> 3.31e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000102007_PLP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PELI1 </td>
   <td style="text-align:left;"> 5.74 </td>
   <td style="text-align:left;"> 9.35e-09 </td>
   <td style="text-align:left;"> 3.38e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000197329_PELI1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TIGIT </td>
   <td style="text-align:left;"> -15.58 </td>
   <td style="text-align:left;"> 9.57e-55 </td>
   <td style="text-align:left;"> 3.49e-50 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_TIGIT_TIGIT.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="32"> nTconv </td>
   <td style="text-align:left;"> CD27 </td>
   <td style="text-align:left;"> -73.83 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD27_CD27.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM35B </td>
   <td style="text-align:left;"> -4.89 </td>
   <td style="text-align:left;"> 1.03e-06 </td>
   <td style="text-align:left;"> 3.71e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000243749_TMEM35B.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> S100A11 </td>
   <td style="text-align:left;"> -5.32 </td>
   <td style="text-align:left;"> 1.04e-07 </td>
   <td style="text-align:left;"> 3.75e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163191_S100A11.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> HERPUD1 </td>
   <td style="text-align:left;"> 4.87 </td>
   <td style="text-align:left;"> 1.09e-06 </td>
   <td style="text-align:left;"> 3.91e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000051108_HERPUD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TUBB </td>
   <td style="text-align:left;"> -5.31 </td>
   <td style="text-align:left;"> 1.10e-07 </td>
   <td style="text-align:left;"> 3.97e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000196230_TUBB.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TC2N </td>
   <td style="text-align:left;"> -4.87 </td>
   <td style="text-align:left;"> 1.12e-06 </td>
   <td style="text-align:left;"> 4.02e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000165929_TC2N.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM128 </td>
   <td style="text-align:left;"> -4.84 </td>
   <td style="text-align:left;"> 1.32e-06 </td>
   <td style="text-align:left;"> 4.74e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000132406_TMEM128.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ODC1 </td>
   <td style="text-align:left;"> 4.83 </td>
   <td style="text-align:left;"> 1.36e-06 </td>
   <td style="text-align:left;"> 4.88e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000115758_ODC1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> WDR1 </td>
   <td style="text-align:left;"> -5.24 </td>
   <td style="text-align:left;"> 1.57e-07 </td>
   <td style="text-align:left;"> 5.67e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000071127_WDR1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> XBP1 </td>
   <td style="text-align:left;"> 9.99 </td>
   <td style="text-align:left;"> 1.63e-23 </td>
   <td style="text-align:left;"> 5.95e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000100219_XBP1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> M6PR </td>
   <td style="text-align:left;"> -5.64 </td>
   <td style="text-align:left;"> 1.67e-08 </td>
   <td style="text-align:left;"> 6.03e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000003056_M6PR.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ICAM2 </td>
   <td style="text-align:left;"> -6.02 </td>
   <td style="text-align:left;"> 1.70e-09 </td>
   <td style="text-align:left;"> 6.15e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000108622_ICAM2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV8-2 </td>
   <td style="text-align:left;"> -5.22 </td>
   <td style="text-align:left;"> 1.83e-07 </td>
   <td style="text-align:left;"> 6.60e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211786_TRAV8-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ACTG1 </td>
   <td style="text-align:left;"> -8.48 </td>
   <td style="text-align:left;"> 2.17e-17 </td>
   <td style="text-align:left;"> 7.90e-13 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000184009_ACTG1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PTGER2 </td>
   <td style="text-align:left;"> 5.18 </td>
   <td style="text-align:left;"> 2.18e-07 </td>
   <td style="text-align:left;"> 7.83e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000125384_PTGER2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> APBB1IP </td>
   <td style="text-align:left;"> -6.31 </td>
   <td style="text-align:left;"> 2.74e-10 </td>
   <td style="text-align:left;"> 9.92e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000077420_APBB1IP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV7-9 </td>
   <td style="text-align:left;"> -11.23 </td>
   <td style="text-align:left;"> 2.83e-29 </td>
   <td style="text-align:left;"> 1.03e-24 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000278030_TRBV7-9.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SIT1 </td>
   <td style="text-align:left;"> -5.13 </td>
   <td style="text-align:left;"> 2.93e-07 </td>
   <td style="text-align:left;"> 1.05e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000137078_SIT1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV27 </td>
   <td style="text-align:left;"> 5.08 </td>
   <td style="text-align:left;"> 3.77e-07 </td>
   <td style="text-align:left;"> 1.36e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211809_TRAV27.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP2 </td>
   <td style="text-align:left;"> -5.89 </td>
   <td style="text-align:left;"> 3.84e-09 </td>
   <td style="text-align:left;"> 1.39e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000106560_GIMAP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV11-2 </td>
   <td style="text-align:left;"> -5.89 </td>
   <td style="text-align:left;"> 3.95e-09 </td>
   <td style="text-align:left;"> 1.43e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000241657_TRBV11-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CISH </td>
   <td style="text-align:left;"> 5.47 </td>
   <td style="text-align:left;"> 4.42e-08 </td>
   <td style="text-align:left;"> 1.59e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000114737_CISH.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ARPC5 </td>
   <td style="text-align:left;"> -8.11 </td>
   <td style="text-align:left;"> 5.09e-16 </td>
   <td style="text-align:left;"> 1.85e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162704_ARPC5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAC </td>
   <td style="text-align:left;"> -6.9 </td>
   <td style="text-align:left;"> 5.27e-12 </td>
   <td style="text-align:left;"> 1.91e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000277734_TRAC.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD226 </td>
   <td style="text-align:left;"> 13.9 </td>
   <td style="text-align:left;"> 6.23e-44 </td>
   <td style="text-align:left;"> 2.27e-39 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD226_CD226.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PASK </td>
   <td style="text-align:left;"> 10.53 </td>
   <td style="text-align:left;"> 6.50e-26 </td>
   <td style="text-align:left;"> 2.37e-21 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000115687_PASK.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DDIT4 </td>
   <td style="text-align:left;"> 8.88 </td>
   <td style="text-align:left;"> 6.98e-19 </td>
   <td style="text-align:left;"> 2.54e-14 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000168209_DDIT4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA1 </td>
   <td style="text-align:left;"> -8.87 </td>
   <td style="text-align:left;"> 7.01e-19 </td>
   <td style="text-align:left;"> 2.55e-14 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000135046_ANXA1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PNP </td>
   <td style="text-align:left;"> -6.17 </td>
   <td style="text-align:left;"> 7.02e-10 </td>
   <td style="text-align:left;"> 2.54e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198805_PNP.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DENND2D </td>
   <td style="text-align:left;"> -6.85 </td>
   <td style="text-align:left;"> 7.14e-12 </td>
   <td style="text-align:left;"> 2.59e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162777_DENND2D.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> RGCC </td>
   <td style="text-align:left;"> -6.16 </td>
   <td style="text-align:left;"> 7.50e-10 </td>
   <td style="text-align:left;"> 2.72e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000102760_RGCC.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TNFAIP8L2 </td>
   <td style="text-align:left;"> -5.34 </td>
   <td style="text-align:left;"> 9.49e-08 </td>
   <td style="text-align:left;"> 3.42e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163154_TNFAIP8L2.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="35"> nTreg </td>
   <td style="text-align:left;"> CD27 </td>
   <td style="text-align:left;"> -50.46 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> 0.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD27_CD27.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SBDS </td>
   <td style="text-align:left;"> 5.3 </td>
   <td style="text-align:left;"> 1.14e-07 </td>
   <td style="text-align:left;"> 4.11e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000126524_SBDS.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LINC00402 </td>
   <td style="text-align:left;"> -4.86 </td>
   <td style="text-align:left;"> 1.19e-06 </td>
   <td style="text-align:left;"> 4.29e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000235532_LINC00402.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DUSP4 </td>
   <td style="text-align:left;"> 4.85 </td>
   <td style="text-align:left;"> 1.23e-06 </td>
   <td style="text-align:left;"> 4.42e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000120875_DUSP4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SCAMP3 </td>
   <td style="text-align:left;"> -4.84 </td>
   <td style="text-align:left;"> 1.33e-06 </td>
   <td style="text-align:left;"> 4.77e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000116521_SCAMP3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ICAM2 </td>
   <td style="text-align:left;"> -6.76 </td>
   <td style="text-align:left;"> 1.33e-11 </td>
   <td style="text-align:left;"> 4.84e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000108622_ICAM2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DENND2D </td>
   <td style="text-align:left;"> -6.06 </td>
   <td style="text-align:left;"> 1.37e-09 </td>
   <td style="text-align:left;"> 4.97e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162777_DENND2D.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PELI1 </td>
   <td style="text-align:left;"> 5.27 </td>
   <td style="text-align:left;"> 1.39e-07 </td>
   <td style="text-align:left;"> 5.02e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000197329_PELI1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CHCHD1 </td>
   <td style="text-align:left;"> -5.25 </td>
   <td style="text-align:left;"> 1.52e-07 </td>
   <td style="text-align:left;"> 5.48e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000172586_CHCHD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> M6PR </td>
   <td style="text-align:left;"> -6.03 </td>
   <td style="text-align:left;"> 1.66e-09 </td>
   <td style="text-align:left;"> 6.00e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000003056_M6PR.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MRPL18 </td>
   <td style="text-align:left;"> -5.64 </td>
   <td style="text-align:left;"> 1.67e-08 </td>
   <td style="text-align:left;"> 6.02e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000112110_MRPL18.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> HERPUD1 </td>
   <td style="text-align:left;"> 5.98 </td>
   <td style="text-align:left;"> 2.29e-09 </td>
   <td style="text-align:left;"> 8.29e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000051108_HERPUD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SOCS3 </td>
   <td style="text-align:left;"> -5.97 </td>
   <td style="text-align:left;"> 2.41e-09 </td>
   <td style="text-align:left;"> 8.70e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000184557_SOCS3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD58 </td>
   <td style="text-align:left;"> 7.6 </td>
   <td style="text-align:left;"> 2.91e-14 </td>
   <td style="text-align:left;"> 1.06e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD58_CD58.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> GIMAP2 </td>
   <td style="text-align:left;"> -8.17 </td>
   <td style="text-align:left;"> 3.21e-16 </td>
   <td style="text-align:left;"> 1.17e-11 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000106560_GIMAP2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> MAP3K8 </td>
   <td style="text-align:left;"> 5.1 </td>
   <td style="text-align:left;"> 3.38e-07 </td>
   <td style="text-align:left;"> 1.21e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000107968_MAP3K8.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TMEM223 </td>
   <td style="text-align:left;"> -5.52 </td>
   <td style="text-align:left;"> 3.43e-08 </td>
   <td style="text-align:left;"> 1.24e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000168569_TMEM223.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAT1 </td>
   <td style="text-align:left;"> -6.62 </td>
   <td style="text-align:left;"> 3.50e-11 </td>
   <td style="text-align:left;"> 1.27e-06 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000163519_TRAT1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SATB1 </td>
   <td style="text-align:left;"> -5.5 </td>
   <td style="text-align:left;"> 3.80e-08 </td>
   <td style="text-align:left;"> 1.37e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182568_SATB1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ASF1A </td>
   <td style="text-align:left;"> -5.07 </td>
   <td style="text-align:left;"> 3.88e-07 </td>
   <td style="text-align:left;"> 1.40e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000111875_ASF1A.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA2 </td>
   <td style="text-align:left;"> 6.92 </td>
   <td style="text-align:left;"> 4.52e-12 </td>
   <td style="text-align:left;"> 1.64e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000182718_ANXA2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD226 </td>
   <td style="text-align:left;"> 16.62 </td>
   <td style="text-align:left;"> 4.70e-62 </td>
   <td style="text-align:left;"> 1.72e-57 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD226_CD226.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> SIT1 </td>
   <td style="text-align:left;"> -6.21 </td>
   <td style="text-align:left;"> 5.28e-10 </td>
   <td style="text-align:left;"> 1.91e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000137078_SIT1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ZFP90 </td>
   <td style="text-align:left;"> -5.82 </td>
   <td style="text-align:left;"> 5.99e-09 </td>
   <td style="text-align:left;"> 2.17e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000184939_ZFP90.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV8-2 </td>
   <td style="text-align:left;"> -7.2 </td>
   <td style="text-align:left;"> 6.03e-13 </td>
   <td style="text-align:left;"> 2.19e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211786_TRAV8-2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> RAB29 </td>
   <td style="text-align:left;"> -5.81 </td>
   <td style="text-align:left;"> 6.13e-09 </td>
   <td style="text-align:left;"> 2.22e-04 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000117280_RAB29.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LY9 </td>
   <td style="text-align:left;"> -4.97 </td>
   <td style="text-align:left;"> 6.80e-07 </td>
   <td style="text-align:left;"> 2.44e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000122224_LY9.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ARPC5 </td>
   <td style="text-align:left;"> -7.18 </td>
   <td style="text-align:left;"> 6.98e-13 </td>
   <td style="text-align:left;"> 2.53e-08 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000162704_ARPC5.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> DDIT4 </td>
   <td style="text-align:left;"> 11.15 </td>
   <td style="text-align:left;"> 7.27e-29 </td>
   <td style="text-align:left;"> 2.65e-24 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000168209_DDIT4.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRBV7-9 </td>
   <td style="text-align:left;"> -7.48 </td>
   <td style="text-align:left;"> 7.46e-14 </td>
   <td style="text-align:left;"> 2.71e-09 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000278030_TRBV7-9.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ANXA2R </td>
   <td style="text-align:left;"> -4.94 </td>
   <td style="text-align:left;"> 7.92e-07 </td>
   <td style="text-align:left;"> 2.85e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000177721_ANXA2R.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV35 </td>
   <td style="text-align:left;"> -5.34 </td>
   <td style="text-align:left;"> 9.12e-08 </td>
   <td style="text-align:left;"> 3.29e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211814_TRAV35.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PSMB8-AS1 </td>
   <td style="text-align:left;"> -5.34 </td>
   <td style="text-align:left;"> 9.22e-08 </td>
   <td style="text-align:left;"> 3.32e-03 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000204261_PSMB8-AS1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> ACTG1 </td>
   <td style="text-align:left;"> -10.05 </td>
   <td style="text-align:left;"> 9.39e-24 </td>
   <td style="text-align:left;"> 3.42e-19 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000184009_ACTG1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TRAV21 </td>
   <td style="text-align:left;"> -4.9 </td>
   <td style="text-align:left;"> 9.59e-07 </td>
   <td style="text-align:left;"> 3.44e-02 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000211801_TRAV21.pdf) </td>
  </tr>
</tbody>
</table>

## 4. Additional genes/proteins






<table class=" lightable-paper lightable-striped" style="font-family: Helvetica; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> celltype </th>
   <th style="text-align:left;"> hgnc_symbol </th>
   <th style="text-align:left;"> z </th>
   <th style="text-align:left;"> pv </th>
   <th style="text-align:left;"> fwer </th>
   <th style="text-align:left;"> .link </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="6"> mTconv </td>
   <td style="text-align:left;"> PDCD1 </td>
   <td style="text-align:left;"> 4.36 </td>
   <td style="text-align:left;"> 1.30e-05 </td>
   <td style="text-align:left;"> 4.64e-01 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000188389_PDCD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LAG3 </td>
   <td style="text-align:left;"> 2.09 </td>
   <td style="text-align:left;"> 3.69e-02 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000089692_LAG3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD279 </td>
   <td style="text-align:left;"> 6.93 </td>
   <td style="text-align:left;"> 4.30e-12 </td>
   <td style="text-align:left;"> 1.56e-07 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD279_CD279.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX </td>
   <td style="text-align:left;"> 0.78 </td>
   <td style="text-align:left;"> 4.38e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198846_TOX.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD366 </td>
   <td style="text-align:left;"> 4.53 </td>
   <td style="text-align:left;"> 5.93e-06 </td>
   <td style="text-align:left;"> 2.12e-01 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD366_CD366.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX2 </td>
   <td style="text-align:left;"> -0.19 </td>
   <td style="text-align:left;"> 8.49e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000124191_TOX2.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="6"> mTreg </td>
   <td style="text-align:left;"> TOX2 </td>
   <td style="text-align:left;"> -1.51 </td>
   <td style="text-align:left;"> 1.32e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000124191_TOX2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PDCD1 </td>
   <td style="text-align:left;"> -0.88 </td>
   <td style="text-align:left;"> 3.77e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000188389_PDCD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> LAG3 </td>
   <td style="text-align:left;"> 11.41 </td>
   <td style="text-align:left;"> 3.85e-30 </td>
   <td style="text-align:left;"> 1.40e-25 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000089692_LAG3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD279 </td>
   <td style="text-align:left;"> 13.41 </td>
   <td style="text-align:left;"> 5.14e-41 </td>
   <td style="text-align:left;"> 1.88e-36 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD279_CD279.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX </td>
   <td style="text-align:left;"> 2.76 </td>
   <td style="text-align:left;"> 5.77e-03 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198846_TOX.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD366 </td>
   <td style="text-align:left;"> 22.74 </td>
   <td style="text-align:left;"> 1.81e-114 </td>
   <td style="text-align:left;"> 6.60e-110 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD366_CD366.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="5"> nTconv </td>
   <td style="text-align:left;"> LAG3 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000089692_LAG3.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX </td>
   <td style="text-align:left;"> 1.23 </td>
   <td style="text-align:left;"> 2.20e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198846_TOX.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD366 </td>
   <td style="text-align:left;"> 2.03 </td>
   <td style="text-align:left;"> 4.22e-02 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD366_CD366.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX2 </td>
   <td style="text-align:left;"> -0.49 </td>
   <td style="text-align:left;"> 6.24e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000124191_TOX2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD279 </td>
   <td style="text-align:left;"> 1.71 </td>
   <td style="text-align:left;"> 8.68e-02 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD279_CD279.pdf) </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align: top !important;" rowspan="5"> nTreg </td>
   <td style="text-align:left;"> CD366 </td>
   <td style="text-align:left;"> 1.61 </td>
   <td style="text-align:left;"> 1.07e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD366_CD366.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> CD279 </td>
   <td style="text-align:left;"> 6.3 </td>
   <td style="text-align:left;"> 2.89e-10 </td>
   <td style="text-align:left;"> 1.05e-05 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_CD279_CD279.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX2 </td>
   <td style="text-align:left;"> 0.95 </td>
   <td style="text-align:left;"> 3.43e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000124191_TOX2.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> PDCD1 </td>
   <td style="text-align:left;"> 0.76 </td>
   <td style="text-align:left;"> 4.46e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000188389_PDCD1.pdf) </td>
  </tr>
  <tr>
   
   <td style="text-align:left;"> TOX </td>
   <td style="text-align:left;"> -0.06 </td>
   <td style="text-align:left;"> 9.55e-01 </td>
   <td style="text-align:left;"> 1.00e+00 </td>
   <td style="text-align:left;"> [PDF](Fig/STEP3//example/Fig_DEG_example_ENSG00000198846_TOX.pdf) </td>
  </tr>
</tbody>
</table>
