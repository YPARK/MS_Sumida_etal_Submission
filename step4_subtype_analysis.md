---
title: "Step 4: subtype analysis"
author: Yongjin Park
theme: jekyll-theme-minimal
date: "2024-02-07"
bibliography: "MS_Ref.bib"
output:
  html_document:
    toc: true
    keep_md: true
    self_contained: true
---




* Marker genes well-known for mTconv cell subtype classification


```r
.markers <-
    c("CD3", "CD3E", "CD14", "LYZ",
      "CCR4","CCR6","CXCR3","CXCR5",
      "ABCA1","GBP4","CDCA7L","ITM2C","NR1D1",
      "CTSH","GATA3","FHIT","CD40LG","C1orf162",
      "MGATA4A","GZMK","IFNGR2",
      "CD183","CD184","CD185","CD196","CD195","CD194",
      "RORC", "TBX21", "HLADR", "CD74", "TCF7", "LEF1",
      "SELL", "CCR7", "CCR8", "IKZF2", "TIGIT", "CD226",
      "BATF", "ANXA2", "BRD9", "HPGD", "LMNA", "TNFRSF4",
      "FOXP3", "FOXP1", "PDCD1", "CD279", "CTLA4", "LAG3",
      "HAVCR2", "CD366", "KLRB1", "FOSL2", "S100A4", "GMAP7",
      "JUN", "IL7R", "MYC", "IL32", "ISG20", "MALAT1",
      "GSDMD", "HDAC1", "GIMAP4", "APOBEC3G", "CD2", "CD28", "CD6",
      "CDKN2A", "CORO1A", "FAS", "FLI1", "GPR25", "MT2A", "KEAP1",
      "IL12RB", "SIRT2", "TNFRSF14", "TRAF3IP3", "IRF2", "PSPH",
      "CD278", "B2M", "RPS26", "MAP1S", "SGK1", "BACH2",
      "HLA-C", "HLA-B", "HLA-E", "HLA-DR", "HLA-DRA", "HLA-DRB1",
      "S1PR4", "KLF2", "SATB1", "TSC22D3", "IL2RA", "CD25") %>%
    unique
```

* Goal: Identify cellular states/subtypes in memory T cells


```r
.hash.hdr <- "result/step1/hash"
.hash.data <- fileset.list(.hash.hdr)
.hash.info <- read.hash(.hash.data)
```


```r
annot.dt <- fread("Tab/step2_celltype.txt.gz")
```



## 1. Memory T conventional


```r
.full.data <- fileset.list("result/step2/matrix_final")
.mkdir("result/step4/")
.data <- fileset.list("result/step4/mtconv")

if.needed(.data, {
    .tags <- unique(annot.dt[celltype == "mTconv"]$tag)
    .data <-
        rcpp_mmutil_copy_selected_columns(.full.data$mtx,
                                          .full.data$row,
                                          .full.data$col,
                                          .tags,
                                          "result/step4/mtconv")
})
```

### A. Clustering cells by batch-balancing k-nearest neighbour graph


```r
.file <- "result/step4/mtconv_bbknn.rds"

if.needed(.file, {
    .batches <- take.batch.info(.data)

    .bbknn <-
        rcpp_mmutil_bbknn_mtx(.data$mtx,
                              r_batches = .batches, # batch label
                              RANK = 30,            # PCs
                              knn = 50,             # 20 nn per batch
                              RECIPROCAL_MATCH = T, # crucial
                              EM_ITER = 20,         # EM steps
                              NUM_THREADS = 16,
                              TAKE_LN = T,
                              USE_SINGULAR_VALUES = F)

    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```


```r
.file <- "Tab/step4_mtconv_leiden.txt.gz"
if.needed(.file, {
    .tags <- readLines(.data$col)
    .leiden <- run.leiden(.bbknn$knn.adj, .tags, res=.3, nrepeat = 100, min.size = 10)
    fwrite(.leiden, .file)
})
.leiden <- fread(.file)
```

[**DOWNLOAD:** mTconv Leiden results](Tab/step4_mtconv_leiden.txt.gz)


```r
.file <- "Tab/step4_tumap_mtconv.txt.gz"
if.needed(.file, {

    set.seed(1)
    .umap <- uwot::tumap(.bbknn$factors.adjusted,
                         learning_rate=.1,
                         n_epochs=3000,
                         n_sgd_threads=16,
                         verbose=T,
                         init="lvrandom",
                         scale=T)

    .tags <- readLines(.data$col)

    colnames(.umap) <- "UMAP" %&% 1:ncol(.umap)

    .umap.dt <-
        data.table(.umap, tag = .tags) %>%
        left_join(.leiden) %>%
        na.omit()

    fwrite(.umap.dt, .file)
})
.umap.dt <- fread(.file)
```

[**DOWNLOAD:** mTconv UMAP results](Tab/step4_tumap_mtconv.txt.gz)


```r
.file <- "Tab/step4_tsne_mtconv.txt.gz"
if.needed(.file, {
    .tsne <- Rtsne::Rtsne(.bbknn$factors.adjusted,
                          check_duplicates = FALSE,
                          verbose = T,
                          num_threads = 16)

    .tags <- readLines(.data$col)

    colnames(.tsne$Y) <- "tSNE" %&% 1:ncol(.tsne$Y)

    .tsne.dt <- data.table(.tsne$Y, tag = .tags) %>%
        left_join(.leiden) %>%
        na.omit()

    fwrite(.tsne.dt, .file)
})
.tsne.dt <- fread(.file)
```

[**DOWNLOAD:** mTconv tSNE results](Tab/step4_tsne_mtconv.txt.gz)


### B. What are the cell-cluster-specific marker genes?


```r
.mkdir("Tab/")
.file <- "Tab/step4_mtconv_gene_stat.txt.gz"
if.needed(.file, {
    x <- bbknn.x(.data, .bbknn)
    marker.stat <- take.marker.stats(x, .leiden)
    fwrite(marker.stat, .file, sep = "\t", col.names = T)
})
marker.stat <- fread(.file, sep = "\t")
```

[**DOWNLOAD:** mTconv marker gene statistics](Tab/step4_mtconv_gene_stat.txt.gz)



### C. Non-linear embedding to confirm the cell clusters of mTconv cells


```r
.cells <-
    left_join(.umap.dt, .tsne.dt) %>%
    left_join(.leiden) %>%
    left_join(.hash.info) %>%
    na.omit()

.lab <-
    .cells[,
           .(UMAP1=median(UMAP1),
             UMAP2=median(UMAP2),
             tSNE1=median(tSNE1),
             tSNE2=median(tSNE2)),
           by = .(component, membership)]

.cols <- .more.colors(nrow(.lab), nc.pal=12)

p1 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(membership))) +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

p2 <-
    .gg.plot(.cells, aes(tSNE1, tSNE2, color=as.factor(membership))) +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

plt <- p1 | p2
print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtconv-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtconv.pdf)


```r
.cols <- .more.colors(10, nc.pal=7, .palette="Set1")

p1 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")

p2 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")

plt <- p1 | p2
print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtconv_sub-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtconv_sub.pdf)

### D. Summary heatmap

**NOTE** The colors are standardized `log1p` expression across genes and cells.


```r
x.melt <- bbknn.x.melt(.data, .bbknn, .markers)
.dt <- x.melt %>% left_join(.cells) %>% na.omit()
.sum.subj <- .dt[, .(x = median(x)), by = .(gene, subject, membership)]
.sum.subj[, x := scale(x), by = .(gene)]
```


```r
.sum <-
    .sum.subj[, .(x = median(x)), by = .(gene, membership)] %>%
    mutate(col = `gene`, row = membership, weight = x) %>%
    col.order(1:10, TRUE) %>%
    as.data.table()

plt <-
    .gg.plot(.sum, aes(row, col, fill=pmin(pmax(weight, -1.5), 1.5)))+
    geom_tile(linewidth=.1, color="black") +
    scale_fill_distiller("", palette = "RdBu", direction = -1) +
    theme(legend.key.width = unit(.2,"lines")) +
    theme(legend.key.height = unit(.5,"lines")) +
    xlab("cell clusters") + ylab("features")
print(plt)
```

![](Fig/STEP4/Fig_mtconv_sum_membership-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_mtconv_sum_membership.pdf)






```r
.dt <- copy(.sum.subj) %>%
    mutate(gene = factor(`gene`, .marker.order)) %>%
    mutate(t = subject %&% "." %&% membership)

plt <-
    .gg.plot(.dt, aes(`t`, `gene`, fill=pmin(pmax(`x`, -1.5), 1.5))) +
    facet_grid(. ~ membership, space="free", scales="free")+
    geom_tile(linewidth=.1, color="black") +
    scale_fill_distiller("", palette = "RdBu", direction = -1) +
    theme(legend.key.width = unit(.2,"lines")) +
    theme(legend.key.height = unit(.5,"lines")) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    xlab("subjects") + ylab("features")

print(plt)
```

![](Fig/STEP4/Fig_mtconv_sum_subj_member-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_mtconv_sum_subj_member.pdf)

#### UMAP for each marker gene


```r
for(g in unique(x.melt$gene)) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(UMAP1, UMAP2, color=pmax(pmin(x, 3), -3))

    plt <-
        .gg.plot(.dt[order(`x`)], .aes) +
        xlab("UMAP1") + ylab("UMAP2") +
        ggrastr::rasterise(geom_point(stroke = 0, size=.7), dpi=300) +
        theme(legend.key.width = unit(.2,"lines")) +
        theme(legend.key.height = unit(.5,"lines")) +
        scale_color_distiller(g, palette = "RdBu", direction = -1) +
        ggtitle(g)

    print(plt)
    .file <- fig.dir %&% "/Fig_mtconv_gene_umap" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5)
}
```

![](Fig/STEP4/Fig_mtconv_gene_umap-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD14.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD183.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-3.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD184.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-4.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD185.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-5.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD194.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-6.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD195.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-7.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD196.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-8.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD226.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-9.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD25.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-10.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD278.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-11.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD279.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-12.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD366.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-13.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-14.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIL32.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-15.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTRAF3IP3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-16.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD6.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-17.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD74.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-18.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFAS.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-19.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapBRD9.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-20.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIKZF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-21.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFOXP3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-22.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapSIRT2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-23.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTBX21.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-24.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFOSL2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-25.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapKEAP1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-26.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTCF7.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-27.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapLAG3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-28.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapLYZ.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-29.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD40LG.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-30.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCORO1A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-31.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCTSH.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-32.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGSDMD.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-33.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGATA3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-34.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapKLRB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-35.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapBACH2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-36.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCCR6.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-37.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGZMK.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-38.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFOXP1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-39.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHDAC1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-40.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-41.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapSGK1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-42.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapMT2A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-43.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapS1PR4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-44.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCCR7.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-45.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapNR1D1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-46.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapKLF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-47.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapMAP1S.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-48.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGIMAP4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-49.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIL2RA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-50.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHAVCR2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-51.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapITM2C.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-52.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapMYC.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-53.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapLEF1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-54.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapC1orf162.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-55.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapRORC.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-56.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapPSPH.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-57.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCDKN2A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-58.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFLI1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-59.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapBATF.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-60.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTSC22D3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-61.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTNFRSF14.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-62.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIFNGR2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-63.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCXCR5.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-64.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapLMNA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-65.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGBP4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-66.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCTLA4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-67.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHPGD.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-68.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCDCA7L.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-69.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapABCA1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-70.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapB2M.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-71.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIRF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-72.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapIL7R.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-73.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapGPR25.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-74.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapISG20.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-75.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapJUN.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-76.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD28.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-77.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCCR8.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-78.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTIGIT.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-79.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapSATB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-80.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapANXA2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-81.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCCR4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-82.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCXCR3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-83.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapTNFRSF4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-84.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapPDCD1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-85.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapSELL.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-86.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapFHIT.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-87.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-DRB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-88.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapS100A4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-89.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapRPS26.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-90.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapCD3E.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-91.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-DRA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-92.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-C.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-93.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-E.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-94.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-B.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-95.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapAPOBEC3G.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-96.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapMALAT1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_umap-97.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_umapHLA-DR.pdf)

#### tSNE for each marker gene


```r
for(g in unique(x.melt$gene)) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(tSNE1, tSNE2, color=pmax(pmin(x, 3), -3))

    plt <-
        .gg.plot(.dt[order(`x`)], .aes) +
        xlab("TSNE1") + ylab("TSNE2") +
        ggrastr::rasterise(geom_point(stroke = 0, size=.7), dpi=300) +
        theme(legend.key.width = unit(.2,"lines")) +
        theme(legend.key.height = unit(.5,"lines")) +
        scale_color_distiller(g, palette = "RdBu", direction = -1) +
        ggtitle(g)

    print(plt)
    .file <- fig.dir %&% "/Fig_mtconv_gene_tsne" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5)
}
```

![](Fig/STEP4/Fig_mtconv_gene_tsne-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD14.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD183.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-3.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD184.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-4.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD185.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-5.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD194.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-6.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD195.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-7.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD196.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-8.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD226.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-9.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD25.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-10.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD278.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-11.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD279.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-12.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD366.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-13.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-14.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIL32.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-15.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTRAF3IP3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-16.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD6.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-17.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD74.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-18.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFAS.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-19.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneBRD9.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-20.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIKZF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-21.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFOXP3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-22.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneSIRT2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-23.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTBX21.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-24.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFOSL2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-25.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneKEAP1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-26.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTCF7.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-27.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneLAG3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-28.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneLYZ.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-29.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD40LG.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-30.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCORO1A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-31.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCTSH.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-32.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGSDMD.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-33.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGATA3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-34.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneKLRB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-35.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneBACH2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-36.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCCR6.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-37.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGZMK.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-38.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFOXP1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-39.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHDAC1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-40.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-41.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneSGK1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-42.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneMT2A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-43.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneS1PR4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-44.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCCR7.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-45.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneNR1D1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-46.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneKLF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-47.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneMAP1S.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-48.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGIMAP4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-49.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIL2RA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-50.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHAVCR2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-51.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneITM2C.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-52.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneMYC.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-53.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneLEF1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-54.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneC1orf162.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-55.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneRORC.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-56.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsnePSPH.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-57.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCDKN2A.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-58.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFLI1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-59.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneBATF.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-60.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTSC22D3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-61.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTNFRSF14.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-62.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIFNGR2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-63.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCXCR5.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-64.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneLMNA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-65.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGBP4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-66.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCTLA4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-67.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHPGD.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-68.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCDCA7L.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-69.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneABCA1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-70.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneB2M.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-71.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIRF2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-72.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneIL7R.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-73.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneGPR25.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-74.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneISG20.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-75.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneJUN.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-76.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD28.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-77.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCCR8.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-78.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTIGIT.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-79.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneSATB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-80.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneANXA2.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-81.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCCR4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-82.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCXCR3.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-83.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneTNFRSF4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-84.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsnePDCD1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-85.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneSELL.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-86.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneFHIT.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-87.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-DRB1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-88.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneS100A4.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-89.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneRPS26.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-90.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneCD3E.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-91.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-DRA.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-92.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-C.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-93.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-E.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-94.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-B.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-95.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneAPOBEC3G.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-96.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneMALAT1.pdf)

![](Fig/STEP4/Fig_mtconv_gene_tsne-97.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtconv_gene_tsneHLA-DR.pdf)

### E. Basic statistics


```r
.stat <-
    .cells[,
           .(N = .N),
           by=.(batch, membership, component, disease)] %>%
    mutate(membership = as.factor(membership)) %>% 
    .sum.stat.batch()

plt <- .plt.sum.stat(.stat) + ggtitle("mTconv")
print(plt)
```

![](Fig/STEP4/Fig_count_mtconv_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_mtconv_tot.pdf)


```r
.stat.tot <-
    .cells[,
           .(N = .N),
           by=.(membership, disease)] %>%
    mutate(membership = as.factor(membership)) %>% 
    .sum.stat.tot() %>%
    mutate(batch = "(N=" %&% num.int(sum(.stat$N)) %&% ")")

plt <- .plt.sum.stat(.stat.tot) + ggtitle("mTconv")
print(plt)
```

![](Fig/STEP4/Fig_count_merged_mtconv_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_merged_mtconv_tot.pdf)


## 2. Memory Treg cells


```r
.full.data <- fileset.list("result/step2/matrix_final")
.mkdir("result/step4/")
.data <- fileset.list("result/step4/mtreg")

if.needed(.data, {
    .tags <- unique(annot.dt[celltype == "mTreg"]$tag)
    .data <-
        rcpp_mmutil_copy_selected_columns(.full.data$mtx,
                                          .full.data$row,
                                          .full.data$col,
                                          .tags,
                                          "result/step4/mtreg")
})
```

### A. Clustering cells by batch-balancing k-nearest neighbour graph


```r
.file <- "result/step4/mtreg_bbknn.rds"

if.needed(.file, {
    .batches <- take.batch.info(.data)

    .bbknn <-
        rcpp_mmutil_bbknn_mtx(.data$mtx,
                              r_batches = .batches, # batch label
                              RANK = 30,            # PCs
                              knn = 50,             # 20 nn per batch
                              RECIPROCAL_MATCH = T, # crucial
                              EM_ITER = 20,         # EM steps
                              NUM_THREADS = 16,
                              TAKE_LN = T,
                              USE_SINGULAR_VALUES = F)

    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```


```r
.file <- "Tab/step4_mtreg_leiden.txt.gz"
if.needed(.file, {
    .tags <- readLines(.data$col)
    .leiden <- run.leiden(.bbknn$knn.adj, .tags, res=.3, nrepeat = 100, min.size = 10)
    fwrite(.leiden, .file)
})
.leiden <- fread(.file)
```

[**DOWNLOAD:** mTreg Leiden results](Tab/step4_mtreg_leiden.txt.gz)


```r
.file <- "Tab/step4_tumap_mtreg.txt.gz"
if.needed(.file, {

    set.seed(1)
    .umap <- uwot::tumap(.bbknn$factors.adjusted,
                         learning_rate=.1,
                         n_epochs=3000,
                         n_sgd_threads=16,
                         verbose=T,
                         init="lvrandom",
                         scale=T)

    .tags <- readLines(.data$col)

    colnames(.umap) <- "UMAP" %&% 1:ncol(.umap)

    .umap.dt <-
        data.table(.umap, tag = .tags) %>%
        left_join(.leiden) %>%
        na.omit()

    fwrite(.umap.dt, .file)
})
.umap.dt <- fread(.file)
```

[**DOWNLOAD:** mtreg UMAP results](Tab/step4_tumap_mtreg.txt.gz)


```r
.file <- "Tab/step4_tsne_mtreg.txt.gz"
if.needed(.file, {

    .tsne <- Rtsne::Rtsne(.bbknn$factors.adjusted,
                          check_duplicates = FALSE,
                          verbose = T,
                          num_threads = 16)

    .tags <- readLines(.data$col)

    colnames(.tsne$Y) <- "tSNE" %&% 1:ncol(.tsne$Y)

    .tsne.dt <- data.table(.tsne$Y, tag = .tags) %>%
        left_join(.leiden) %>%
        na.omit()

    fwrite(.tsne.dt, .file)
})
.tsne.dt <- fread(.file)
```

[**DOWNLOAD:** mtreg tSNE results](Tab/step4_tsne_mtreg.txt.gz)


### B. What are the cell-cluster-specific marker genes?


```r
.mkdir("Tab/")
.file <- "Tab/step4_mtreg_gene_stat.txt.gz"
if.needed(.file, {
    x <- bbknn.x(.data, .bbknn)
    marker.stat <- take.marker.stats(x, .leiden)
    fwrite(marker.stat, .file, sep = "\t", col.names = T)
})
marker.stat <- fread(.file, sep = "\t")
```

[**DOWNLOAD:** mTreg marker gene statistics](Tab/step4_mtreg_gene_stat.txt.gz)

### C. Non-linear embedding to confirm the cell clusters of mtreg cells


```r
.cells <-
    left_join(.umap.dt, .tsne.dt) %>%
    left_join(.leiden) %>%
    left_join(.hash.info) %>%
    na.omit()

.lab <-
    .cells[,
           .(UMAP1=median(UMAP1),
             UMAP2=median(UMAP2),
             tSNE1=median(tSNE1),
             tSNE2=median(tSNE2)),
           by = .(component, membership)]

.cols <- .more.colors(nrow(.lab), nc.pal=12)

p1 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(membership))) +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

p2 <-
    .gg.plot(.cells, aes(tSNE1, tSNE2, color=as.factor(membership))) +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

plt <- p1 | p2
print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtreg-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtreg.pdf)


```r
.cols <- .more.colors(10, nc.pal=7, .palette="Set1")

p1 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")

p2 <-
    .gg.plot(.cells, aes(UMAP1, UMAP2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")

plt <- p1 | p2
print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtreg_sub-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtreg_sub.pdf)

### D. Summary heatmap

**NOTE** The colors are standardized `log1p` expression across genes and cells.


```r
x.melt <- bbknn.x.melt(.data, .bbknn, .markers)
.dt <- x.melt %>% left_join(.cells) %>% na.omit()
.sum.subj <- .dt[, .(x = median(x)), by = .(gene, subject, membership)]
.sum.subj[, x := scale(x), by = .(gene)]
```


```r
.sum <-
    .sum.subj[, .(x = median(x)), by = .(gene, membership)] %>%
    mutate(col = `gene`, row = membership, weight = x) %>%
    col.order(1:10, TRUE) %>%
    as.data.table()

plt <-
    .gg.plot(.sum, aes(row, col, fill=pmin(pmax(weight, -1.5), 1.5)))+
    geom_tile(linewidth=.1, color="black") +
    scale_fill_distiller("", palette = "RdBu", direction = -1) +
    theme(legend.key.width = unit(.2,"lines")) +
    theme(legend.key.height = unit(.5,"lines")) +
    xlab("cell clusters") + ylab("features")
print(plt)
```

![](Fig/STEP4/Fig_mtreg_sum_membership-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_mtreg_sum_membership.pdf)




```r
.dt <- copy(.sum.subj) %>%
    mutate(gene = factor(`gene`, .marker.order)) %>%
    mutate(t = subject %&% "." %&% membership)

plt <-
    .gg.plot(.dt, aes(`t`, `gene`, fill=pmin(pmax(`x`, -1.5), 1.5))) +
    facet_grid(. ~ membership, space="free", scales="free")+
    geom_tile(linewidth=.1, color="black") +
    scale_fill_distiller("", palette = "RdBu", direction = -1) +
    theme(legend.key.width = unit(.2,"lines")) +
    theme(legend.key.height = unit(.5,"lines")) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    xlab("subjects") + ylab("features")

print(plt)
```

![](Fig/STEP4/Fig_mtreg_sum_subj_member-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_mtreg_sum_subj_member.pdf)

#### UMAP for each marker gene


```r
for(g in unique(x.melt$gene)) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(UMAP1, UMAP2, color=pmax(pmin(x, 3), -3))

    plt <-
        .gg.plot(.dt[order(`x`)], .aes) +
        xlab("UMAP1") + ylab("UMAP2") +
        ggrastr::rasterise(geom_point(stroke = 0, size=.7), dpi=300) +
        theme(legend.key.width = unit(.2,"lines")) +
        theme(legend.key.height = unit(.5,"lines")) +
        scale_color_distiller(g, palette = "RdBu", direction = -1) +
        ggtitle(g)

    print(plt)
    .file <- fig.dir %&% "/Fig_mtreg_gene_umap" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5)
}
```

![](Fig/STEP4/Fig_mtreg_gene_umap-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD183.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-3.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD184.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-4.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD185.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-5.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD194.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-6.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD195.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-7.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD196.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-8.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD226.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-9.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-10.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD278.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-11.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD279.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-12.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD366.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-13.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-14.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIL32.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-15.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTRAF3IP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-16.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-17.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD74.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-18.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFAS.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-19.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapBRD9.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-20.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIKZF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-21.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFOXP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-22.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapSIRT2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-23.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTBX21.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-24.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFOSL2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-25.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapKEAP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-26.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTCF7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-27.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapLAG3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-28.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapLYZ.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-29.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD40LG.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-30.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCORO1A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-31.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCTSH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-32.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGSDMD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-33.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGATA3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-34.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapKLRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-35.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapBACH2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-36.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCCR6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-37.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGZMK.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-38.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFOXP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-39.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHDAC1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-40.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-41.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapSGK1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-42.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapMT2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-43.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapS1PR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-44.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCCR7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-45.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapNR1D1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-46.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapKLF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-47.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapMAP1S.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-48.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGIMAP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-49.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIL2RA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-50.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHAVCR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-51.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapITM2C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-52.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapMYC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-53.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapLEF1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-54.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapC1orf162.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-55.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapRORC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-56.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapPSPH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-57.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCDKN2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-58.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFLI1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-59.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapBATF.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-60.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTSC22D3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-61.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTNFRSF14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-62.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIFNGR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-63.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCXCR5.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-64.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapLMNA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-65.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGBP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-66.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCTLA4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-67.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHPGD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-68.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCDCA7L.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-69.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapABCA1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-70.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapB2M.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-71.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIRF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-72.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapIL7R.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-73.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapGPR25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-74.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapISG20.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-75.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapJUN.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-76.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD28.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-77.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCCR8.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-78.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTIGIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-79.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapSATB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-80.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapANXA2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-81.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCCR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-82.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCXCR3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-83.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapTNFRSF4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-84.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapPDCD1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-85.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapSELL.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-86.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapFHIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-87.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-DRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-88.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapS100A4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-89.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapRPS26.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-90.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapCD3E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-91.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-DRA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-92.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-93.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-94.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-B.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-95.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapAPOBEC3G.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-96.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapMALAT1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-97.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umapHLA-DR.pdf)

#### tSNE for each marker gene


```r
for(g in unique(x.melt$gene)) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(tSNE1, tSNE2, color=pmax(pmin(x, 3), -3))

    plt <-
        .gg.plot(.dt[order(`x`)], .aes) +
        xlab("TSNE1") + ylab("TSNE2") +
        ggrastr::rasterise(geom_point(stroke = 0, size=.7), dpi=300) +
        theme(legend.key.width = unit(.2,"lines")) +
        theme(legend.key.height = unit(.5,"lines")) +
        scale_color_distiller(g, palette = "RdBu", direction = -1) +
        ggtitle(g)

    print(plt)
    .file <- fig.dir %&% "/Fig_mtreg_gene_tsne" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5)
}
```

![](Fig/STEP4/Fig_mtreg_gene_tsne-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD183.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-3.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD184.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-4.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD185.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-5.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD194.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-6.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD195.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-7.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD196.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-8.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD226.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-9.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-10.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD278.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-11.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD279.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-12.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD366.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-13.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-14.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIL32.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-15.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTRAF3IP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-16.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-17.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD74.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-18.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFAS.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-19.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneBRD9.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-20.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIKZF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-21.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFOXP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-22.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneSIRT2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-23.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTBX21.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-24.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFOSL2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-25.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneKEAP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-26.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTCF7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-27.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneLAG3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-28.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneLYZ.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-29.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD40LG.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-30.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCORO1A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-31.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCTSH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-32.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGSDMD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-33.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGATA3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-34.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneKLRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-35.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneBACH2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-36.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCCR6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-37.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGZMK.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-38.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFOXP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-39.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHDAC1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-40.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-41.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneSGK1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-42.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneMT2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-43.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneS1PR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-44.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCCR7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-45.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneNR1D1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-46.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneKLF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-47.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneMAP1S.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-48.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGIMAP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-49.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIL2RA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-50.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHAVCR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-51.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneITM2C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-52.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneMYC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-53.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneLEF1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-54.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneC1orf162.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-55.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneRORC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-56.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsnePSPH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-57.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCDKN2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-58.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFLI1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-59.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneBATF.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-60.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTSC22D3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-61.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTNFRSF14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-62.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIFNGR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-63.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCXCR5.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-64.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneLMNA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-65.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGBP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-66.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCTLA4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-67.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHPGD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-68.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCDCA7L.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-69.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneABCA1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-70.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneB2M.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-71.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIRF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-72.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneIL7R.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-73.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneGPR25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-74.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneISG20.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-75.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneJUN.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-76.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD28.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-77.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCCR8.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-78.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTIGIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-79.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneSATB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-80.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneANXA2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-81.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCCR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-82.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCXCR3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-83.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneTNFRSF4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-84.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsnePDCD1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-85.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneSELL.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-86.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneFHIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-87.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-DRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-88.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneS100A4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-89.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneRPS26.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-90.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneCD3E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-91.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-DRA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-92.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-93.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-94.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-B.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-95.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneAPOBEC3G.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-96.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneMALAT1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_tsne-97.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_tsneHLA-DR.pdf)

### E. Basic statistics


```r
.stat <-
    .cells[,
           .(N = .N),
           by=.(batch, membership, component, disease)] %>%
    mutate(membership = as.factor(membership)) %>% 
    .sum.stat.batch()

plt <- .plt.sum.stat(.stat) + ggtitle("mtreg")
print(plt)
```

![](Fig/STEP4/Fig_count_mtreg_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_mtreg_tot.pdf)


```r
.stat.tot <-
    .cells[,
           .(N = .N),
           by=.(membership, disease)] %>%
    mutate(membership = as.factor(membership)) %>% 
    .sum.stat.tot() %>%
    mutate(batch = "(N=" %&% num.int(sum(.stat$N)) %&% ")")

plt <- .plt.sum.stat(.stat.tot) + ggtitle("mtreg")
print(plt)
```

![](Fig/STEP4/Fig_count_merged_mtreg_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_merged_mtreg_tot.pdf)

