---
title: "Step 4: subtype analysis"
author: Yongjin Park
theme: jekyll-theme-minimal
date: "2024-02-08"
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

### Perform outlier Q/C to detect and remove batch-specific cells


```r
.file <- "result/step4/mtconv_svd.rds"
if.needed(.file, {
    .svd <- rcpp_mmutil_svd(.data$mtx, RANK=30, TAKE_LN=T, NUM_THREADS = 16, EM_ITER = 20)
    saveRDS(.svd, .file)
})
.svd <- readRDS(.file)
```

#### Can we identify troublesome batch-specific principal components?

##### Show PCs computed on surface protein data


```r
V <- sweep(.svd$V, 2, .svd$D, `*`)
rownames(V) <- readLines(.data$col)
plots <- lapply(1:9, .plot.vd.fun, VD = .pca.df(V))
plt <- wrap_plots(plots, ncol = 3)
print(plt)
```

![](Fig/STEP4/Fig_mtconv_PC-1.png)<!-- -->

* We can adjust systemic batch-specific effects by batch-balancing k-nearest neighbour method. This shouldn't be much trouble.


```r
.file <- "result/step4/mtconv_bbknn.rds"
if.needed(.file, {

    .batches <- take.batch.info(.data)

    .bbknn <-
        rcpp_mmutil_bbknn(r_svd_u = .svd$U,
                          r_svd_v = .svd$V,
                          r_svd_d = .svd$D,
                          r_batches = .batches, # batch label
                          knn = 50,             # 20 nn per batch
                          RECIPROCAL_MATCH = T, # crucial
                          NUM_THREADS = 16,
                          USE_SINGULAR_VALUES = T)

    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```


```r
VD <- .bbknn$factors.adjusted
rownames(VD) <- readLines(.data$col)
plots <- lapply(1:9, .plot.vd.fun, VD = .pca.df(VD))
plt <- wrap_plots(plots, ncol = 3)
print(plt)
```

![](Fig/STEP4/Fig_mtconv_PC_adj-1.png)<!-- -->

### A. Clustering cells by batch-balancing k-nearest neighbour graph


```r
.file <- "Tab/step4_mtconv_leiden_raw.txt.gz"
if.needed(.file, {
    .tags <- readLines(.data$col)
    .leiden <- run.leiden(.bbknn$knn.adj, .tags, res=.7, nrepeat = 100, min.size = 100)
    fwrite(.leiden, .file)
})
.leiden.raw <- fread(.file)
```

#### Filter out small cluster cells if we have enriched specific batches


```r
.file <- "Tab/step4_mtconv_leiden.txt.gz"
if.needed(.file, {
    .leiden <- qc.leiden(.leiden.raw, .cutoff = .75)
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
                         init="spectral",
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
                          num_threads = 16,
                          perplexity = 100)

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

##### Confirm batch/individual-specific effects


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

#### UMAP for each marker gene (normalized expression)


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

    .file <- fig.dir %&% "/Fig_mtconv_gene_umap_" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5, cat.link = F)
    cat("[" %&% g %&% "](" %&% .file %&% ") ")
}
```

[CD14](Fig/STEP4//Fig_mtconv_gene_umap_CD14.pdf) [CD183](Fig/STEP4//Fig_mtconv_gene_umap_CD183.pdf) [CD184](Fig/STEP4//Fig_mtconv_gene_umap_CD184.pdf) [CD185](Fig/STEP4//Fig_mtconv_gene_umap_CD185.pdf) [CD194](Fig/STEP4//Fig_mtconv_gene_umap_CD194.pdf) [CD195](Fig/STEP4//Fig_mtconv_gene_umap_CD195.pdf) [CD196](Fig/STEP4//Fig_mtconv_gene_umap_CD196.pdf) [CD226](Fig/STEP4//Fig_mtconv_gene_umap_CD226.pdf) [CD25](Fig/STEP4//Fig_mtconv_gene_umap_CD25.pdf) [CD278](Fig/STEP4//Fig_mtconv_gene_umap_CD278.pdf) [CD279](Fig/STEP4//Fig_mtconv_gene_umap_CD279.pdf) [CD366](Fig/STEP4//Fig_mtconv_gene_umap_CD366.pdf) [CD3](Fig/STEP4//Fig_mtconv_gene_umap_CD3.pdf) [IL32](Fig/STEP4//Fig_mtconv_gene_umap_IL32.pdf) [TRAF3IP3](Fig/STEP4//Fig_mtconv_gene_umap_TRAF3IP3.pdf) [CD6](Fig/STEP4//Fig_mtconv_gene_umap_CD6.pdf) [CD74](Fig/STEP4//Fig_mtconv_gene_umap_CD74.pdf) [FAS](Fig/STEP4//Fig_mtconv_gene_umap_FAS.pdf) [BRD9](Fig/STEP4//Fig_mtconv_gene_umap_BRD9.pdf) [IKZF2](Fig/STEP4//Fig_mtconv_gene_umap_IKZF2.pdf) [FOXP3](Fig/STEP4//Fig_mtconv_gene_umap_FOXP3.pdf) [SIRT2](Fig/STEP4//Fig_mtconv_gene_umap_SIRT2.pdf) [TBX21](Fig/STEP4//Fig_mtconv_gene_umap_TBX21.pdf) [FOSL2](Fig/STEP4//Fig_mtconv_gene_umap_FOSL2.pdf) [KEAP1](Fig/STEP4//Fig_mtconv_gene_umap_KEAP1.pdf) [TCF7](Fig/STEP4//Fig_mtconv_gene_umap_TCF7.pdf) [LAG3](Fig/STEP4//Fig_mtconv_gene_umap_LAG3.pdf) [LYZ](Fig/STEP4//Fig_mtconv_gene_umap_LYZ.pdf) [CD40LG](Fig/STEP4//Fig_mtconv_gene_umap_CD40LG.pdf) [CORO1A](Fig/STEP4//Fig_mtconv_gene_umap_CORO1A.pdf) [CTSH](Fig/STEP4//Fig_mtconv_gene_umap_CTSH.pdf) [GSDMD](Fig/STEP4//Fig_mtconv_gene_umap_GSDMD.pdf) [GATA3](Fig/STEP4//Fig_mtconv_gene_umap_GATA3.pdf) [KLRB1](Fig/STEP4//Fig_mtconv_gene_umap_KLRB1.pdf) [BACH2](Fig/STEP4//Fig_mtconv_gene_umap_BACH2.pdf) [CCR6](Fig/STEP4//Fig_mtconv_gene_umap_CCR6.pdf) [GZMK](Fig/STEP4//Fig_mtconv_gene_umap_GZMK.pdf) [FOXP1](Fig/STEP4//Fig_mtconv_gene_umap_FOXP1.pdf) [HDAC1](Fig/STEP4//Fig_mtconv_gene_umap_HDAC1.pdf) [CD2](Fig/STEP4//Fig_mtconv_gene_umap_CD2.pdf) [SGK1](Fig/STEP4//Fig_mtconv_gene_umap_SGK1.pdf) [MT2A](Fig/STEP4//Fig_mtconv_gene_umap_MT2A.pdf) [S1PR4](Fig/STEP4//Fig_mtconv_gene_umap_S1PR4.pdf) [CCR7](Fig/STEP4//Fig_mtconv_gene_umap_CCR7.pdf) [NR1D1](Fig/STEP4//Fig_mtconv_gene_umap_NR1D1.pdf) [KLF2](Fig/STEP4//Fig_mtconv_gene_umap_KLF2.pdf) [MAP1S](Fig/STEP4//Fig_mtconv_gene_umap_MAP1S.pdf) [GIMAP4](Fig/STEP4//Fig_mtconv_gene_umap_GIMAP4.pdf) [IL2RA](Fig/STEP4//Fig_mtconv_gene_umap_IL2RA.pdf) [HAVCR2](Fig/STEP4//Fig_mtconv_gene_umap_HAVCR2.pdf) [ITM2C](Fig/STEP4//Fig_mtconv_gene_umap_ITM2C.pdf) [MYC](Fig/STEP4//Fig_mtconv_gene_umap_MYC.pdf) [LEF1](Fig/STEP4//Fig_mtconv_gene_umap_LEF1.pdf) [C1orf162](Fig/STEP4//Fig_mtconv_gene_umap_C1orf162.pdf) [RORC](Fig/STEP4//Fig_mtconv_gene_umap_RORC.pdf) [PSPH](Fig/STEP4//Fig_mtconv_gene_umap_PSPH.pdf) [CDKN2A](Fig/STEP4//Fig_mtconv_gene_umap_CDKN2A.pdf) [FLI1](Fig/STEP4//Fig_mtconv_gene_umap_FLI1.pdf) [BATF](Fig/STEP4//Fig_mtconv_gene_umap_BATF.pdf) [TSC22D3](Fig/STEP4//Fig_mtconv_gene_umap_TSC22D3.pdf) [TNFRSF14](Fig/STEP4//Fig_mtconv_gene_umap_TNFRSF14.pdf) [IFNGR2](Fig/STEP4//Fig_mtconv_gene_umap_IFNGR2.pdf) [CXCR5](Fig/STEP4//Fig_mtconv_gene_umap_CXCR5.pdf) [LMNA](Fig/STEP4//Fig_mtconv_gene_umap_LMNA.pdf) [GBP4](Fig/STEP4//Fig_mtconv_gene_umap_GBP4.pdf) [CTLA4](Fig/STEP4//Fig_mtconv_gene_umap_CTLA4.pdf) [HPGD](Fig/STEP4//Fig_mtconv_gene_umap_HPGD.pdf) [CDCA7L](Fig/STEP4//Fig_mtconv_gene_umap_CDCA7L.pdf) [ABCA1](Fig/STEP4//Fig_mtconv_gene_umap_ABCA1.pdf) [B2M](Fig/STEP4//Fig_mtconv_gene_umap_B2M.pdf) [IRF2](Fig/STEP4//Fig_mtconv_gene_umap_IRF2.pdf) [IL7R](Fig/STEP4//Fig_mtconv_gene_umap_IL7R.pdf) [GPR25](Fig/STEP4//Fig_mtconv_gene_umap_GPR25.pdf) [ISG20](Fig/STEP4//Fig_mtconv_gene_umap_ISG20.pdf) [JUN](Fig/STEP4//Fig_mtconv_gene_umap_JUN.pdf) [CD28](Fig/STEP4//Fig_mtconv_gene_umap_CD28.pdf) [CCR8](Fig/STEP4//Fig_mtconv_gene_umap_CCR8.pdf) [TIGIT](Fig/STEP4//Fig_mtconv_gene_umap_TIGIT.pdf) [SATB1](Fig/STEP4//Fig_mtconv_gene_umap_SATB1.pdf) [ANXA2](Fig/STEP4//Fig_mtconv_gene_umap_ANXA2.pdf) [CCR4](Fig/STEP4//Fig_mtconv_gene_umap_CCR4.pdf) [CXCR3](Fig/STEP4//Fig_mtconv_gene_umap_CXCR3.pdf) [TNFRSF4](Fig/STEP4//Fig_mtconv_gene_umap_TNFRSF4.pdf) [PDCD1](Fig/STEP4//Fig_mtconv_gene_umap_PDCD1.pdf) [SELL](Fig/STEP4//Fig_mtconv_gene_umap_SELL.pdf) [FHIT](Fig/STEP4//Fig_mtconv_gene_umap_FHIT.pdf) [HLA-DRB1](Fig/STEP4//Fig_mtconv_gene_umap_HLA-DRB1.pdf) [S100A4](Fig/STEP4//Fig_mtconv_gene_umap_S100A4.pdf) [RPS26](Fig/STEP4//Fig_mtconv_gene_umap_RPS26.pdf) [CD3E](Fig/STEP4//Fig_mtconv_gene_umap_CD3E.pdf) [HLA-DRA](Fig/STEP4//Fig_mtconv_gene_umap_HLA-DRA.pdf) [HLA-C](Fig/STEP4//Fig_mtconv_gene_umap_HLA-C.pdf) [HLA-E](Fig/STEP4//Fig_mtconv_gene_umap_HLA-E.pdf) [HLA-B](Fig/STEP4//Fig_mtconv_gene_umap_HLA-B.pdf) [APOBEC3G](Fig/STEP4//Fig_mtconv_gene_umap_APOBEC3G.pdf) [MALAT1](Fig/STEP4//Fig_mtconv_gene_umap_MALAT1.pdf) [HLA-DR](Fig/STEP4//Fig_mtconv_gene_umap_HLA-DR.pdf) 

#### tSNE for each marker gene (normalized expression)


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

    .file <- fig.dir %&% "/Fig_mtconv_gene_tsne_" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5, cat.link = F)
    cat("[" %&% g %&% "](" %&% .file %&% ") ")
}
```

[CD14](Fig/STEP4//Fig_mtconv_gene_tsne_CD14.pdf) [CD183](Fig/STEP4//Fig_mtconv_gene_tsne_CD183.pdf) [CD184](Fig/STEP4//Fig_mtconv_gene_tsne_CD184.pdf) [CD185](Fig/STEP4//Fig_mtconv_gene_tsne_CD185.pdf) [CD194](Fig/STEP4//Fig_mtconv_gene_tsne_CD194.pdf) [CD195](Fig/STEP4//Fig_mtconv_gene_tsne_CD195.pdf) [CD196](Fig/STEP4//Fig_mtconv_gene_tsne_CD196.pdf) [CD226](Fig/STEP4//Fig_mtconv_gene_tsne_CD226.pdf) [CD25](Fig/STEP4//Fig_mtconv_gene_tsne_CD25.pdf) [CD278](Fig/STEP4//Fig_mtconv_gene_tsne_CD278.pdf) [CD279](Fig/STEP4//Fig_mtconv_gene_tsne_CD279.pdf) [CD366](Fig/STEP4//Fig_mtconv_gene_tsne_CD366.pdf) [CD3](Fig/STEP4//Fig_mtconv_gene_tsne_CD3.pdf) [IL32](Fig/STEP4//Fig_mtconv_gene_tsne_IL32.pdf) [TRAF3IP3](Fig/STEP4//Fig_mtconv_gene_tsne_TRAF3IP3.pdf) [CD6](Fig/STEP4//Fig_mtconv_gene_tsne_CD6.pdf) [CD74](Fig/STEP4//Fig_mtconv_gene_tsne_CD74.pdf) [FAS](Fig/STEP4//Fig_mtconv_gene_tsne_FAS.pdf) [BRD9](Fig/STEP4//Fig_mtconv_gene_tsne_BRD9.pdf) [IKZF2](Fig/STEP4//Fig_mtconv_gene_tsne_IKZF2.pdf) [FOXP3](Fig/STEP4//Fig_mtconv_gene_tsne_FOXP3.pdf) [SIRT2](Fig/STEP4//Fig_mtconv_gene_tsne_SIRT2.pdf) [TBX21](Fig/STEP4//Fig_mtconv_gene_tsne_TBX21.pdf) [FOSL2](Fig/STEP4//Fig_mtconv_gene_tsne_FOSL2.pdf) [KEAP1](Fig/STEP4//Fig_mtconv_gene_tsne_KEAP1.pdf) [TCF7](Fig/STEP4//Fig_mtconv_gene_tsne_TCF7.pdf) [LAG3](Fig/STEP4//Fig_mtconv_gene_tsne_LAG3.pdf) [LYZ](Fig/STEP4//Fig_mtconv_gene_tsne_LYZ.pdf) [CD40LG](Fig/STEP4//Fig_mtconv_gene_tsne_CD40LG.pdf) [CORO1A](Fig/STEP4//Fig_mtconv_gene_tsne_CORO1A.pdf) [CTSH](Fig/STEP4//Fig_mtconv_gene_tsne_CTSH.pdf) [GSDMD](Fig/STEP4//Fig_mtconv_gene_tsne_GSDMD.pdf) [GATA3](Fig/STEP4//Fig_mtconv_gene_tsne_GATA3.pdf) [KLRB1](Fig/STEP4//Fig_mtconv_gene_tsne_KLRB1.pdf) [BACH2](Fig/STEP4//Fig_mtconv_gene_tsne_BACH2.pdf) [CCR6](Fig/STEP4//Fig_mtconv_gene_tsne_CCR6.pdf) [GZMK](Fig/STEP4//Fig_mtconv_gene_tsne_GZMK.pdf) [FOXP1](Fig/STEP4//Fig_mtconv_gene_tsne_FOXP1.pdf) [HDAC1](Fig/STEP4//Fig_mtconv_gene_tsne_HDAC1.pdf) [CD2](Fig/STEP4//Fig_mtconv_gene_tsne_CD2.pdf) [SGK1](Fig/STEP4//Fig_mtconv_gene_tsne_SGK1.pdf) [MT2A](Fig/STEP4//Fig_mtconv_gene_tsne_MT2A.pdf) [S1PR4](Fig/STEP4//Fig_mtconv_gene_tsne_S1PR4.pdf) [CCR7](Fig/STEP4//Fig_mtconv_gene_tsne_CCR7.pdf) [NR1D1](Fig/STEP4//Fig_mtconv_gene_tsne_NR1D1.pdf) [KLF2](Fig/STEP4//Fig_mtconv_gene_tsne_KLF2.pdf) [MAP1S](Fig/STEP4//Fig_mtconv_gene_tsne_MAP1S.pdf) [GIMAP4](Fig/STEP4//Fig_mtconv_gene_tsne_GIMAP4.pdf) [IL2RA](Fig/STEP4//Fig_mtconv_gene_tsne_IL2RA.pdf) [HAVCR2](Fig/STEP4//Fig_mtconv_gene_tsne_HAVCR2.pdf) [ITM2C](Fig/STEP4//Fig_mtconv_gene_tsne_ITM2C.pdf) [MYC](Fig/STEP4//Fig_mtconv_gene_tsne_MYC.pdf) [LEF1](Fig/STEP4//Fig_mtconv_gene_tsne_LEF1.pdf) [C1orf162](Fig/STEP4//Fig_mtconv_gene_tsne_C1orf162.pdf) [RORC](Fig/STEP4//Fig_mtconv_gene_tsne_RORC.pdf) [PSPH](Fig/STEP4//Fig_mtconv_gene_tsne_PSPH.pdf) [CDKN2A](Fig/STEP4//Fig_mtconv_gene_tsne_CDKN2A.pdf) [FLI1](Fig/STEP4//Fig_mtconv_gene_tsne_FLI1.pdf) [BATF](Fig/STEP4//Fig_mtconv_gene_tsne_BATF.pdf) [TSC22D3](Fig/STEP4//Fig_mtconv_gene_tsne_TSC22D3.pdf) [TNFRSF14](Fig/STEP4//Fig_mtconv_gene_tsne_TNFRSF14.pdf) [IFNGR2](Fig/STEP4//Fig_mtconv_gene_tsne_IFNGR2.pdf) [CXCR5](Fig/STEP4//Fig_mtconv_gene_tsne_CXCR5.pdf) [LMNA](Fig/STEP4//Fig_mtconv_gene_tsne_LMNA.pdf) [GBP4](Fig/STEP4//Fig_mtconv_gene_tsne_GBP4.pdf) [CTLA4](Fig/STEP4//Fig_mtconv_gene_tsne_CTLA4.pdf) [HPGD](Fig/STEP4//Fig_mtconv_gene_tsne_HPGD.pdf) [CDCA7L](Fig/STEP4//Fig_mtconv_gene_tsne_CDCA7L.pdf) [ABCA1](Fig/STEP4//Fig_mtconv_gene_tsne_ABCA1.pdf) [B2M](Fig/STEP4//Fig_mtconv_gene_tsne_B2M.pdf) [IRF2](Fig/STEP4//Fig_mtconv_gene_tsne_IRF2.pdf) [IL7R](Fig/STEP4//Fig_mtconv_gene_tsne_IL7R.pdf) [GPR25](Fig/STEP4//Fig_mtconv_gene_tsne_GPR25.pdf) [ISG20](Fig/STEP4//Fig_mtconv_gene_tsne_ISG20.pdf) [JUN](Fig/STEP4//Fig_mtconv_gene_tsne_JUN.pdf) [CD28](Fig/STEP4//Fig_mtconv_gene_tsne_CD28.pdf) [CCR8](Fig/STEP4//Fig_mtconv_gene_tsne_CCR8.pdf) [TIGIT](Fig/STEP4//Fig_mtconv_gene_tsne_TIGIT.pdf) [SATB1](Fig/STEP4//Fig_mtconv_gene_tsne_SATB1.pdf) [ANXA2](Fig/STEP4//Fig_mtconv_gene_tsne_ANXA2.pdf) [CCR4](Fig/STEP4//Fig_mtconv_gene_tsne_CCR4.pdf) [CXCR3](Fig/STEP4//Fig_mtconv_gene_tsne_CXCR3.pdf) [TNFRSF4](Fig/STEP4//Fig_mtconv_gene_tsne_TNFRSF4.pdf) [PDCD1](Fig/STEP4//Fig_mtconv_gene_tsne_PDCD1.pdf) [SELL](Fig/STEP4//Fig_mtconv_gene_tsne_SELL.pdf) [FHIT](Fig/STEP4//Fig_mtconv_gene_tsne_FHIT.pdf) [HLA-DRB1](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-DRB1.pdf) [S100A4](Fig/STEP4//Fig_mtconv_gene_tsne_S100A4.pdf) [RPS26](Fig/STEP4//Fig_mtconv_gene_tsne_RPS26.pdf) [CD3E](Fig/STEP4//Fig_mtconv_gene_tsne_CD3E.pdf) [HLA-DRA](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-DRA.pdf) [HLA-C](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-C.pdf) [HLA-E](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-E.pdf) [HLA-B](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-B.pdf) [APOBEC3G](Fig/STEP4//Fig_mtconv_gene_tsne_APOBEC3G.pdf) [MALAT1](Fig/STEP4//Fig_mtconv_gene_tsne_MALAT1.pdf) [HLA-DR](Fig/STEP4//Fig_mtconv_gene_tsne_HLA-DR.pdf) 

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

### Perform outlier Q/C to detect and remove batch-specific cells


```r
.file <- "result/step4/mtreg_svd.rds"
if.needed(.file, {
    .svd <- rcpp_mmutil_svd(.data$mtx, RANK=30, TAKE_LN=T, NUM_THREADS = 16, EM_ITER = 20)
    saveRDS(.svd, .file)
})
.svd <- readRDS(.file)
```

#### Can we identify troublesome batch-specific principal components?

##### Show PCs computed on surface protein data


```r
V <- sweep(.svd$V, 2, .svd$D, `*`)
rownames(V) <- readLines(.data$col)
plots <- lapply(1:9, .plot.vd.fun, VD = .pca.df(V))
plt <- wrap_plots(plots, ncol = 3)
print(plt)
```

![](Fig/STEP4/Fig_mtreg_PC-1.png)<!-- -->

* We can adjust systemic batch-specific effects by batch-balancing k-nearest neighbour method. This shouldn't be much trouble.


```r
.file <- "result/step4/mtreg_bbknn.rds"
if.needed(.file, {

    .batches <- take.batch.info(.data)

    .bbknn <-
        rcpp_mmutil_bbknn(r_svd_u = .svd$U,
                          r_svd_v = .svd$V,
                          r_svd_d = .svd$D,
                          r_batches = .batches, # batch label
                          knn = 50,             # 20 nn per batch
                          RECIPROCAL_MATCH = T, # crucial
                          NUM_THREADS = 16,
                          USE_SINGULAR_VALUES = T)

    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```


```r
VD <- .bbknn$factors.adjusted
rownames(VD) <- readLines(.data$col)
plots <- lapply(1:9, .plot.vd.fun, VD = .pca.df(VD))
plt <- wrap_plots(plots, ncol = 3)
print(plt)
```

![](Fig/STEP4/Fig_mtreg_PC_adj-1.png)<!-- -->

### A. Clustering cells by batch-balancing k-nearest neighbour graph



```r
.file <- "Tab/step4_mtreg_leiden_raw.txt.gz"
if.needed(.file, {
    .tags <- readLines(.data$col)
    .leiden <- run.leiden(.bbknn$knn.adj, .tags, res=.7, nrepeat = 100, min.size = 100)
    fwrite(.leiden, .file)
})
.leiden.raw <- fread(.file)
```

#### Filter out small cluster cells if we have enriched specific batches


```r
.file <- "Tab/step4_mtreg_leiden.txt.gz"
if.needed(.file, {
    .leiden <- qc.leiden(.leiden.raw, .cutoff = .75)
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
                         init="spectral",
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

##### Confirm batch/individual-specific effects


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
    col.order(1:15, TRUE) %>%
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

    .file <- fig.dir %&% "/Fig_mtreg_gene_umap_" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5, cat.link = F)
    cat("[" %&% g %&% "](" %&% .file %&% ") ")
}
```

[CD14](Fig/STEP4//Fig_mtreg_gene_umap_CD14.pdf) [CD183](Fig/STEP4//Fig_mtreg_gene_umap_CD183.pdf) [CD184](Fig/STEP4//Fig_mtreg_gene_umap_CD184.pdf) [CD185](Fig/STEP4//Fig_mtreg_gene_umap_CD185.pdf) [CD194](Fig/STEP4//Fig_mtreg_gene_umap_CD194.pdf) [CD195](Fig/STEP4//Fig_mtreg_gene_umap_CD195.pdf) [CD196](Fig/STEP4//Fig_mtreg_gene_umap_CD196.pdf) [CD226](Fig/STEP4//Fig_mtreg_gene_umap_CD226.pdf) [CD25](Fig/STEP4//Fig_mtreg_gene_umap_CD25.pdf) [CD278](Fig/STEP4//Fig_mtreg_gene_umap_CD278.pdf) [CD279](Fig/STEP4//Fig_mtreg_gene_umap_CD279.pdf) [CD366](Fig/STEP4//Fig_mtreg_gene_umap_CD366.pdf) [CD3](Fig/STEP4//Fig_mtreg_gene_umap_CD3.pdf) [IL32](Fig/STEP4//Fig_mtreg_gene_umap_IL32.pdf) [TRAF3IP3](Fig/STEP4//Fig_mtreg_gene_umap_TRAF3IP3.pdf) [CD6](Fig/STEP4//Fig_mtreg_gene_umap_CD6.pdf) [CD74](Fig/STEP4//Fig_mtreg_gene_umap_CD74.pdf) [FAS](Fig/STEP4//Fig_mtreg_gene_umap_FAS.pdf) [BRD9](Fig/STEP4//Fig_mtreg_gene_umap_BRD9.pdf) [IKZF2](Fig/STEP4//Fig_mtreg_gene_umap_IKZF2.pdf) [FOXP3](Fig/STEP4//Fig_mtreg_gene_umap_FOXP3.pdf) [SIRT2](Fig/STEP4//Fig_mtreg_gene_umap_SIRT2.pdf) [TBX21](Fig/STEP4//Fig_mtreg_gene_umap_TBX21.pdf) [FOSL2](Fig/STEP4//Fig_mtreg_gene_umap_FOSL2.pdf) [KEAP1](Fig/STEP4//Fig_mtreg_gene_umap_KEAP1.pdf) [TCF7](Fig/STEP4//Fig_mtreg_gene_umap_TCF7.pdf) [LAG3](Fig/STEP4//Fig_mtreg_gene_umap_LAG3.pdf) [LYZ](Fig/STEP4//Fig_mtreg_gene_umap_LYZ.pdf) [CD40LG](Fig/STEP4//Fig_mtreg_gene_umap_CD40LG.pdf) [CORO1A](Fig/STEP4//Fig_mtreg_gene_umap_CORO1A.pdf) [CTSH](Fig/STEP4//Fig_mtreg_gene_umap_CTSH.pdf) [GSDMD](Fig/STEP4//Fig_mtreg_gene_umap_GSDMD.pdf) [GATA3](Fig/STEP4//Fig_mtreg_gene_umap_GATA3.pdf) [KLRB1](Fig/STEP4//Fig_mtreg_gene_umap_KLRB1.pdf) [BACH2](Fig/STEP4//Fig_mtreg_gene_umap_BACH2.pdf) [CCR6](Fig/STEP4//Fig_mtreg_gene_umap_CCR6.pdf) [GZMK](Fig/STEP4//Fig_mtreg_gene_umap_GZMK.pdf) [FOXP1](Fig/STEP4//Fig_mtreg_gene_umap_FOXP1.pdf) [HDAC1](Fig/STEP4//Fig_mtreg_gene_umap_HDAC1.pdf) [CD2](Fig/STEP4//Fig_mtreg_gene_umap_CD2.pdf) [SGK1](Fig/STEP4//Fig_mtreg_gene_umap_SGK1.pdf) [MT2A](Fig/STEP4//Fig_mtreg_gene_umap_MT2A.pdf) [S1PR4](Fig/STEP4//Fig_mtreg_gene_umap_S1PR4.pdf) [CCR7](Fig/STEP4//Fig_mtreg_gene_umap_CCR7.pdf) [NR1D1](Fig/STEP4//Fig_mtreg_gene_umap_NR1D1.pdf) [KLF2](Fig/STEP4//Fig_mtreg_gene_umap_KLF2.pdf) [MAP1S](Fig/STEP4//Fig_mtreg_gene_umap_MAP1S.pdf) [GIMAP4](Fig/STEP4//Fig_mtreg_gene_umap_GIMAP4.pdf) [IL2RA](Fig/STEP4//Fig_mtreg_gene_umap_IL2RA.pdf) [HAVCR2](Fig/STEP4//Fig_mtreg_gene_umap_HAVCR2.pdf) [ITM2C](Fig/STEP4//Fig_mtreg_gene_umap_ITM2C.pdf) [MYC](Fig/STEP4//Fig_mtreg_gene_umap_MYC.pdf) [LEF1](Fig/STEP4//Fig_mtreg_gene_umap_LEF1.pdf) [C1orf162](Fig/STEP4//Fig_mtreg_gene_umap_C1orf162.pdf) [RORC](Fig/STEP4//Fig_mtreg_gene_umap_RORC.pdf) [PSPH](Fig/STEP4//Fig_mtreg_gene_umap_PSPH.pdf) [CDKN2A](Fig/STEP4//Fig_mtreg_gene_umap_CDKN2A.pdf) [FLI1](Fig/STEP4//Fig_mtreg_gene_umap_FLI1.pdf) [BATF](Fig/STEP4//Fig_mtreg_gene_umap_BATF.pdf) [TSC22D3](Fig/STEP4//Fig_mtreg_gene_umap_TSC22D3.pdf) [TNFRSF14](Fig/STEP4//Fig_mtreg_gene_umap_TNFRSF14.pdf) [IFNGR2](Fig/STEP4//Fig_mtreg_gene_umap_IFNGR2.pdf) [CXCR5](Fig/STEP4//Fig_mtreg_gene_umap_CXCR5.pdf) [LMNA](Fig/STEP4//Fig_mtreg_gene_umap_LMNA.pdf) [GBP4](Fig/STEP4//Fig_mtreg_gene_umap_GBP4.pdf) [CTLA4](Fig/STEP4//Fig_mtreg_gene_umap_CTLA4.pdf) [HPGD](Fig/STEP4//Fig_mtreg_gene_umap_HPGD.pdf) [CDCA7L](Fig/STEP4//Fig_mtreg_gene_umap_CDCA7L.pdf) [ABCA1](Fig/STEP4//Fig_mtreg_gene_umap_ABCA1.pdf) [B2M](Fig/STEP4//Fig_mtreg_gene_umap_B2M.pdf) [IRF2](Fig/STEP4//Fig_mtreg_gene_umap_IRF2.pdf) [IL7R](Fig/STEP4//Fig_mtreg_gene_umap_IL7R.pdf) [GPR25](Fig/STEP4//Fig_mtreg_gene_umap_GPR25.pdf) [ISG20](Fig/STEP4//Fig_mtreg_gene_umap_ISG20.pdf) [JUN](Fig/STEP4//Fig_mtreg_gene_umap_JUN.pdf) [CD28](Fig/STEP4//Fig_mtreg_gene_umap_CD28.pdf) [CCR8](Fig/STEP4//Fig_mtreg_gene_umap_CCR8.pdf) [TIGIT](Fig/STEP4//Fig_mtreg_gene_umap_TIGIT.pdf) [SATB1](Fig/STEP4//Fig_mtreg_gene_umap_SATB1.pdf) [ANXA2](Fig/STEP4//Fig_mtreg_gene_umap_ANXA2.pdf) [CCR4](Fig/STEP4//Fig_mtreg_gene_umap_CCR4.pdf) [CXCR3](Fig/STEP4//Fig_mtreg_gene_umap_CXCR3.pdf) [TNFRSF4](Fig/STEP4//Fig_mtreg_gene_umap_TNFRSF4.pdf) [PDCD1](Fig/STEP4//Fig_mtreg_gene_umap_PDCD1.pdf) [SELL](Fig/STEP4//Fig_mtreg_gene_umap_SELL.pdf) [FHIT](Fig/STEP4//Fig_mtreg_gene_umap_FHIT.pdf) [HLA-DRB1](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DRB1.pdf) [S100A4](Fig/STEP4//Fig_mtreg_gene_umap_S100A4.pdf) [RPS26](Fig/STEP4//Fig_mtreg_gene_umap_RPS26.pdf) [CD3E](Fig/STEP4//Fig_mtreg_gene_umap_CD3E.pdf) [HLA-DRA](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DRA.pdf) [HLA-C](Fig/STEP4//Fig_mtreg_gene_umap_HLA-C.pdf) [HLA-E](Fig/STEP4//Fig_mtreg_gene_umap_HLA-E.pdf) [HLA-B](Fig/STEP4//Fig_mtreg_gene_umap_HLA-B.pdf) [APOBEC3G](Fig/STEP4//Fig_mtreg_gene_umap_APOBEC3G.pdf) [MALAT1](Fig/STEP4//Fig_mtreg_gene_umap_MALAT1.pdf) [HLA-DR](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DR.pdf) 

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

    .file <- fig.dir %&% "/Fig_mtreg_gene_tsne_" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5, cat.link = F)
    cat("[" %&% g %&% "](" %&% .file %&% ") ")
}
```

[CD14](Fig/STEP4//Fig_mtreg_gene_tsne_CD14.pdf) [CD183](Fig/STEP4//Fig_mtreg_gene_tsne_CD183.pdf) [CD184](Fig/STEP4//Fig_mtreg_gene_tsne_CD184.pdf) [CD185](Fig/STEP4//Fig_mtreg_gene_tsne_CD185.pdf) [CD194](Fig/STEP4//Fig_mtreg_gene_tsne_CD194.pdf) [CD195](Fig/STEP4//Fig_mtreg_gene_tsne_CD195.pdf) [CD196](Fig/STEP4//Fig_mtreg_gene_tsne_CD196.pdf) [CD226](Fig/STEP4//Fig_mtreg_gene_tsne_CD226.pdf) [CD25](Fig/STEP4//Fig_mtreg_gene_tsne_CD25.pdf) [CD278](Fig/STEP4//Fig_mtreg_gene_tsne_CD278.pdf) [CD279](Fig/STEP4//Fig_mtreg_gene_tsne_CD279.pdf) [CD366](Fig/STEP4//Fig_mtreg_gene_tsne_CD366.pdf) [CD3](Fig/STEP4//Fig_mtreg_gene_tsne_CD3.pdf) [IL32](Fig/STEP4//Fig_mtreg_gene_tsne_IL32.pdf) [TRAF3IP3](Fig/STEP4//Fig_mtreg_gene_tsne_TRAF3IP3.pdf) [CD6](Fig/STEP4//Fig_mtreg_gene_tsne_CD6.pdf) [CD74](Fig/STEP4//Fig_mtreg_gene_tsne_CD74.pdf) [FAS](Fig/STEP4//Fig_mtreg_gene_tsne_FAS.pdf) [BRD9](Fig/STEP4//Fig_mtreg_gene_tsne_BRD9.pdf) [IKZF2](Fig/STEP4//Fig_mtreg_gene_tsne_IKZF2.pdf) [FOXP3](Fig/STEP4//Fig_mtreg_gene_tsne_FOXP3.pdf) [SIRT2](Fig/STEP4//Fig_mtreg_gene_tsne_SIRT2.pdf) [TBX21](Fig/STEP4//Fig_mtreg_gene_tsne_TBX21.pdf) [FOSL2](Fig/STEP4//Fig_mtreg_gene_tsne_FOSL2.pdf) [KEAP1](Fig/STEP4//Fig_mtreg_gene_tsne_KEAP1.pdf) [TCF7](Fig/STEP4//Fig_mtreg_gene_tsne_TCF7.pdf) [LAG3](Fig/STEP4//Fig_mtreg_gene_tsne_LAG3.pdf) [LYZ](Fig/STEP4//Fig_mtreg_gene_tsne_LYZ.pdf) [CD40LG](Fig/STEP4//Fig_mtreg_gene_tsne_CD40LG.pdf) [CORO1A](Fig/STEP4//Fig_mtreg_gene_tsne_CORO1A.pdf) [CTSH](Fig/STEP4//Fig_mtreg_gene_tsne_CTSH.pdf) [GSDMD](Fig/STEP4//Fig_mtreg_gene_tsne_GSDMD.pdf) [GATA3](Fig/STEP4//Fig_mtreg_gene_tsne_GATA3.pdf) [KLRB1](Fig/STEP4//Fig_mtreg_gene_tsne_KLRB1.pdf) [BACH2](Fig/STEP4//Fig_mtreg_gene_tsne_BACH2.pdf) [CCR6](Fig/STEP4//Fig_mtreg_gene_tsne_CCR6.pdf) [GZMK](Fig/STEP4//Fig_mtreg_gene_tsne_GZMK.pdf) [FOXP1](Fig/STEP4//Fig_mtreg_gene_tsne_FOXP1.pdf) [HDAC1](Fig/STEP4//Fig_mtreg_gene_tsne_HDAC1.pdf) [CD2](Fig/STEP4//Fig_mtreg_gene_tsne_CD2.pdf) [SGK1](Fig/STEP4//Fig_mtreg_gene_tsne_SGK1.pdf) [MT2A](Fig/STEP4//Fig_mtreg_gene_tsne_MT2A.pdf) [S1PR4](Fig/STEP4//Fig_mtreg_gene_tsne_S1PR4.pdf) [CCR7](Fig/STEP4//Fig_mtreg_gene_tsne_CCR7.pdf) [NR1D1](Fig/STEP4//Fig_mtreg_gene_tsne_NR1D1.pdf) [KLF2](Fig/STEP4//Fig_mtreg_gene_tsne_KLF2.pdf) [MAP1S](Fig/STEP4//Fig_mtreg_gene_tsne_MAP1S.pdf) [GIMAP4](Fig/STEP4//Fig_mtreg_gene_tsne_GIMAP4.pdf) [IL2RA](Fig/STEP4//Fig_mtreg_gene_tsne_IL2RA.pdf) [HAVCR2](Fig/STEP4//Fig_mtreg_gene_tsne_HAVCR2.pdf) [ITM2C](Fig/STEP4//Fig_mtreg_gene_tsne_ITM2C.pdf) [MYC](Fig/STEP4//Fig_mtreg_gene_tsne_MYC.pdf) [LEF1](Fig/STEP4//Fig_mtreg_gene_tsne_LEF1.pdf) [C1orf162](Fig/STEP4//Fig_mtreg_gene_tsne_C1orf162.pdf) [RORC](Fig/STEP4//Fig_mtreg_gene_tsne_RORC.pdf) [PSPH](Fig/STEP4//Fig_mtreg_gene_tsne_PSPH.pdf) [CDKN2A](Fig/STEP4//Fig_mtreg_gene_tsne_CDKN2A.pdf) [FLI1](Fig/STEP4//Fig_mtreg_gene_tsne_FLI1.pdf) [BATF](Fig/STEP4//Fig_mtreg_gene_tsne_BATF.pdf) [TSC22D3](Fig/STEP4//Fig_mtreg_gene_tsne_TSC22D3.pdf) [TNFRSF14](Fig/STEP4//Fig_mtreg_gene_tsne_TNFRSF14.pdf) [IFNGR2](Fig/STEP4//Fig_mtreg_gene_tsne_IFNGR2.pdf) [CXCR5](Fig/STEP4//Fig_mtreg_gene_tsne_CXCR5.pdf) [LMNA](Fig/STEP4//Fig_mtreg_gene_tsne_LMNA.pdf) [GBP4](Fig/STEP4//Fig_mtreg_gene_tsne_GBP4.pdf) [CTLA4](Fig/STEP4//Fig_mtreg_gene_tsne_CTLA4.pdf) [HPGD](Fig/STEP4//Fig_mtreg_gene_tsne_HPGD.pdf) [CDCA7L](Fig/STEP4//Fig_mtreg_gene_tsne_CDCA7L.pdf) [ABCA1](Fig/STEP4//Fig_mtreg_gene_tsne_ABCA1.pdf) [B2M](Fig/STEP4//Fig_mtreg_gene_tsne_B2M.pdf) [IRF2](Fig/STEP4//Fig_mtreg_gene_tsne_IRF2.pdf) [IL7R](Fig/STEP4//Fig_mtreg_gene_tsne_IL7R.pdf) [GPR25](Fig/STEP4//Fig_mtreg_gene_tsne_GPR25.pdf) [ISG20](Fig/STEP4//Fig_mtreg_gene_tsne_ISG20.pdf) [JUN](Fig/STEP4//Fig_mtreg_gene_tsne_JUN.pdf) [CD28](Fig/STEP4//Fig_mtreg_gene_tsne_CD28.pdf) [CCR8](Fig/STEP4//Fig_mtreg_gene_tsne_CCR8.pdf) [TIGIT](Fig/STEP4//Fig_mtreg_gene_tsne_TIGIT.pdf) [SATB1](Fig/STEP4//Fig_mtreg_gene_tsne_SATB1.pdf) [ANXA2](Fig/STEP4//Fig_mtreg_gene_tsne_ANXA2.pdf) [CCR4](Fig/STEP4//Fig_mtreg_gene_tsne_CCR4.pdf) [CXCR3](Fig/STEP4//Fig_mtreg_gene_tsne_CXCR3.pdf) [TNFRSF4](Fig/STEP4//Fig_mtreg_gene_tsne_TNFRSF4.pdf) [PDCD1](Fig/STEP4//Fig_mtreg_gene_tsne_PDCD1.pdf) [SELL](Fig/STEP4//Fig_mtreg_gene_tsne_SELL.pdf) [FHIT](Fig/STEP4//Fig_mtreg_gene_tsne_FHIT.pdf) [HLA-DRB1](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-DRB1.pdf) [S100A4](Fig/STEP4//Fig_mtreg_gene_tsne_S100A4.pdf) [RPS26](Fig/STEP4//Fig_mtreg_gene_tsne_RPS26.pdf) [CD3E](Fig/STEP4//Fig_mtreg_gene_tsne_CD3E.pdf) [HLA-DRA](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-DRA.pdf) [HLA-C](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-C.pdf) [HLA-E](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-E.pdf) [HLA-B](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-B.pdf) [APOBEC3G](Fig/STEP4//Fig_mtreg_gene_tsne_APOBEC3G.pdf) [MALAT1](Fig/STEP4//Fig_mtreg_gene_tsne_MALAT1.pdf) [HLA-DR](Fig/STEP4//Fig_mtreg_gene_tsne_HLA-DR.pdf) 

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

