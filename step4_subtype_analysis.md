---
title: "Step 4: subtype analysis"
author: Yongjin Park
theme: jekyll-theme-minimal
date: "2024-01-29"
bibliography: "MS_Ref.bib"
output:
  html_document:
    toc: true
    keep_md: true
    self_contained: true
---





```r
parse.tag <- function(x) {
    x[, c("barcode", "batch") := tstrsplit(tag, split="[_]")];
    x[, barcode := gsub(`barcode`, pattern="-[0-9]$", replacement="")];
    x[, batch := as.integer(`batch`)]
    x;
}
parse.gene <- function(x){
    x[, c("ensembl_gene_id","hgnc_symbol") := tstrsplit(`gene`,split="_")]
    x[, i := 1:.N]
    return(x)
}

bbknn.x <- function(.data, .bbknn, .bbknn.umap, .subset = NULL) {
    .cols <-
        fread(.data$col, header=F, col.names = "tag") %>%
        parse.tag() %>%
        mutate(j = 1:n()) %>%
        left_join(.bbknn.umap) %>%
        na.omit() %>%
        as.data.table()

    .rows <-
        fread(.data$row, col.names="gene", header=F) %>%
        parse.gene()

    if(!is.null(.subset)){
        .rows <- .rows %>%
            filter(hgnc_symbol %in% .subset |
                   ensembl_gene_id %in% .subset)
    }

    vv <- .bbknn$factors.adjusted[.cols$j, , drop = F]
    X <- sweep(.bbknn$U, 2, .bbknn$D, `*`) %*% t(vv)

    ## double scale
    X <- apply(X, 2, scale)
    X.t <- apply(t(X), 2, scale)
    X <- t(X.t)

    ret <- X[.rows$i, , drop = F]
    rownames(ret) <- .rows$hgnc_symbol
    colnames(ret) <- .cols$tag
    return(ret)
}

bbknn.x.melt <- function(.data, .bbknn, .bbknn.umap, .subset = NULL) {

    X <- bbknn.x(.data, .bbknn, .bbknn.umap, .subset)
    ret <- reshape2::melt(X) %>% as.data.table()
    ret[, tag := Var2]
    ret[, gene := Var1]
    ret[, Var1 := NULL]
    ret[, Var2 := NULL]
    ret[, x := value]
    ret[, value := NULL]
    return(ret)
}

take.marker.stats <- function(x, dict){

    require(matrixTests)

    ret <- data.table()
    for(k in unique(dict$membership)){
        tag.k <- .bbknn.umap[membership == k]$tag
        x.k <- x[, colnames(x) %in% tag.k, drop = F]
        y.k <- x[, !(colnames(x) %in% tag.k), drop = F]
        ## wilcox within vs. outside
        wilcox.k <- row_wilcoxon_twosample(x.k, y.k)
        ## t-test
        ttest.k <- row_t_welch(x.k, y.k)
        ret.k <- data.table(
            gene = rownames(x),
            membership = k,
            wilcox.p = wilcox.k$pvalue,
            ttest.t = ttest.k$statistic,
            mean = ttest.k$mean.x,
            mean.outside = ttest.k$mean.y,
            sd = sqrt(ttest.k$var.x),
            sd.outside = sqrt(ttest.k$var.y),
            ttest.p = ttest.k$pvalue)
        ret <- rbind(ret, ret.k)
    }
    return(ret)
}
```

# Results

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
.sample.info <-
    readxl::read_xlsx("data/Hashing list MS Treg project.xlsx", 1) %>%
    mutate(disease=substr(`subject`,1,2)) %>%
    mutate(hash = gsub(pattern="#", replacement="", `hash`)) %>%
    mutate(hash = as.integer(hash)) %>%
    as.data.table()
```




```r
annot.dt <- fread("Tab/step2_celltype.txt.gz")
```

## 1. Memory T conventional


```r
.full.data <- fileset.list("result/step1/matrix_final")
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

### Clustering cells by batch-balancing k-nearest neighbour graph


```r
.file <- "result/step4/mtconv_bbknn.rds"
if.needed(.file, {

    .batches <-
        fread(.data$col, col.names="tag", header=F) %>%
        parse.tag() %>%
        left_join(.hash.info)

    bb <- .batches %>%
        mutate(b = batch %&% "_" %&% disease) %>%
        select(b) %>%
        unlist()

    .bbknn <- rcpp_mmutil_bbknn_svd(.data$mtx, bb,
                                    knn = 30, RANK = 50,
                                    EM_ITER = 20,
                                    TAKE_LN = TRUE,
                                    NUM_THREADS = 16)
    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```


```r
.file <- "result/step4/mtconv_bbknn_leiden.txt.gz"
if.needed(.file, {
    .cells <- readLines(.data$col)
    .bbknn.umap <- run.bbknn.umap(.bbknn$knn,
                                  .cells,
                                  symmetrize=T,
                                  min_dist=0.01,
                                  spread=5,
                                  res=.5)
    fwrite(.bbknn.umap, .file)
})

.bbknn.umap <- fread(.file)
.bbknn.umap[, membership := as.factor(membership)]
```



[**DOWNLOAD:** mTconv clustering results](Tab/step4_mtconv_cluster.txt.gz)

### What are the cell-cluster-specific marker genes?


```r
.mkdir("Tab/")
.file <- "Tab/step4_mtconv_gene_stat.txt.gz"
if.needed(.file, {
    x <- bbknn.x(.data, .bbknn, .bbknn.umap)
    marker.stat <- take.marker.stats(x, .bbknn.umap)
    fwrite(marker.stat, .file, sep = "\t", col.names = T)
})
marker.stat <- fread(.file, sep = "\t")
```

[**DOWNLOAD:** mTconv marker gene statistics](Tab/step4_mtconv_gene_stat.txt.gz)

### Non-linear embedding to confirm the cell clusters of mTconv cells


```r
.lab <-
    .cells[,
           .(umap1=median(umap1), umap2=median(umap2)),
           by = .(membership)]
.cols <- .more.colors(nrow(.lab), nc.pal=12)
plt <-
    .gg.plot(.bbknn.umap, aes(umap1, umap2, color=membership)) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtconv-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtconv.pdf)


```r
.cols <- .more.colors(10, nc.pal=7, .palette="Set1")
plt <-
    .gg.plot(.cells, aes(umap1, umap2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")

print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtconv_sub-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_bbknn_mtconv_sub.pdf)


```r
x.melt <- bbknn.x.melt(.data, .bbknn, .bbknn.umap, .markers)
```

#### Summary heatmap


```r
.dt <- x.melt %>% left_join(.cells)
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
.marker.order <- sort(unique(.sum$col))
```


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

**NOTE** The colors are standardized `log1p` expression across genes and cells.


```r
for(g in unique(x.melt$gene)) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(umap1, umap2, color=pmax(pmin(x, 3), -3))

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

### Basic statistics




```r
.stat <-
    .cells[,
           .(N = .N),
           by=.(batch, membership, disease)] %>%
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
    .sum.stat.tot %>%
    mutate(batch = "(N=" %&% num.int(sum(.stat$N)) %&% ")")

plt <- .plt.sum.stat(.stat.tot) + ggtitle("mTconv")
print(plt)
```

![](Fig/STEP4/Fig_count_merged_mtconv_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_merged_mtconv_tot.pdf)

## 2. Memory Treg cells


```r
.full.data <- fileset.list("result/step1/matrix_final")
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

### Clustering cells by batch-balancing k-nearest neighbour graph


```r
.file <- "result/step4/mtreg_bbknn.rds"
if.needed(.file, {

    .batches <-
        fread(.data$col, col.names="tag", header=F) %>%
        parse.tag() %>%
        left_join(.hash.info)

    bb <- .batches %>%
        mutate(b = batch %&% "_" %&% disease) %>%
        select(b) %>%
        unlist()

    .bbknn <- rcpp_mmutil_bbknn_svd(.data$mtx, bb,
                                    knn = 30, RANK = 50,
                                    EM_ITER = 20,
                                    TAKE_LN = TRUE,
                                    NUM_THREADS = 16)
    saveRDS(.bbknn, .file)
})
.bbknn <- readRDS(.file)
```



```r
.file <- "result/step4/mtreg_bbknn_leiden.txt.gz"
.cells <- readLines(.data$col)

if.needed(.file, {
    .bbknn.umap <- run.bbknn.umap(.bbknn$knn,
                                  .cells,
                                  symmetrize=T,
                                  min_dist=0.01,
                                  spread=5,
                                  res=.5)
    fwrite(.bbknn.umap, .file)
})

.bbknn.umap <- fread(.file)
.bbknn.umap[, membership := as.factor(membership)]
```



### Merge clusters 3 and 6

We merged the clusters \#3 and \#6 since they show similar gene expression patterns for known marker genes. Merge cells in the cluster \#6 to \#3.


```r
.bbknn.umap[membership == 6, membership := 3]
.bbknn.umap[, membership := factor(membership, 1:5)]
.cells[membership == 6, membership := 3]
.cells[, membership := factor(membership, 1:5)]
```

[**DOWNLOAD:** mTreg clustering results](Tab/step4_mtreg_cluster.txt.gz)

### What are the cell-cluster-specific marker genes?


```r
.mkdir("Tab/")
.file <- "Tab/step4_mtreg_gene_stat.txt.gz"
if.needed(.file, {
    x <- bbknn.x(.data, .bbknn, .bbknn.umap)
    marker.stat <- take.marker.stats(x, .bbknn.umap)
    fwrite(marker.stat, .file, sep = "\t", col.names = T)
})
marker.stat <- fread(.file, sep = "\t")
```

[**DOWNLOAD:** mTreg marker gene statistics](Tab/step4_mtreg_gene_stat.txt.gz)

### Non-linear embedding to confirm the cell clusters of mTreg cells


```r
.lab <-
    .cells[,
           .(umap1=median(umap1), umap2=median(umap2)),
           by = .(membership)]

.cols <- .more.colors(nrow(.lab), nc.pal=12)

plt <-
    .gg.plot(.bbknn.umap, aes(umap1, umap2, color=membership)) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    geom_text(aes(label=membership), data=.lab, size=4, color="black") +
    scale_color_manual(values = .cols, guide="none")

print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtreg-1.png)<!-- -->

```r
.file <- fig.dir %&% "/Fig_bbknn_mtreg.pdf"
.gg.save(filename = .file, plot = plt, width=3, height=3)
```



[PDF](Fig/STEP4//Fig_bbknn_mtreg.pdf)

#### Subject-level variation is not the major source of this clustering pattern.


```r
.cols <- .more.colors(10, nc.pal=7, .palette="Set1")
plt <-
    .gg.plot(.cells, aes(umap1, umap2, color=as.factor(subject))) +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::rasterise(geom_point(stroke=0, alpha=.8, size=.7), dpi=300) +
    scale_color_manual(values = .cols, guide="none")
print(plt)
```

![](Fig/STEP4/Fig_bbknn_mtreg_sub-1.png)<!-- -->

```r
.file <- fig.dir %&% "/Fig_bbknn_mtreg_sub.pdf"
.gg.save(filename = .file, plot = plt, width=3, height=3)
```



[PDF](Fig/STEP4//Fig_bbknn_mtreg_sub.pdf)

#### Confirm cell type identity based on known marker genes?


```r
x.melt <- bbknn.x.melt(.data, .bbknn, .bbknn.umap, .markers)
```


```r
.dt <- x.melt %>% left_join(.cells)
.sum.subj <- .dt[, .(x = median(x)), by = .(gene, subject, membership)]
.sum.subj[, x := scale(x), by = .(gene)]
```

#### Summary heatmap


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

```r
.file <- fig.dir %&% "/Fig_mtreg_sum_membership.pdf"
.gg.save(filename = .file, plot = plt, width=2, height=6)
```



[PDF](Fig/STEP4//Fig_mtreg_sum_membership.pdf)


```r
.marker.order <- sort(unique(.sum$col))
```


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

```r
.file <- fig.dir %&% "/Fig_mtreg_sum_subj_memberhip.pdf"
.gg.save(filename = .file, plot = plt, width=6, height=6)
```



[PDF](Fig/STEP4//Fig_mtreg_sum_subj_memberhip.pdf)

**NOTE** The colors are standardized `log1p` expression across genes and cells.


```r
.genes <- unique(x.melt$gene)
for(g in .genes) {
    .dt <- left_join(x.melt[gene == g], .cells)
    .aes <- aes(umap1, umap2, color=pmax(pmin(x, 3), -3))

    plt <-
        .gg.plot(.dt[order(`x`)], .aes) +
        xlab("UMAP1") + ylab("UMAP2") +
        ggrastr::rasterise(geom_point(stroke = 0, size=.7), dpi=300) +
        theme(legend.key.width = unit(.2,"lines")) +
        theme(legend.key.height = unit(.5,"lines")) +
        scale_color_distiller(g, palette = "RdBu", direction = -1) +
        ggtitle(g)

    print(plt)
    .file <- fig.dir %&% "/Fig_mtreg_gene_umap_" %&% g %&% ".pdf"
    .gg.save(filename = .file, plot = plt, width=3, height=2.5)
}
```

![](Fig/STEP4/Fig_mtreg_gene_umap-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD183.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-3.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD184.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-4.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD185.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-5.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD194.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-6.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD195.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-7.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD196.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-8.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD226.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-9.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-10.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD278.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-11.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD279.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-12.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD366.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-13.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-14.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IL32.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-15.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TRAF3IP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-16.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-17.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD74.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-18.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FAS.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-19.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_BRD9.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-20.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IKZF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-21.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FOXP3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-22.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_SIRT2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-23.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TBX21.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-24.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FOSL2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-25.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_KEAP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-26.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TCF7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-27.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_LAG3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-28.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_LYZ.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-29.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD40LG.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-30.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CORO1A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-31.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CTSH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-32.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GSDMD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-33.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GATA3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-34.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_KLRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-35.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_BACH2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-36.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CCR6.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-37.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GZMK.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-38.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FOXP1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-39.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HDAC1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-40.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-41.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_SGK1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-42.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_MT2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-43.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_S1PR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-44.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CCR7.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-45.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_NR1D1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-46.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_KLF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-47.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_MAP1S.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-48.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GIMAP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-49.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IL2RA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-50.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HAVCR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-51.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_ITM2C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-52.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_MYC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-53.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_LEF1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-54.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_C1orf162.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-55.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_RORC.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-56.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_PSPH.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-57.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CDKN2A.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-58.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FLI1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-59.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_BATF.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-60.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TSC22D3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-61.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TNFRSF14.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-62.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IFNGR2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-63.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CXCR5.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-64.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_LMNA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-65.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GBP4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-66.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CTLA4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-67.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HPGD.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-68.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CDCA7L.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-69.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_ABCA1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-70.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_B2M.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-71.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IRF2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-72.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_IL7R.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-73.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_GPR25.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-74.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_ISG20.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-75.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_JUN.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-76.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD28.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-77.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CCR8.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-78.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TIGIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-79.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_SATB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-80.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_ANXA2.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-81.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CCR4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-82.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CXCR3.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-83.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_TNFRSF4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-84.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_PDCD1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-85.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_SELL.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-86.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_FHIT.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-87.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DRB1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-88.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_S100A4.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-89.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_RPS26.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-90.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_CD3E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-91.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DRA.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-92.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-C.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-93.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-E.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-94.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-B.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-95.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_APOBEC3G.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-96.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_MALAT1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap-97.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_HLA-DR.pdf)


```r
.stat <-
    .cells[,
           .(N = .N),
           by=.(batch, membership, disease)] %>%
    .sum.stat.batch()

plt <- .plt.sum.stat(.stat) + ggtitle("mTreg")
print(plt)
```

![](Fig/STEP4/Fig_count_mtreg_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_mtreg_tot.pdf)


```r
.stat.tot <-
    .cells[,
           .(N = .N),
           by=.(membership, disease)] %>%
    .sum.stat.tot %>%
    mutate(batch = "(N=" %&% num.int(sum(.stat$N)) %&% ")")

plt <- .plt.sum.stat(.stat.tot) + ggtitle("mTreg")
print(plt)
```

![](Fig/STEP4/Fig_count_merged_mtreg_tot-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_count_merged_mtreg_tot.pdf)


### PRDM1 vs. SGK1 (also standardized across cells)

![](Fig/STEP4/Fig_mtreg_gene_umap_other-1.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_PRDM1.pdf)

![](Fig/STEP4/Fig_mtreg_gene_umap_other-2.png)<!-- -->

[PDF](Fig/STEP4//Fig_mtreg_gene_umap_SGK1.pdf)

### What about the PRDM1 long and short isoforms?

![](Fig/STEP4/Fig_mtreg_PRDM1_iso_umap-1.png)<!-- -->


[PDF](Fig/STEP4//Fig_mtreg_gene_umap_prdm1_iso.pdf)

### Can we find new marker genes/proteins for mTreg clusters?


```r
## Actually, p-values are meaningless because we did clustering based
## on expression and testing the significance of clustering results.
.stat <-
    marker.stat[p.adjust(wilcox.p) < .01 &
                p.adjust(ttest.p) < .01 &
                ttest.t > 0]

.top.markers <-
    .stat[order(`ttest.t`, `sd`, decreasing = T),
          head(.SD, 10),
          by = .(membership)]

top.stat <- .stat[gene %in% .top.markers$gene]
```
