---
title: "Step 5: differential Expression analysis for sub clusters"
date: "2024-02-08"
author: Yongjin Park
bibliography: "MS_Ref.bib"
output:
  html_document:
    toc: true
    keep_md: true
    self_contained: true
---




```r
.hash.hdr <- "result/step1/hash"
.hash.data <- fileset.list(.hash.hdr)
.hash.info <- read.hash(.hash.data)
```



# 1. mTconv subtype DEG analysis


```r
annot.dt <-
    fread("Tab/step4_mtconv_leiden.txt.gz") %>%
    left_join(.hash.info) %>%
    na.omit()
```


```r
.data <- fileset.list("result/step1/matrix_final")
.mkdir("result/step5/deg/")
.deg.data <- fileset.list("result/step5/deg/mtconv_hc_ms")

if.needed(.deg.data, {
    .deg.data <-
        rcpp_mmutil_copy_selected_columns(.data$mtx,
                                          .data$row,
                                          .data$col,
                                          unique(annot.dt$tag),
                                          "result/step5/deg/mtconv_hc_ms")
})

.file <- "result/step5/deg/mtconv_hc_ms.rds"

if.needed(.file, {

    .membership <-
        annot.dt[, .(tag, membership)] %>%
        as.data.frame()

    .cell2indv <- annot.dt[, .(tag, subject)] %>%
        unique %>%
        as.data.frame()

    .indv2exp <- .cell2indv %>%
        select(subject) %>%
        mutate(disease = substr(`subject`, 1, 2)) %>%
        as.data.frame()

    .deg.stat <-
        make.cocoa(.deg.data, .membership, .cell2indv, .indv2exp,
                   .rank = 10, .em.iter = 20, .em.tol = 1e-8, .take.ln = TRUE,
                   knn = 30, impute.by.knn = TRUE, num.threads = 8)

   saveRDS(.deg.stat, .file)
})

.deg.stat <- readRDS(.file)

.mems <- unique(annot.dt$membership)
.indvs <- unique(annot.dt$subject)

.hc.ms.dt <-
    list(tot = sort.col(.deg.stat$sum, .mems, .indvs),
         cfa = sort.col(.deg.stat$resid.ln.mu, .mems, .indvs),
         cfa.sd = sort.col(.deg.stat$resid.ln.mu.sd, .mems, .indvs)) %>%
    combine.statistics()

.hc.ms.dt[, c("subject","membership") := tstrsplit(`Var2`,split="_")]
.hc.ms.dt[, disease := substr(`subject`,1,2)]
.hc.ms.dt[, gene := `Var1`]
.hc.ms.dt[, `Var1` := NULL]
.hc.ms.dt[, `Var2` := NULL]

hc.ms.deg <- summarize.deg(.hc.ms.dt) %>% parse.gene()
```

[**DOWNLOAD:** mTconv DEG MS vs HC](Tab/DEG_mtconv_MS_vs_HC.txt.gz)

### Found 755 unique genes strongly perturbed by MS with FDR 5%

* Up-regulated: 499

* Down-regulated:  288

* Total pairs of genes and clusters: 65,226


```r
count.deg <- function(.dt, fdr.cutoff = .05) {
    .dt[fdr < fdr.cutoff &
        sign(ADD) == sign(ADE) &
        sign(ADC) == sign(ADE),
        .(n = .N),
        by = .(membership,
               direction = if_else(z > 0, "up", "down"))
        ] %>%
        mutate(direction = factor(direction, c("up", "down"))) %>%
        group_by(membership) %>%
        arrange(desc(direction)) %>%
        mutate(nc = cumsum(n)) %>%
        ungroup
}
```

![](Fig/STEP5/Fig_mTconv_DEG_count-1.png)<!-- -->

#### Examples




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


```r
.bulk.genes <-
    bulk.mTconv.dt[bulk.qv < .2, .(hgnc_symbol)] %>% 
    unique() %>%
    unlist()

hc.ms.deg[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`gene`, split="_")]

.genes.show <-
    hc.ms.deg[fdr < 0.05 &
              hgnc_symbol %in% .bulk.genes &
              sign(ADD) == sign(ADE) &
              sign(ADC) == sign(ADE)] %>%
    select(gene) %>%
    unique %>%
    unlist %>%
    as.character
```

[CD226_CD226](Fig/STEP5//Fig_mTconv_DEG_example_CD226_CD226.pdf) [CD27_CD27](Fig/STEP5//Fig_mTconv_DEG_example_CD27_CD27.pdf) [CD58_CD58](Fig/STEP5//Fig_mTconv_DEG_example_CD58_CD58.pdf) [ENSG00000009790_TRAF3IP3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000009790_TRAF3IP3.pdf) [ENSG00000034677_RNF19A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000034677_RNF19A.pdf) [ENSG00000038274_MAT2B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000038274_MAT2B.pdf) [ENSG00000057657_PRDM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000057657_PRDM1.pdf) [ENSG00000058668_ATP2B4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000058668_ATP2B4.pdf) [ENSG00000065978_YBX1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000065978_YBX1.pdf) [ENSG00000070756_PABPC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000070756_PABPC1.pdf) [ENSG00000071073_MGAT4A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000071073_MGAT4A.pdf) [ENSG00000073849_ST6GAL1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000073849_ST6GAL1.pdf) [ENSG00000076043_REXO2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000076043_REXO2.pdf) [ENSG00000077420_APBB1IP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000077420_APBB1IP.pdf) [ENSG00000081237_PTPRC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000081237_PTPRC.pdf) [ENSG00000083799_CYLD](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000083799_CYLD.pdf) [ENSG00000085491_SLC25A24](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000085491_SLC25A24.pdf) [ENSG00000099622_CIRBP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000099622_CIRBP.pdf) [ENSG00000099783_HNRNPM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000099783_HNRNPM.pdf) [ENSG00000099985_OSM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000099985_OSM.pdf) [ENSG00000100100_PIK3IP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100100_PIK3IP1.pdf) [ENSG00000100129_EIF3L](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100129_EIF3L.pdf) [ENSG00000100219_XBP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100219_XBP1.pdf) [ENSG00000100225_FBXO7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100225_FBXO7.pdf) [ENSG00000100316_RPL3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100316_RPL3.pdf) [ENSG00000100644_HIF1A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100644_HIF1A.pdf) [ENSG00000100664_EIF5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100664_EIF5.pdf) [ENSG00000101166_PRELID3B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000101166_PRELID3B.pdf) [ENSG00000101608_MYL12A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000101608_MYL12A.pdf) [ENSG00000101654_RNMT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000101654_RNMT.pdf) [ENSG00000102007_PLP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102007_PLP2.pdf) [ENSG00000102245_CD40LG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102245_CD40LG.pdf) [ENSG00000102409_BEX4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102409_BEX4.pdf) [ENSG00000102871_TRADD](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102871_TRADD.pdf) [ENSG00000104529_EEF1D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104529_EEF1D.pdf) [ENSG00000104660_LEPROTL1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104660_LEPROTL1.pdf) [ENSG00000106460_TMEM106B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000106460_TMEM106B.pdf) [ENSG00000108298_RPL19](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000108298_RPL19.pdf) [ENSG00000109861_CTSC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109861_CTSC.pdf) [ENSG00000109971_HSPA8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109971_HSPA8.pdf) [ENSG00000110852_CLEC2B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000110852_CLEC2B.pdf) [ENSG00000112110_MRPL18](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000112110_MRPL18.pdf) [ENSG00000113407_TARS](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000113407_TARS.pdf) [ENSG00000114737_CISH](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000114737_CISH.pdf) [ENSG00000115687_PASK](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000115687_PASK.pdf) [ENSG00000115758_ODC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000115758_ODC1.pdf) [ENSG00000116717_GADD45A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000116717_GADD45A.pdf) [ENSG00000116791_CRYZ](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000116791_CRYZ.pdf) [ENSG00000116906_GNPAT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000116906_GNPAT.pdf) [ENSG00000117318_ID3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000117318_ID3.pdf) [ENSG00000118515_SGK1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000118515_SGK1.pdf) [ENSG00000119718_EIF2B2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000119718_EIF2B2.pdf) [ENSG00000120656_TAF12](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000120656_TAF12.pdf) [ENSG00000120875_DUSP4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000120875_DUSP4.pdf) [ENSG00000120913_PDLIM2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000120913_PDLIM2.pdf) [ENSG00000121067_SPOP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000121067_SPOP.pdf) [ENSG00000121895_TMEM156](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000121895_TMEM156.pdf) [ENSG00000122026_RPL21](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122026_RPL21.pdf) [ENSG00000122188_LAX1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122188_LAX1.pdf) [ENSG00000122406_RPL5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122406_RPL5.pdf) [ENSG00000123609_NMI](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000123609_NMI.pdf) [ENSG00000124177_CHD6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000124177_CHD6.pdf) [ENSG00000125148_MT2A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000125148_MT2A.pdf) [ENSG00000128340_RAC2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000128340_RAC2.pdf) [ENSG00000129315_CCNT1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000129315_CCNT1.pdf) [ENSG00000131143_COX4I1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000131143_COX4I1.pdf) [ENSG00000132406_TMEM128](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132406_TMEM128.pdf) [ENSG00000132475_H3F3B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132475_H3F3B.pdf) [ENSG00000132965_ALOX5AP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132965_ALOX5AP.pdf) [ENSG00000133561_GIMAP6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133561_GIMAP6.pdf) [ENSG00000133574_GIMAP4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133574_GIMAP4.pdf) [ENSG00000133872_SARAF](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133872_SARAF.pdf) [ENSG00000134308_YWHAQ](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000134308_YWHAQ.pdf) [ENSG00000134954_ETS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000134954_ETS1.pdf) [ENSG00000135046_ANXA1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000135046_ANXA1.pdf) [ENSG00000135486_HNRNPA1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000135486_HNRNPA1.pdf) [ENSG00000136732_GYPC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000136732_GYPC.pdf) [ENSG00000136819_C9orf78](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000136819_C9orf78.pdf) [ENSG00000136827_TOR1A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000136827_TOR1A.pdf) [ENSG00000137154_RPS6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000137154_RPS6.pdf) [ENSG00000137876_RSL24D1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000137876_RSL24D1.pdf) [ENSG00000138166_DUSP5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000138166_DUSP5.pdf) [ENSG00000138433_CIR1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000138433_CIR1.pdf) [ENSG00000138757_G3BP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000138757_G3BP2.pdf) [ENSG00000139187_KLRG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000139187_KLRG1.pdf) [ENSG00000139211_AMIGO2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000139211_AMIGO2.pdf) [ENSG00000140988_RPS2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000140988_RPS2.pdf) [ENSG00000143183_TMCO1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143183_TMCO1.pdf) [ENSG00000143333_RGS16](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143333_RGS16.pdf) [ENSG00000143570_SLC39A1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143570_SLC39A1.pdf) [ENSG00000143772_ITPKB](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143772_ITPKB.pdf) [ENSG00000143891_GALM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143891_GALM.pdf) [ENSG00000144746_ARL6IP5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000144746_ARL6IP5.pdf) [ENSG00000145741_BTF3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145741_BTF3.pdf) [ENSG00000145779_TNFAIP8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145779_TNFAIP8.pdf) [ENSG00000145860_RNF145](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145860_RNF145.pdf) [ENSG00000146457_WTAP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000146457_WTAP.pdf) [ENSG00000147403_RPL10](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000147403_RPL10.pdf) [ENSG00000147604_RPL7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000147604_RPL7.pdf) [ENSG00000148154_UGCG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148154_UGCG.pdf) [ENSG00000148303_RPL7A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148303_RPL7A.pdf) [ENSG00000148634_HERC4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148634_HERC4.pdf) [ENSG00000148834_GSTO1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148834_GSTO1.pdf) [ENSG00000149273_RPS3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000149273_RPS3.pdf) [ENSG00000150637_CD226](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000150637_CD226.pdf) [ENSG00000152518_ZFP36L2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000152518_ZFP36L2.pdf) [ENSG00000156508_EEF1A1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000156508_EEF1A1.pdf) [ENSG00000160593_JAML](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000160593_JAML.pdf) [ENSG00000161016_RPL8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000161016_RPL8.pdf) [ENSG00000162517_PEF1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162517_PEF1.pdf) [ENSG00000162704_ARPC5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162704_ARPC5.pdf) [ENSG00000162777_DENND2D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162777_DENND2D.pdf) [ENSG00000163191_S100A11](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163191_S100A11.pdf) [ENSG00000163682_RPL9](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163682_RPL9.pdf) [ENSG00000163939_PBRM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163939_PBRM1.pdf) [ENSG00000164024_METAP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000164024_METAP1.pdf) [ENSG00000164104_HMGB2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000164104_HMGB2.pdf) [ENSG00000164442_CITED2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000164442_CITED2.pdf) [ENSG00000166913_YWHAB](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000166913_YWHAB.pdf) [ENSG00000167658_EEF2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000167658_EEF2.pdf) [ENSG00000167996_FTH1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000167996_FTH1.pdf) [ENSG00000168028_RPSA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168028_RPSA.pdf) [ENSG00000168209_DDIT4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168209_DDIT4.pdf) [ENSG00000168300_PCMTD1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168300_PCMTD1.pdf) [ENSG00000168421_RHOH](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168421_RHOH.pdf) [ENSG00000168685_IL7R](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168685_IL7R.pdf) [ENSG00000169100_SLC25A6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000169100_SLC25A6.pdf) [ENSG00000169239_CA5B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000169239_CA5B.pdf) [ENSG00000170430_MGMT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000170430_MGMT.pdf) [ENSG00000171566_PLRG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171566_PLRG1.pdf) [ENSG00000171792_RHNO1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171792_RHNO1.pdf) [ENSG00000171863_RPS7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171863_RPS7.pdf) [ENSG00000173726_TOMM20](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000173726_TOMM20.pdf) [ENSG00000173812_EIF1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000173812_EIF1.pdf) [ENSG00000174444_RPL4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174444_RPL4.pdf) [ENSG00000174500_GCSAM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174500_GCSAM.pdf) [ENSG00000174720_LARP7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174720_LARP7.pdf) [ENSG00000174748_RPL15](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174748_RPL15.pdf) [ENSG00000174851_YIF1A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174851_YIF1A.pdf) [ENSG00000175390_EIF3F](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000175390_EIF3F.pdf) [ENSG00000177917_ARL6IP6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000177917_ARL6IP6.pdf) [ENSG00000178149_DALRD3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000178149_DALRD3.pdf) [ENSG00000179144_GIMAP7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000179144_GIMAP7.pdf) [ENSG00000179820_MYADM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000179820_MYADM.pdf) [ENSG00000180611_MB21D2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000180611_MB21D2.pdf) [ENSG00000181163_NPM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000181163_NPM1.pdf) [ENSG00000182472_CAPN12](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000182472_CAPN12.pdf) [ENSG00000182718_ANXA2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000182718_ANXA2.pdf) [ENSG00000183813_CCR4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000183813_CCR4.pdf) [ENSG00000184009_ACTG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184009_ACTG1.pdf) [ENSG00000184508_HDDC3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184508_HDDC3.pdf) [ENSG00000184588_PDE4B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184588_PDE4B.pdf) [ENSG00000185338_SOCS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000185338_SOCS1.pdf) [ENSG00000187514_PTMA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000187514_PTMA.pdf) [ENSG00000188846_RPL14](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000188846_RPL14.pdf) [ENSG00000189077_TMEM120A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000189077_TMEM120A.pdf) [ENSG00000196329_GIMAP5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000196329_GIMAP5.pdf) [ENSG00000196531_NACA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000196531_NACA.pdf) [ENSG00000197111_PCBP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197111_PCBP2.pdf) [ENSG00000197329_PELI1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197329_PELI1.pdf) [ENSG00000198242_RPL23A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198242_RPL23A.pdf) [ENSG00000198763_MT-ND2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198763_MT-ND2.pdf) [ENSG00000198830_HMGN2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198830_HMGN2.pdf) [ENSG00000198938_MT-CO3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198938_MT-CO3.pdf) [ENSG00000204628_RACK1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000204628_RACK1.pdf) [ENSG00000211694_TRGV10](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211694_TRGV10.pdf) [ENSG00000211750_TRBV24-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211750_TRBV24-1.pdf) [ENSG00000211786_TRAV8-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211786_TRAV8-2.pdf) [ENSG00000211787_TRAV8-3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211787_TRAV8-3.pdf) [ENSG00000211788_TRAV13-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211788_TRAV13-1.pdf) [ENSG00000211789_TRAV12-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211789_TRAV12-2.pdf) [ENSG00000211790_TRAV8-4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211790_TRAV8-4.pdf) [ENSG00000211797_TRAV17](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211797_TRAV17.pdf) [ENSG00000211801_TRAV21](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211801_TRAV21.pdf) [ENSG00000211803_TRAV23DV6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211803_TRAV23DV6.pdf) [ENSG00000211806_TRAV25](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211806_TRAV25.pdf) [ENSG00000211807_TRAV26-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211807_TRAV26-1.pdf) [ENSG00000211809_TRAV27](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211809_TRAV27.pdf) [ENSG00000211810_TRAV29DV5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211810_TRAV29DV5.pdf) [ENSG00000211815_TRAV36DV7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211815_TRAV36DV7.pdf) [ENSG00000211817_TRAV38-2DV8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211817_TRAV38-2DV8.pdf) [ENSG00000213145_CRIP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213145_CRIP1.pdf) [ENSG00000213402_PTPRCAP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213402_PTPRCAP.pdf) [ENSG00000213626_LBH](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213626_LBH.pdf) [ENSG00000213719_CLIC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213719_CLIC1.pdf) [ENSG00000224032_EPB41L4A-AS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000224032_EPB41L4A-AS1.pdf) [ENSG00000230124_ACBD6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000230124_ACBD6.pdf) [ENSG00000241657_TRBV11-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000241657_TRBV11-2.pdf) [ENSG00000243927_MRPS6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000243927_MRPS6.pdf) [ENSG00000274752_TRBV12-3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000274752_TRBV12-3.pdf) [ENSG00000277734_TRAC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000277734_TRAC.pdf) [ENSG00000278030_TRBV7-9](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000278030_TRBV7-9.pdf) [TIGIT_TIGIT](Fig/STEP5//Fig_mTconv_DEG_example_TIGIT_TIGIT.pdf) 

# 2. mTreg subtype DEG analysis


```r
annot.dt <-
    fread("Tab/step4_mtreg_leiden.txt.gz") %>%
    left_join(.hash.info) %>%
    na.omit()
```


```r
.data <- fileset.list("result/step1/matrix_final")
.mkdir("result/step5/deg/")
.deg.data <- fileset.list("result/step5/deg/mtreg_hc_ms")

if.needed(.deg.data, {
    .deg.data <-
        rcpp_mmutil_copy_selected_columns(.data$mtx,
                                          .data$row,
                                          .data$col,
                                          unique(annot.dt$tag),
                                          "result/step5/deg/mtreg_hc_ms")
})

.file <- "result/step5/deg/mtreg_hc_ms.rds"

if.needed(.file, {

    .membership <-
        annot.dt[, .(tag, membership)] %>%
        as.data.frame()

    .cell2indv <- annot.dt[, .(tag, subject)] %>%
        unique %>%
        as.data.frame()

    .indv2exp <- .cell2indv %>%
        select(subject) %>%
        mutate(disease = substr(`subject`, 1, 2)) %>%
        as.data.frame()

    .deg.stat <-
        make.cocoa(.deg.data, .membership, .cell2indv, .indv2exp,
                   .rank = 10, .em.iter = 20, .em.tol = 1e-8, .take.ln = TRUE,
                   knn = 30, impute.by.knn = TRUE, num.threads = 8)

   saveRDS(.deg.stat, .file)
})

.deg.stat <- readRDS(.file)

.mems <- unique(annot.dt$membership)
.indvs <- unique(annot.dt$subject)

.hc.ms.dt <-
    list(tot = sort.col(.deg.stat$sum, .mems, .indvs),
         cfa = sort.col(.deg.stat$resid.ln.mu, .mems, .indvs),
         cfa.sd = sort.col(.deg.stat$resid.ln.mu.sd, .mems, .indvs)) %>%
    combine.statistics()

.hc.ms.dt[, c("subject","membership") := tstrsplit(`Var2`,split="_")]
.hc.ms.dt[, disease := substr(`subject`,1,2)]
.hc.ms.dt[, gene := `Var1`]
.hc.ms.dt[, `Var1` := NULL]
.hc.ms.dt[, `Var2` := NULL]

hc.ms.deg <- summarize.deg(.hc.ms.dt)
hc.ms.deg[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`gene`, split="_")]
```

[**DOWNLOAD:** mtreg DEG MS vs HC](Tab/DEG_mtreg_MS_vs_HC.txt.gz)

### Found 959 unique genes strongly perturbed by MS with FDR 5%

* Up-regulated: 292

* Down-regulated:  731

* Total pairs of genes and clusters: 94,900

![](Fig/STEP5/Fig_mTreg_DEG_count-1.png)<!-- -->


```r
.bulk.genes <-
    bulk.mTreg.dt[bulk.qv < .2, .(hgnc_symbol)] %>% 
    unique() %>%
    unlist()

hc.ms.deg[, c("ensembl_gene_id", "hgnc_symbol") := tstrsplit(`gene`, split="_")]

.genes.show <-
    hc.ms.deg[fdr < 0.05 &
              hgnc_symbol %in% .bulk.genes &
              sign(ADD) == sign(ADE) &
              sign(ADC) == sign(ADE)] %>%
    select(gene) %>%
    unique %>%
    unlist %>%
    as.character
```

[ENSG00000005486_RHBDD2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000005486_RHBDD2.pdf) [ENSG00000035403_VCL](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000035403_VCL.pdf) [ENSG00000051108_HERPUD1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000051108_HERPUD1.pdf) [ENSG00000057657_PRDM1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000057657_PRDM1.pdf) [ENSG00000069399_BCL3](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000069399_BCL3.pdf) [ENSG00000071073_MGAT4A](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000071073_MGAT4A.pdf) [ENSG00000077420_APBB1IP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000077420_APBB1IP.pdf) [ENSG00000087502_ERGIC2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000087502_ERGIC2.pdf) [ENSG00000088986_DYNLL1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000088986_DYNLL1.pdf) [ENSG00000089220_PEBP1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000089220_PEBP1.pdf) [ENSG00000104660_LEPROTL1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000104660_LEPROTL1.pdf) [ENSG00000108622_ICAM2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000108622_ICAM2.pdf) [ENSG00000109971_HSPA8](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000109971_HSPA8.pdf) [ENSG00000110852_CLEC2B](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000110852_CLEC2B.pdf) [ENSG00000118197_DDX59](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000118197_DDX59.pdf) [ENSG00000122085_MTERF4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000122085_MTERF4.pdf) [ENSG00000125746_EML2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000125746_EML2.pdf) [ENSG00000125898_FAM110A](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000125898_FAM110A.pdf) [ENSG00000133030_MPRIP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000133030_MPRIP.pdf) [ENSG00000133816_MICAL2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000133816_MICAL2.pdf) [ENSG00000135604_STX11](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000135604_STX11.pdf) [ENSG00000138668_HNRNPD](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000138668_HNRNPD.pdf) [ENSG00000142227_EMP3](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000142227_EMP3.pdf) [ENSG00000163599_CTLA4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000163599_CTLA4.pdf) [ENSG00000169756_LIMS1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000169756_LIMS1.pdf) [ENSG00000173113_TRMT112](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000173113_TRMT112.pdf) [ENSG00000181924_COA4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000181924_COA4.pdf) [ENSG00000203780_FANK1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000203780_FANK1.pdf) [ENSG00000211778_TRAV4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000211778_TRAV4.pdf) [ENSG00000211815_TRAV36DV7](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000211815_TRAV36DV7.pdf) [ENSG00000213626_LBH](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000213626_LBH.pdf) [ENSG00000234745_HLA-B](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000234745_HLA-B.pdf) [ENSG00000242485_MRPL20](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000242485_MRPL20.pdf) 
