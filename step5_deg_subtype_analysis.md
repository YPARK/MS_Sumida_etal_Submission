---
title: "Step 5: differential Expression analysis for sub clusters"
date: "2024-02-07"
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
                   .rank = 15, .em.iter = 20, .em.tol = 1e-8, .take.ln = TRUE,
                   knn = 50, impute.by.knn = TRUE, num.threads = 8)

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

### Found 939 unique genes strongly perturbed by MS with FDR 5%

* Up-regulated: 411

* Down-regulated:  687

* Total pairs of genes and clusters: 34,274


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

[PDF](Fig/STEP5//Fig_mTconv_DEG_count.pdf)

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

[CD226_CD226](Fig/STEP5//Fig_mTconv_DEG_example_CD226_CD226.pdf) [CD27_CD27](Fig/STEP5//Fig_mTconv_DEG_example_CD27_CD27.pdf) [CD58_CD58](Fig/STEP5//Fig_mTconv_DEG_example_CD58_CD58.pdf) [ENSG00000000419_DPM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000000419_DPM1.pdf) [ENSG00000009790_TRAF3IP3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000009790_TRAF3IP3.pdf) [ENSG00000009844_VTA1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000009844_VTA1.pdf) [ENSG00000023909_GCLM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000023909_GCLM.pdf) [ENSG00000029639_TFB1M](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000029639_TFB1M.pdf) [ENSG00000034677_RNF19A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000034677_RNF19A.pdf) [ENSG00000043462_LCP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000043462_LCP2.pdf) [ENSG00000057657_PRDM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000057657_PRDM1.pdf) [ENSG00000065978_YBX1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000065978_YBX1.pdf) [ENSG00000070756_PABPC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000070756_PABPC1.pdf) [ENSG00000071073_MGAT4A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000071073_MGAT4A.pdf) [ENSG00000073849_ST6GAL1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000073849_ST6GAL1.pdf) [ENSG00000075239_ACAT1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000075239_ACAT1.pdf) [ENSG00000075426_FOSL2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000075426_FOSL2.pdf) [ENSG00000077420_APBB1IP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000077420_APBB1IP.pdf) [ENSG00000081087_OSTM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000081087_OSTM1.pdf) [ENSG00000081237_PTPRC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000081237_PTPRC.pdf) [ENSG00000088986_DYNLL1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000088986_DYNLL1.pdf) [ENSG00000089157_RPLP0](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000089157_RPLP0.pdf) [ENSG00000090263_MRPS33](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000090263_MRPS33.pdf) [ENSG00000092094_OSGEP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000092094_OSGEP.pdf) [ENSG00000099622_CIRBP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000099622_CIRBP.pdf) [ENSG00000099985_OSM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000099985_OSM.pdf) [ENSG00000100201_DDX17](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100201_DDX17.pdf) [ENSG00000100219_XBP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100219_XBP1.pdf) [ENSG00000100316_RPL3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100316_RPL3.pdf) [ENSG00000100650_SRSF5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100650_SRSF5.pdf) [ENSG00000100664_EIF5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000100664_EIF5.pdf) [ENSG00000101608_MYL12A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000101608_MYL12A.pdf) [ENSG00000101654_RNMT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000101654_RNMT.pdf) [ENSG00000102007_PLP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102007_PLP2.pdf) [ENSG00000102245_CD40LG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102245_CD40LG.pdf) [ENSG00000102409_BEX4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102409_BEX4.pdf) [ENSG00000102760_RGCC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000102760_RGCC.pdf) [ENSG00000103018_CYB5B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000103018_CYB5B.pdf) [ENSG00000103653_CSK](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000103653_CSK.pdf) [ENSG00000104408_EIF3E](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104408_EIF3E.pdf) [ENSG00000104529_EEF1D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104529_EEF1D.pdf) [ENSG00000104904_OAZ1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104904_OAZ1.pdf) [ENSG00000104998_IL27RA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000104998_IL27RA.pdf) [ENSG00000106153_CHCHD2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000106153_CHCHD2.pdf) [ENSG00000106460_TMEM106B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000106460_TMEM106B.pdf) [ENSG00000106560_GIMAP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000106560_GIMAP2.pdf) [ENSG00000108298_RPL19](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000108298_RPL19.pdf) [ENSG00000108622_ICAM2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000108622_ICAM2.pdf) [ENSG00000109046_WSB1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109046_WSB1.pdf) [ENSG00000109062_SLC9A3R1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109062_SLC9A3R1.pdf) [ENSG00000109861_CTSC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109861_CTSC.pdf) [ENSG00000109971_HSPA8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000109971_HSPA8.pdf) [ENSG00000110876_SELPLG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000110876_SELPLG.pdf) [ENSG00000112110_MRPL18](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000112110_MRPL18.pdf) [ENSG00000114209_PDCD10](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000114209_PDCD10.pdf) [ENSG00000114850_SSR3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000114850_SSR3.pdf) [ENSG00000114942_EEF1B2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000114942_EEF1B2.pdf) [ENSG00000115687_PASK](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000115687_PASK.pdf) [ENSG00000115758_ODC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000115758_ODC1.pdf) [ENSG00000116157_GPX7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000116157_GPX7.pdf) [ENSG00000116679_IVNS1ABP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000116679_IVNS1ABP.pdf) [ENSG00000117280_RAB29](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000117280_RAB29.pdf) [ENSG00000117318_ID3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000117318_ID3.pdf) [ENSG00000117410_ATP6V0B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000117410_ATP6V0B.pdf) [ENSG00000118515_SGK1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000118515_SGK1.pdf) [ENSG00000118640_VAMP8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000118640_VAMP8.pdf) [ENSG00000119718_EIF2B2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000119718_EIF2B2.pdf) [ENSG00000120875_DUSP4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000120875_DUSP4.pdf) [ENSG00000120913_PDLIM2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000120913_PDLIM2.pdf) [ENSG00000121067_SPOP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000121067_SPOP.pdf) [ENSG00000121774_KHDRBS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000121774_KHDRBS1.pdf) [ENSG00000121858_TNFSF10](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000121858_TNFSF10.pdf) [ENSG00000122026_RPL21](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122026_RPL21.pdf) [ENSG00000122085_MTERF4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122085_MTERF4.pdf) [ENSG00000122224_LY9](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122224_LY9.pdf) [ENSG00000122406_RPL5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122406_RPL5.pdf) [ENSG00000122862_SRGN](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000122862_SRGN.pdf) [ENSG00000123353_ORMDL2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000123353_ORMDL2.pdf) [ENSG00000123416_TUBA1B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000123416_TUBA1B.pdf) [ENSG00000125089_SH3TC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000125089_SH3TC1.pdf) [ENSG00000125148_MT2A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000125148_MT2A.pdf) [ENSG00000125347_IRF1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000125347_IRF1.pdf) [ENSG00000126524_SBDS](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000126524_SBDS.pdf) [ENSG00000126768_TIMM17B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000126768_TIMM17B.pdf) [ENSG00000127838_PNKD](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000127838_PNKD.pdf) [ENSG00000128340_RAC2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000128340_RAC2.pdf) [ENSG00000129521_EGLN3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000129521_EGLN3.pdf) [ENSG00000131669_NINJ1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000131669_NINJ1.pdf) [ENSG00000132406_TMEM128](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132406_TMEM128.pdf) [ENSG00000132475_H3F3B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132475_H3F3B.pdf) [ENSG00000132591_ERAL1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132591_ERAL1.pdf) [ENSG00000132646_PCNA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132646_PCNA.pdf) [ENSG00000132780_NASP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132780_NASP.pdf) [ENSG00000132965_ALOX5AP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000132965_ALOX5AP.pdf) [ENSG00000133561_GIMAP6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133561_GIMAP6.pdf) [ENSG00000133574_GIMAP4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133574_GIMAP4.pdf) [ENSG00000133639_BTG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133639_BTG1.pdf) [ENSG00000133872_SARAF](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000133872_SARAF.pdf) [ENSG00000134375_TIMM17A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000134375_TIMM17A.pdf) [ENSG00000134954_ETS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000134954_ETS1.pdf) [ENSG00000135046_ANXA1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000135046_ANXA1.pdf) [ENSG00000135316_SYNCRIP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000135316_SYNCRIP.pdf) [ENSG00000136213_CHST12](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000136213_CHST12.pdf) [ENSG00000136490_LIMD2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000136490_LIMD2.pdf) [ENSG00000137078_SIT1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000137078_SIT1.pdf) [ENSG00000137154_RPS6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000137154_RPS6.pdf) [ENSG00000137876_RSL24D1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000137876_RSL24D1.pdf) [ENSG00000138166_DUSP5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000138166_DUSP5.pdf) [ENSG00000138767_CNOT6L](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000138767_CNOT6L.pdf) [ENSG00000139187_KLRG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000139187_KLRG1.pdf) [ENSG00000140988_RPS2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000140988_RPS2.pdf) [ENSG00000142227_EMP3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000142227_EMP3.pdf) [ENSG00000143106_PSMA5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143106_PSMA5.pdf) [ENSG00000143119_CD53](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143119_CD53.pdf) [ENSG00000143333_RGS16](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000143333_RGS16.pdf) [ENSG00000145216_FIP1L1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145216_FIP1L1.pdf) [ENSG00000145741_BTF3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145741_BTF3.pdf) [ENSG00000145779_TNFAIP8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145779_TNFAIP8.pdf) [ENSG00000145860_RNF145](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000145860_RNF145.pdf) [ENSG00000146112_PPP1R18](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000146112_PPP1R18.pdf) [ENSG00000146386_ABRACL](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000146386_ABRACL.pdf) [ENSG00000147403_RPL10](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000147403_RPL10.pdf) [ENSG00000147604_RPL7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000147604_RPL7.pdf) [ENSG00000148154_UGCG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148154_UGCG.pdf) [ENSG00000148303_RPL7A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000148303_RPL7A.pdf) [ENSG00000149273_RPS3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000149273_RPS3.pdf) [ENSG00000149923_PPP4C](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000149923_PPP4C.pdf) [ENSG00000150637_CD226](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000150637_CD226.pdf) [ENSG00000152518_ZFP36L2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000152518_ZFP36L2.pdf) [ENSG00000155660_PDIA4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000155660_PDIA4.pdf) [ENSG00000156508_EEF1A1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000156508_EEF1A1.pdf) [ENSG00000157303_SUSD3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000157303_SUSD3.pdf) [ENSG00000160654_CD3G](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000160654_CD3G.pdf) [ENSG00000161016_RPL8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000161016_RPL8.pdf) [ENSG00000162600_OMA1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162600_OMA1.pdf) [ENSG00000162704_ARPC5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162704_ARPC5.pdf) [ENSG00000162777_DENND2D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000162777_DENND2D.pdf) [ENSG00000163154_TNFAIP8L2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163154_TNFAIP8L2.pdf) [ENSG00000163191_S100A11](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163191_S100A11.pdf) [ENSG00000163468_CCT3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163468_CCT3.pdf) [ENSG00000163519_TRAT1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163519_TRAT1.pdf) [ENSG00000163606_CD200R1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163606_CD200R1.pdf) [ENSG00000163682_RPL9](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163682_RPL9.pdf) [ENSG00000163961_RNF168](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000163961_RNF168.pdf) [ENSG00000164172_MOCS2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000164172_MOCS2.pdf) [ENSG00000164967_RPP25L](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000164967_RPP25L.pdf) [ENSG00000165527_ARF6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000165527_ARF6.pdf) [ENSG00000166012_TAF1D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000166012_TAF1D.pdf) [ENSG00000166925_TSC22D4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000166925_TSC22D4.pdf) [ENSG00000167658_EEF2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000167658_EEF2.pdf) [ENSG00000167851_CD300A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000167851_CD300A.pdf) [ENSG00000167996_FTH1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000167996_FTH1.pdf) [ENSG00000168028_RPSA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168028_RPSA.pdf) [ENSG00000168209_DDIT4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168209_DDIT4.pdf) [ENSG00000168310_IRF2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168310_IRF2.pdf) [ENSG00000168522_FNTA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168522_FNTA.pdf) [ENSG00000168685_IL7R](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000168685_IL7R.pdf) [ENSG00000169100_SLC25A6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000169100_SLC25A6.pdf) [ENSG00000169220_RGS14](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000169220_RGS14.pdf) [ENSG00000170430_MGMT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000170430_MGMT.pdf) [ENSG00000170915_PAQR8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000170915_PAQR8.pdf) [ENSG00000171492_LRRC8D](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171492_LRRC8D.pdf) [ENSG00000171566_PLRG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171566_PLRG1.pdf) [ENSG00000171863_RPS7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000171863_RPS7.pdf) [ENSG00000173113_TRMT112](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000173113_TRMT112.pdf) [ENSG00000173812_EIF1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000173812_EIF1.pdf) [ENSG00000174444_RPL4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174444_RPL4.pdf) [ENSG00000174500_GCSAM](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174500_GCSAM.pdf) [ENSG00000174720_LARP7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174720_LARP7.pdf) [ENSG00000174748_RPL15](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174748_RPL15.pdf) [ENSG00000174953_DHX36](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000174953_DHX36.pdf) [ENSG00000175166_PSMD2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000175166_PSMD2.pdf) [ENSG00000176049_JAKMIP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000176049_JAKMIP2.pdf) [ENSG00000177854_TMEM187](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000177854_TMEM187.pdf) [ENSG00000177917_ARL6IP6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000177917_ARL6IP6.pdf) [ENSG00000179144_GIMAP7](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000179144_GIMAP7.pdf) [ENSG00000181163_NPM1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000181163_NPM1.pdf) [ENSG00000182718_ANXA2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000182718_ANXA2.pdf) [ENSG00000182866_LCK](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000182866_LCK.pdf) [ENSG00000183155_RABIF](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000183155_RABIF.pdf) [ENSG00000183527_PSMG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000183527_PSMG1.pdf) [ENSG00000184009_ACTG1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184009_ACTG1.pdf) [ENSG00000184557_SOCS3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184557_SOCS3.pdf) [ENSG00000184588_PDE4B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000184588_PDE4B.pdf) [ENSG00000185338_SOCS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000185338_SOCS1.pdf) [ENSG00000187091_PLCD1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000187091_PLCD1.pdf) [ENSG00000187514_PTMA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000187514_PTMA.pdf) [ENSG00000188846_RPL14](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000188846_RPL14.pdf) [ENSG00000196230_TUBB](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000196230_TUBB.pdf) [ENSG00000196329_GIMAP5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000196329_GIMAP5.pdf) [ENSG00000196531_NACA](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000196531_NACA.pdf) [ENSG00000197111_PCBP2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197111_PCBP2.pdf) [ENSG00000197329_PELI1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197329_PELI1.pdf) [ENSG00000197471_SPN](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197471_SPN.pdf) [ENSG00000197958_RPL12](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000197958_RPL12.pdf) [ENSG00000198034_RPS4X](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198034_RPS4X.pdf) [ENSG00000198242_RPL23A](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198242_RPL23A.pdf) [ENSG00000198355_PIM3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198355_PIM3.pdf) [ENSG00000198763_MT-ND2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198763_MT-ND2.pdf) [ENSG00000198805_PNP](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198805_PNP.pdf) [ENSG00000198938_MT-CO3](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000198938_MT-CO3.pdf) [ENSG00000204261_PSMB8-AS1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000204261_PSMB8-AS1.pdf) [ENSG00000204519_ZNF551](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000204519_ZNF551.pdf) [ENSG00000204628_RACK1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000204628_RACK1.pdf) [ENSG00000211450_SELENOH](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211450_SELENOH.pdf) [ENSG00000211750_TRBV24-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211750_TRBV24-1.pdf) [ENSG00000211786_TRAV8-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211786_TRAV8-2.pdf) [ENSG00000211788_TRAV13-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211788_TRAV13-1.pdf) [ENSG00000211790_TRAV8-4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211790_TRAV8-4.pdf) [ENSG00000211792_TRAV14DV4](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211792_TRAV14DV4.pdf) [ENSG00000211793_TRAV9-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211793_TRAV9-2.pdf) [ENSG00000211797_TRAV17](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211797_TRAV17.pdf) [ENSG00000211801_TRAV21](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211801_TRAV21.pdf) [ENSG00000211806_TRAV25](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211806_TRAV25.pdf) [ENSG00000211807_TRAV26-1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211807_TRAV26-1.pdf) [ENSG00000211809_TRAV27](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211809_TRAV27.pdf) [ENSG00000211810_TRAV29DV5](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211810_TRAV29DV5.pdf) [ENSG00000211817_TRAV38-2DV8](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000211817_TRAV38-2DV8.pdf) [ENSG00000213145_CRIP1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213145_CRIP1.pdf) [ENSG00000213626_LBH](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213626_LBH.pdf) [ENSG00000213658_LAT](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213658_LAT.pdf) [ENSG00000213719_CLIC1](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000213719_CLIC1.pdf) [ENSG00000232956_SNHG15](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000232956_SNHG15.pdf) [ENSG00000235532_LINC00402](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000235532_LINC00402.pdf) [ENSG00000241657_TRBV11-2](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000241657_TRBV11-2.pdf) [ENSG00000243646_IL10RB](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000243646_IL10RB.pdf) [ENSG00000243749_TMEM35B](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000243749_TMEM35B.pdf) [ENSG00000243927_MRPS6](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000243927_MRPS6.pdf) [ENSG00000269335_IKBKG](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000269335_IKBKG.pdf) [ENSG00000277734_TRAC](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000277734_TRAC.pdf) [ENSG00000278030_TRBV7-9](Fig/STEP5//Fig_mTconv_DEG_example_ENSG00000278030_TRBV7-9.pdf) [TIGIT_TIGIT](Fig/STEP5//Fig_mTconv_DEG_example_TIGIT_TIGIT.pdf) 

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
                   .rank = 15, .em.iter = 20, .em.tol = 1e-8, .take.ln = TRUE,
                   knn = 50, impute.by.knn = TRUE, num.threads = 8)

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

### Found 1,011 unique genes strongly perturbed by MS with FDR 5%

* Up-regulated: 263

* Down-regulated:  789

* Total pairs of genes and clusters: 50,116

![](Fig/STEP5/Fig_mTreg_DEG_count-1.png)<!-- -->

[PDF](Fig/STEP5//Fig_mtreg_DEG_count.pdf)


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

[ENSG00000035403_VCL](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000035403_VCL.pdf) [ENSG00000051108_HERPUD1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000051108_HERPUD1.pdf) [ENSG00000057657_PRDM1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000057657_PRDM1.pdf) [ENSG00000071073_MGAT4A](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000071073_MGAT4A.pdf) [ENSG00000077420_APBB1IP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000077420_APBB1IP.pdf) [ENSG00000088986_DYNLL1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000088986_DYNLL1.pdf) [ENSG00000089220_PEBP1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000089220_PEBP1.pdf) [ENSG00000092621_PHGDH](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000092621_PHGDH.pdf) [ENSG00000108590_MED31](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000108590_MED31.pdf) [ENSG00000108622_ICAM2](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000108622_ICAM2.pdf) [ENSG00000108679_LGALS3BP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000108679_LGALS3BP.pdf) [ENSG00000109971_HSPA8](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000109971_HSPA8.pdf) [ENSG00000110852_CLEC2B](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000110852_CLEC2B.pdf) [ENSG00000117592_PRDX6](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000117592_PRDX6.pdf) [ENSG00000121895_TMEM156](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000121895_TMEM156.pdf) [ENSG00000125898_FAM110A](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000125898_FAM110A.pdf) [ENSG00000133030_MPRIP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000133030_MPRIP.pdf) [ENSG00000135486_HNRNPA1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000135486_HNRNPA1.pdf) [ENSG00000135604_STX11](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000135604_STX11.pdf) [ENSG00000142227_EMP3](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000142227_EMP3.pdf) [ENSG00000163599_CTLA4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000163599_CTLA4.pdf) [ENSG00000166848_TERF2IP](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000166848_TERF2IP.pdf) [ENSG00000168569_TMEM223](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000168569_TMEM223.pdf) [ENSG00000169100_SLC25A6](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000169100_SLC25A6.pdf) [ENSG00000169756_LIMS1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000169756_LIMS1.pdf) [ENSG00000173113_TRMT112](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000173113_TRMT112.pdf) [ENSG00000178977_LINC00324](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000178977_LINC00324.pdf) [ENSG00000187257_RSBN1L](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000187257_RSBN1L.pdf) [ENSG00000203780_FANK1](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000203780_FANK1.pdf) [ENSG00000211778_TRAV4](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000211778_TRAV4.pdf) [ENSG00000213626_LBH](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000213626_LBH.pdf) [ENSG00000234745_HLA-B](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000234745_HLA-B.pdf) [ENSG00000242485_MRPL20](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000242485_MRPL20.pdf) [ENSG00000243749_TMEM35B](Fig/STEP5//Fig_mTreg_DEG_example_ENSG00000243749_TMEM35B.pdf) 
