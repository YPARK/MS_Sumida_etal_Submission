options(stringsAsFactors = FALSE)

library(tidyverse)
library(data.table)
library(ggrepel)

`%&%` <- function(a, b) paste0(a, b)

sigmoid <- function(x) 1/(1 + exp(-x))

logit <- function(x) log(x) - log(1 - x)

meta.z <- function(z) mean(z) * sqrt(length(z))

num.sci <- function(x) format(x, digits = 2, scientific = TRUE)

num.round <- function(x, d=2) round(x, digits = d) %>% as.character()

num.int <- function(x) format(x, big.mark = ',')

.gsub.rm <- function(x, pat) gsub(x, pattern = pat, replacement = '')

.unlist <- function(...) unlist(..., use.names = FALSE)

.select <- function(...) dplyr::select(...) %>% .unlist()

.mkdir <- function(...) dir.create(..., recursive=TRUE, showWarnings = FALSE)

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

.gg.save <- function(filename, cat.link=TRUE, ...) {
    if(file.exists(filename)) {
        log.msg('File already exits: %s', filename)
    } else {
        ggplot2::ggsave(..., filename = filename,
                        limitsize = FALSE,
                        units = 'in',
                        dpi = 300,
                        useDingbats = FALSE)
    }
    if(cat.link){
        cat("\n\n[PDF](" %&% filename %&% ")\n\n")
    }
}

.gg.plot <- function(...) {
    ggplot2::ggplot(...) +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.background = element_blank(),
                       plot.margin = unit(c(0,.5,0,.5), 'lines'),
                       panel.background = element_rect(size = 0, fill = 'gray95'),
                       strip.background = element_blank(),
                       legend.background = element_blank(),
                       legend.text = element_text(size = 6),
                       legend.title = element_text(size = 6),
                       axis.title = element_text(size = 8),
                       legend.key.width = unit(1, 'lines'),
                       legend.key.height = unit(.2, 'lines'),
                       legend.key.size = unit(1, 'lines'),
                       axis.line = element_line(color = 'gray20', size = .5),
                       axis.text = element_text(size = 6))
}

load.data <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

col.order <- function(pair.tab, .ro, ret.tab = FALSE) {

    M = pair.tab %>%
        dplyr::select(row, col, weight) %>%
        mutate(row = factor(row, .ro)) %>%
        tidyr::spread(key = col, value = weight, fill = 0)

    co = order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co = colnames(M)[-1][co]
    if(ret.tab) {
        ret = pair.tab %>%
            mutate(row = factor(row, .ro)) %>% 
            mutate(col = factor(col, .co))
    } else {
        ret = .co
    }
    return(ret)
}

order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)
    
    .tab = pair.tab %>% dplyr::select(row, col, weight)

    M = .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr = M[, 1] %>% unlist(use.names = FALSE)
    cc = colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro = row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## log.msg('Sort the rows: %d', length(ro))
    co = order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret = pair.tab %>%
            mutate(row = factor(row, rr[ro])) %>%
            mutate(col = factor(col, cc[co]))
    } else {
        ret = list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}

################################################################
#' Read gene set information
read.geneset <- function(ensg) {

    ANNOT.FILE = "../data/msigdb.txt.gz"

    annot.tab = fread(ANNOT.FILE)[, transcript_start := as.integer(transcript_start)]
    annot.tab = annot.tab[, transcript_end := as.integer(transcript_end)]
    annot.tab = annot.tab[ensembl_gene_id %in% ensg]

    gs.tab = annot.tab[, .(ensembl_gene_id, hgnc_symbol, gs_subcat, gs_name)][, val := 1]

    C.mat = dcast(gs.tab, ensembl_gene_id + hgnc_symbol ~ gs_name,
                  fun.aggregate = length, value.var = "val")

    list(C = C.mat, gs = gs.tab)
}

################################################################
#' Read gene ontology and coding genes
read.ontology <- function() {

    .temp.file <- ".ensembl.ontology.rdata"

    if(file.exists(.temp.file)) {
        load(.temp.file)

    } else {

        ensembl = biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                                   host='uswest.ensembl.org',
                                   path='/biomart/martservice',
                                   dataset='hsapiens_gene_ensembl')

        ensembl.hs = biomaRt::useDataset('hsapiens_gene_ensembl',mart=ensembl)

        .attr = c('ensembl_gene_id',
                  'hgnc_symbol',
                  'chromosome_name',
                  'transcription_start_site',
                  'transcript_start',
                  'transcript_end',
                  'description',
                  'percentage_gene_gc_content')

        .temp = biomaRt::getBM(attributes=.attr,
                               filters=c('biotype'),
                               values=c('protein_coding'),
                               mart=ensembl.hs,
                               useCache = FALSE)

        genes.desc = .temp %>%
            dplyr::select(hgnc_symbol, description) %>%
            unique()

        .temp <- as.data.table(.temp)

        coding.genes = .temp[,
                             .(hgnc_symbol = paste(unique(hgnc_symbol), collapse='|'),
                               transcription_start_site = mean(transcription_start_site),
                               transcript_start = min(transcript_start),
                               transcript_end = max(transcript_end)),
                             by = .(chromosome_name, ensembl_gene_id)]

        ensg.tot = unique(coding.genes$ensembl_gene_id)
        go.attr = c('ensembl_gene_id', 'go_id',
                    'name_1006', 'namespace_1003')

        genes.go = biomaRt::getBM(filters = 'ensembl_gene_id',
                                  values = ensg.tot,
                                  attributes = go.attr,
                                  mart = ensembl.hs,
                                  useCache = FALSE)
        .ensembl.ontology <- 
            list(coding = coding.genes,
                 go = genes.go,
                 desc = genes.desc)

        save(.ensembl.ontology, file=.temp.file)
    }

    return(.ensembl.ontology)
}

################################################################
#' Run enrichment test and report hypergeometric p-value
summary.enrichment.stat <- function(tab) {

    m = max(tab$annot.size)                 ## white balls
    n = max(tab$ntot) - m                   ## black balls
    k = max(tab$n.drawn)                    ## balls drawn
    q = length(unique(tab$ensembl_gene_id)) ## whilte balls drawn
    p = k * m / (m + n)                     ## expected by chance

    ## calculate hypergeometric p-value
    ## phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
    ##    q: vector of quantiles representing the number of white balls
    ##       (drawn without replacement from an urn which contains both
    ##       black and white balls).
    ##    m: the number of white balls in the urn.
    ##    n: the number of black balls in the urn.
    ##    k: the number of balls drawn from the urn.
    ##    p: probability, it must be between 0 and 1.

    pval = phyper(q, m, n, k, lower.tail = FALSE)

    data.table(n.white = m,
               n.drawn = k,
               n.overlap = q,
               n.expected = p,
               pval = pval)
}

