#!/usr/bin/env Rscript

library(tidyverse)
library(R.utils)
library(data.table)
library(Matrix)

batches = list.files("CITE_seq/",
                     pattern="batch",
                     full.names=TRUE)

.files = batches %>%
    sapply(FUN=list.files, pattern="*RData", full.names=TRUE) %>%
    unlist(use.names=FALSE) %>%
    Filter(f = function(x) !is.na(str_locate(x, "raw_")[1])) %>% 
    Filter(f = function(x) !is.na(str_locate(x, "rna_")[1]))

load.data <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

squeeze.mtx <- function(.mtx, .feat) {

    .dt = as.data.table(as.matrix(.mtx))
    .dt = .dt[, .feat := unlist(.feat)] %>% as.data.table()
    .dt = .dt[, lapply(.SD, sum), by=.(.feat)]

    .feat = .dt$.feat
    .dt[, .feat := NULL]
    .dt = as.matrix(.dt)

    rownames(.dt) = .feat
    Matrix(.dt, sparse=TRUE)
}

.convert <- function(ff) {

    ff.prot = str_replace(ff, "rna", "protein")

    b = basename(ff) %>%
        str_split(pattern="_") %>%
        (function(x) x[[1]][4])

    d = basename(dirname(ff))

    mtx.out.file = paste0("batch_", b, ".mtx.gz")

    if(file.exists(mtx.out.file)) return(NULL)

    ## read protein mtx
    mtx.prot = load.data(ff.prot)

    prot.feature =
        tibble(f = rownames(mtx.prot)) %>%
        dplyr::mutate(f = str_remove_all(f, "[-]+Total[A-Z]+")) %>% 
        dplyr::mutate(f = str_remove(f, "anti[-]*human[-]*")) %>%
        tidyr::separate("f", c("f","etc"), sep="[-(]") %>%
        dplyr::mutate(f = paste0("anti_", str_to_upper(f))) %>% 
        dplyr::select(f)

    mtx.prot = squeeze.mtx(mtx.prot, prot.feature)

    cat("Read the proteins ...\n", file=stderr())

    barcode = colnames(mtx.prot)
    
    ## read rna mtx
    mtx.rna = load.data(ff)
    rna.feature = rownames(mtx.rna)
    mtx.rna = squeeze.mtx(mtx.rna, rna.feature)

    stopifnot(all(colnames(mtx.rna) == barcode))

    cat("Read the genes ...\n", file=stderr())

    mtx.joint = rbind(mtx.prot, mtx.rna) 

    feature = rownames(mtx.joint) %>% as_tibble
    barcode = colnames(mtx.joint) %>% as_tibble

    ## write down the antibody features only
    Matrix::writeMM(mtx.joint, str_remove(mtx.out.file, ".gz"))
    gzip(str_remove(mtx.out.file, ".gz"))

    write_tsv(barcode, paste0("batch_", b, ".barcode.tsv.gz"), col_names=FALSE)
    write_tsv(feature, paste0("batch_", b, ".feature.tsv.gz"), col_names=FALSE)

    return(NULL)
}

. = sapply(.files, .convert)

## combine features
if(!file.exists("features.tsv.gz")) {
    list.files(".", "feature.tsv.gz") %>%
        Filter(f=file.exists) %>%
        lapply(read_tsv, col_names="gene", col_types="c") %>%
        bind_rows() %>%
        distinct() %>%
        arrange(gene) %>% 
        write_tsv("features.tsv.gz", col_names=FALSE)
}

################################################################
## take phenotype labels
.read.hash <- function(ff) {

    b = str_split(basename(ff), pattern="_")[[1]][1] %>%
        str_remove("batch")

    read_delim(ff, delim=",") %>%
        select(X1, hash, label, pred_mcplda) %>%
        rename(barcode = X1) %>%
        mutate(batch = b)
}

.read.hash.files <- function(pat = "pc10_sd03.csv.gz") {
    hash.tib = tibble()

    .files = list.files("CITE_seq/classification_results/",
                        pat,
                        full.names = TRUE)

    hash.tib = bind_rows(lapply(.files, .read.hash)) %>%
        bind_rows(hash.tib) %>%
        mutate(hash=str_remove_all(hash, "-")) %>% 
        mutate(hash=str_remove_all(hash, "anti")) %>% 
        mutate(hash=str_remove_all(hash, "TotalC")) %>% 
        mutate(hash=str_remove_all(hash, "human")) %>%
        mutate(hash=str_remove_all(hash, "Hashtag")) %>% 
        mutate(hash=as.integer(hash)) %>%
        mutate(disease=if_else(hash <= 2, "HC", "MS"))
}

if(!file.exists("hashtag.tsv.gz")) {
    hash.tib = .read.hash.files("pc10_sd03.csv.gz")
    write_tsv(hash.tib, "hashtag.tsv.gz")
}
