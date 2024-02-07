degree.cutoff <- function(G, .cutoff = 3) {
    G.sub <- G
    n.remove <- sum(igraph::degree(G.sub) < .cutoff)
    while(n.remove > 0){
        vv <- igraph::V(G.sub)
        .retain <- vv[igraph::degree(G.sub) >= .cutoff]
        G.sub <- igraph::induced_subgraph(G.sub, .retain)
        n.remove <- sum(igraph::degree(G.sub) < .cutoff)
    }
    return(G.sub)
}

comp.cutoff <- function(G, .cutoff = 10){
    .comp <- igraph::components(G)
    .kk <- which(.comp$csize >= .cutoff) # valid component(s)
    .valid <- which(.comp$membership %in% .kk) # valid nodes
    G <- igraph::induced_subgraph(G, igraph::V(G)[.valid])
    return(G)
}

run.leiden <- function(.knn,
                       .cells,
                       res=.3,
                       symmetrize = T,
                       nrepeat=20,
                       do.correction = T,
                       min.size = 100,
                       deg.cutoff = 3){

    stopifnot(is.list(.knn))

    if(symmetrize){
        .dt <- data.table(from = c(.knn[[1]],.knn[[2]]),
                          to = c(.knn[[2]],.knn[[1]]),
                          weight = c(.knn[[3]], .knn[[3]]))
    } else {
        .dt <- data.table::setDT(.knn[1:3])
        colnames(.dt) <- c("from","to","weight")
    }

    .df <-
        .dt[`from` != `to`,
            .(weight = mean(`weight`)),
            by = .(`from`, `to`)] %>%
        as.data.frame()

    G <- igraph::graph_from_data_frame(.df, directed=FALSE)

    if(do.correction){
        G <- degree.cutoff(G, deg.cutoff)
        G <- comp.cutoff(G, min.size)
    }

    .comp <- igraph::components(G)

    ret <- data.table()
    for(k in unique(.comp$membership)){
        C <- list(quality = 0)
        for(r in 1:nrepeat){
            c.r <- igraph::cluster_leiden(G, resolution_parameter = res, objective_function = "modularity")
            message(paste("quality: ", c.r$quality, "\n"))
            if(r == 1 || max(c.r$quality) > max(C$quality)){
                C <- c.r
            }
        }
        
        .idx <- as.integer(C$names)
        .mem <- C$membership
        .clust <- data.table(tag = .cells[.idx],
                             membership = .mem,
                             component = k)
        ret <- rbind(ret, .clust)
    }

    .sz <- ret[, .(.N), by = .(membership, component)]
    .remove <- .sz[`N` < min.size, ]
    if(nrow(.remove) > 0){
        message("Removing ", nrow(.remove), " clusters")
        ret <- ret %>% anti_join(.remove[, .(membership, component)])
    }

    return(ret)
}

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

bbknn.x <- function(.data, .bbknn, .subset = NULL, .rescale = T) {
    .cols <-
        fread(.data$col, header=F, col.names = "tag") %>%
        parse.tag() %>%
        mutate(j = 1:n()) %>%
        as.data.table()

    .rows <-
        fread(.data$row, col.names="gene", header=F) %>%
        parse.gene()

    if(!is.null(.subset)){
        .rows <- .rows %>%
            filter(hgnc_symbol %in% .subset |
                   ensembl_gene_id %in% .subset |
                   gene %in% .subset)
    }

    vv <- .bbknn$factors.adjusted[.cols$j, , drop = F]
    X <- sweep(.bbknn$U, 2, .bbknn$D, `*`) %*% t(vv)

    ## double scale
    if(.rescale){
        X <- apply(X, 2, scale)
        X.t <- apply(t(X), 2, scale)
        X <- t(X.t)
    }
    ret <- X[.rows$i, , drop = F]
    rownames(ret) <- .rows$hgnc_symbol
    colnames(ret) <- .cols$tag
    return(ret)
}

bbknn.x.melt <- function(.data, .bbknn, .bbknn.umap, .subset = NULL) {

    X <- bbknn.x(.data, .bbknn, .subset)
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
        tag.k <- dict[membership == k]$tag
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
        message("finished the cluster ", k, " vs. others\n")
        ret <- rbind(ret, ret.k)
    }
    return(ret)
}

read.hash <- function(.hash.data,
                      .sample.file = "data/Hashing list MS Treg project.xlsx")
{
    require(readxl)
    
    .hash.mtx <- read.dense(.hash.data$mtx)
    .hash.argmax <- apply(.hash.mtx, 2, which.max)
    
    .hash.cells <- .hash.data$col %>%
        fread(header=F, col.names="tag") %>%
        mutate(hash = .hash.argmax) %>% 
        as.data.table() %>% 
        parse.tag()
    
    .sample.info <-
        readxl::read_xlsx(.sample.file, 1) %>%
        rename(Sample = `Cell type`) %>%
        mutate(hash = gsub("#","",`hash`)) %>%
        mutate(hash = as.integer(`hash`)) %>%
        mutate(disease = substr(`subject`, 1, 2))
    
    .hash.info <-
        left_join(.hash.cells, .sample.info)
}
