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

run.bbknn.umap <- function(.knn,
                           .cells,
                           res=.3,
                           symmetrize = T,
                           nrepeat=20,
                           do.correction = T,
                           ...){

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
        G <- comp.cutoff(G, 10)
        G <- degree.cutoff(G, 3)
    }

    C <- list(quality = 0)
    for(r in 1:nrepeat){
        c.r <- igraph::cluster_leiden(G, resolution_parameter = res, objective_function = "modularity")
        message(paste("quality: ", c.r$quality, "\n"))
        if(r == 1 || max(c.r$quality) > max(C$quality)){
            C <- c.r
        }
    }

    .clust <- data.table(tag = .cells[as.integer(C$names)],
                         membership = C$membership)

    A <- igraph::as_adj(G, attr="weight")

    umap.A <- uwot::optimize_graph_layout(A, verbose = TRUE,
                                          n_sgd_threads = "auto",
                                          ...)

    .dt <-
        data.table(tag = .cells[as.integer(rownames(A))],
                   umap1 = umap.A[,1],
                   umap2 = umap.A[,2]) %>% 
        (function(x) left_join(.clust, x, by = "tag")) %>%
        as.data.table()
    return(.dt)
}
