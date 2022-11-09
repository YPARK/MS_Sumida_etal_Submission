library(torch)
library(mmutilR)

#' @param xx torch matrix
#' @return R's native matrix format
as.mat <- function(xx){
    as.matrix(xx$to(dev=torch_device("cpu")))
}

#' @param .mean mean
#' @param .lnvar log(variance)
normal.stoch <- function(.mean, .lnvar) {
    .eps <- torch_randn_like(.lnvar)
    .sig <- torch_exp(.lnvar / 2.)
    .mean + .eps * .sig
}

#' @param .mean mean vector
#' @param .lnvar log(variance)
kl.loss <- function(.mean, .lnvar) {
    -0.5 * torch_sum(1. + .lnvar - torch_pow(.mean, 2.) - torch_exp(.lnvar), dim = -1);
}

#' FC ReLU layer generator
#' @param in_ input dimension
#' @param out_ output dimension
#' @return a generator
#'
build.fc.relu <- function(in_, out_) {
    nn_sequential(nn_linear(in_, out_),
                  nn_relu())
}

#' Build a stack of layers
#' @param in_ input dimension (first layer)
#' @param layers dimensions for the subsequent layers
#' @param generator a shared generator function
#' @param stdizer standardizer
#' @param name a header for the names
#'
build.stack <- function(in_, layers, generator = nn_linear, name = "", stdizer = nn_batch_norm1d) {
    ret <-
        nn_module(
            classname = "stack.layers",
            initialize = function(in_,
                                  .layers,
                                  .generator,
                                  .name) {
                d.prev <- in_
                for(l in 1:length(.layers)) {
                    d.this <- .layers[l]
                    .name.l <- paste0(name, ".a", (l - 1))
                    self$add_module(.name.l,
                                    module = .generator(d.prev, d.this))
                    if(!is.null(stdizer)){
                        .name.l <- paste0(name, ".d", (l - 1))
                        self$add_module(.name.l,
                                        module = stdizer(d.this))
                    }
                    d.prev <- d.this
                }
            },
            forward = function(input) {
                for(module in private$modules_) {
                    input <- module(input)
                }
                input
            })
    return(ret(in_, layers, generator, name))
}

###################
## cite-seq data ##
###################

CITEseqData <- dataset(
    name = "single-cell CITE-Seq data",
    initialize = function(sc.data, MY.DEV = torch_device("cpu")) {
        self$sc.data <- sc.data
        self$DEV <- MY.DEV
        self$info <- mmutilR::rcpp_mmutil_info(self$sc.data$mtx)
        ncells <- self$info$max.col
        .vec <- self$read.vec(sc.data$row)
        self$protein.pos <- which(stringr::str_starts(.vec, "anti_"))
        self$gene.pos <- which(!stringr::str_starts(.vec, "anti_"))
        self$proteins <- .vec[self$protein.pos]
        self$genes <- .vec[self$gene.pos]
        self$mem.idx <- mmutilR::rcpp_mmutil_read_index(sc.data$idx)
        self$read.full()
        message("Initialized scdata")
    },
    read.full = function() {
        .mtx <- self$sc.data$mtx
        .idx <- self$mem.idx
        .loc <- 1:self$info$max.col
        .out <- mmutilR::rcpp_mmutil_read_columns_sparse(.mtx, .idx, .loc)
        .ind <- rbind(.out$col, .out$row) # column is row & vice versa
        .sz <- c(.out$max.col, .out$max.row)
        xx <- torch_sparse_coo_tensor(.ind, .out$val, .sz, dtype = torch_float16())
        xx <- torch_clamp(xx$to_dense(), min = 0, max = 1e4)
        message("Porting to the device")
        self$x.list <- lapply(1:nrow(xx), function(r) {
            cat(r,"\r",file=stderr());flush(stderr())
            xx[r, , drop=FALSE]$to_sparse()$to(dtype = torch_float16(),
                                               device=self$DEV)
        })
        message("Distributed to the device")
    },
    read.vec = function(x) {
        if(stringr::str_ends(x, ".gz$")) {
            con <- gzfile(x)
            ret <- readLines(con)
            close(con)
            return(ret)
        }
        readLines(x)
    },
    .getitem = function(.loc) {
        xx <- torch_cat(self$x.list[.loc], dim=1)$to_dense()$to(device=self$DEV)
        .gene <- xx[, self$gene.pos, drop = FALSE]$to(dtype=torch_float32())
        .prot <- xx[, self$protein.pos, drop = FALSE]$to(dtype=torch_float32())
        list(gene = .gene, protein = .prot)
    },
    .length = function() {
        self$info$max.col
    })

###################
## encoder model ##
###################
build.ETM.encoder <-
    nn_module(
        classname = "bimodal_citseq_encoder",
        initialize = function(d.gene, d.prot, K = 16, encoder.layers = c(16, 16)) {
            ## (0) data normalizer
            self$bn.prot <- nn_batch_norm1d(d.prot)
            self$bn.gene <- nn_batch_norm1d(d.gene)
            ## (1) data -> fc
            self$fc.gene <- build.stack(d.gene, encoder.layers,
                                        generator = build.fc.relu,
                                        name = "ETM.enc.fc.gene",
                                        stdizer = nn_batch_norm1d)
            self$fc.prot <- build.stack(d.prot, encoder.layers,
                                        generator = build.fc.relu,
                                        name = "ETM.enc.fc.prot",
                                        stdizer = nn_batch_norm1d)
            d <- encoder.layers[length(encoder.layers)]
            ## (2) fc -> K
            self$z.prot.mean <- nn_linear(d, K)
            self$z.prot.lnvar <- nn_linear(d, K)
            self$z.gene.mean <- nn_linear(d, K)
            self$z.gene.lnvar <- nn_linear(d, K)
        },
        forward = function(.batch) {
            ## normalization of proteins
            yy <- torch_log1p(.batch$protein) %>%
                nnf_normalize(dim=2, eps=1e-8)
            yy <- self$bn.prot(yy)
            ## normalization of genes
            xx <- torch_log1p(.batch$gene) %>%
                nnf_normalize(dim=2, eps=1e-8)
            xx <- self$bn.gene(xx)
            ## (1) data -> fc
            ss.gene <- self$fc.gene(xx)
            ss.prot <- self$fc.prot(yy)
            ## (2) stochastic
            mm.gene <- self$z.gene.mean(ss.gene)
            lv.gene <- torch_clamp(self$z.gene.lnvar(ss.gene), -4.0, 4.0)
            z.gene <- normal.stoch(mm.gene, lv.gene)
            mm.prot <- self$z.prot.mean(ss.prot)
            lv.prot <- torch_clamp(self$z.prot.lnvar(ss.prot), -4.0, 4.0)
            z.prot <- normal.stoch(mm.prot, lv.prot)
            ## (3) combine them equally
            z <- z.gene * .5 + z.prot * .5
            list(z = z,
                 z.mean = (mm.gene + mm.prot) * .5,
                 z.gene.mean = mm.gene,
                 z.gene.lnvar = lv.gene,
                 z.prot.mean = mm.prot,
                 z.prot.lnvar = lv.prot)
        })

###################
## decoder model ##
###################

build.ETM.decoder <-
    nn_module(
        classname = "bimodal_citeseq_decoder",
        initialize = function(d.gene, d.prot, K, pip0 = 0.1, v0 = 1){
            self$prot.bias <- nn_parameter(torch_randn(1, d.prot))
            self$gene.bias <- nn_parameter(torch_randn(1, d.gene))
            ## model parameters
            self$prot.mean <- nn_parameter(torch_randn(K, d.prot))
            self$prot.lnvar <- nn_parameter(torch_randn(K, d.prot))
            self$prot.logit <- nn_parameter(torch_randn(K, d.prot))
            self$gene.mean <- nn_parameter(torch_randn(K, d.gene))
            self$gene.lnvar <- nn_parameter(torch_randn(K, d.gene))
            self$gene.logit <- nn_parameter(torch_randn(K, d.gene))
            ## helper function
            self$log.softmax <- nn_log_softmax(dim=2)
            ## hyper-parameters
            self$logit.0 <- nn_parameter(torch_logit(pip0, eps=1e-6), requires_grad = FALSE)
            self$lnvar.0 <- nn_parameter(torch_log(v0), requires_grad = FALSE)
        },
        get.mean = function(.logit, .mean) {
            .pip <- torch_sigmoid(.logit)
            .mean * .pip
        },
        get.var = function(.logit, .mean, .lnvar) {
            .pip <- torch_sigmoid(.logit)
            .var <- .pip * (1 - .pip) * torch_square(.mean)
            .var + .pip * torch_exp(.lnvar)
        },
        get.stoch = function(.logit, .mean, .lnvar) {
            .pip <- torch_sigmoid(.logit)
            .mean <- .mean * .pip
            .var <- .pip * (1 - .pip) * torch_square(.mean)
            .var <- .var + .pip * torch_exp(.lnvar)
            .eps <- torch_randn_like(.var)
            .mean + .eps * torch_sqrt(.var)
        },
        get.gene.stoch = function() {
            self$get.stoch(self$gene.logit, self$gene.mean, self$gene.lnvar)
        },
        get.prot.stoch = function() {
            self$get.stoch(self$prot.logit, self$prot.mean, self$prot.lnvar)
        },
        forward = function(zz) {
            hh <- self$soft.max(zz)
            beta.gene <- self$safe.exp(self$get.gene.stoch() + self$gene.bias)
            beta.prot <- self$safe.exp(self$get.prot.stoch() + self$prot.bias)
            .prot <- torch_mm(hh, beta.prot)
            .gene <- torch_mm(hh, beta.gene)
            kl.prot <- self$sparse.kl.loss(self$prot.logit,
                                           self$prot.mean,
                                           self$prot.lnvar)
            kl.gene <- self$sparse.kl.loss(self$gene.logit,
                                           self$gene.mean,
                                           self$gene.lnvar)
            list(protein = .prot, gene = .gene, kl = kl.gene + kl.prot)
        },
        #' D(q(beta) || p(beta))
        #' @param spike.logit logit of the spike parameters
        #' @param slab.mean mean of the slab parameters
        #' @param slab.lnvar log-var of the slab parameters
        sparse.kl.loss = function(spike.logit,
                                  slab.mean,
                                  slab.lnvar) {
            logit.0 <- self$logit.0
            lnvar.0 <- self$lnvar.0
            ## PIP KL between p and p0
            ## p * ln(p / p0) + (1-p) * ln(1-p/1-p0)
            ## = p * ln(p / 1-p) + ln(1-p) +
            ##   p * ln(1-p0 / p0) - ln(1-p0)
            ## = sigmoid(logit) * logit - softplus(logit)
            ##   - sigmoid(logit) * logit0 + softplus(logit0)
            pip.hat <- torch_sigmoid(spike.logit)
            kl.pip.1 <- pip.hat * (spike.logit - logit.0)
            kl.pip <- kl.pip.1 - nnf_softplus(spike.logit) + nnf_softplus(logit.0)
            ## Gaussian KL between N(μ,ν) and N(0, v0)
            v0.inv <- torch_exp(-lnvar.0)
            sq.term <- v0.inv * (torch_square(slab.mean) + torch_exp(slab.lnvar))
            kl.g <- -0.5 * (1. + slab.lnvar - lnvar.0 - sq.term)
            ## Combine both logit and Gaussian KL
            torch_sum(kl.pip + pip.hat * kl.g)
        },
        safe.exp = function(xx) {
            torch_exp(torch_clamp(xx, -10, 10))
        },
        soft.max = function(xx) {
            torch_exp(self$log.softmax(xx))
        })


#####################################
## combine ETM encoder and decoder ##
#####################################

build.ETM <-
    nn_module(
        classname = "bimodal_citeseq_ETM",
        initialize = function(d.gene, d.prot, K, encoder.layers) {
            self$K <- K
            self$encoder <- build.ETM.encoder(d.gene, d.prot, K, encoder.layers)
            self$decoder <- build.ETM.decoder(d.gene, d.prot, K)
            self$log.softmax <- nn_log_softmax(dim=2)
        },
        forward = function(.batch.data) {
            .softmax <- nn_log_softmax(2)
            .enc <- self$encoder(.batch.data)
            .recon <- self$decoder(.enc$z)
            h <- self$soft.max(.enc$z.mean)
            llik.gene <- self$dir.llik(.batch.data$gene, .recon$gene)
            llik.prot <- self$dir.llik(.batch.data$protein, .recon$protein)
            kl.gene <- self$kl.loss(.enc$z.gene.mean, .enc$z.gene.lnvar)
            kl.prot <- self$kl.loss(.enc$z.prot.mean, .enc$z.prot.lnvar)
            list(h = h,
                 recon = .recon,
                 llik.gene = llik.gene,
                 llik.protein = llik.prot,
                 z.kl = kl.gene + kl.prot,
                 z.gene.mean = .enc$z.gene.mean,
                 z.prot.mean = .enc$z.prot.mean,
                 model.kl = .recon$kl)
        },
        soft.max = function(xx) {
            torch_exp(self$log.softmax(xx))
        },
        #' Dirichlet log-likelihood:
        #' lgamma(sum a) - lgamma(sum a + x)
        #' sum lgamma(a + x) - lgamma(a)
        #' @param xx
        #' @param aa
        #' @return log-likelihood
        dir.llik = function(xx, aa, a0 = 1e-2, .tol = 1e-8){
            term1 <- (torch_lgamma(torch_sum(aa + a0, dim=2)) -
                      torch_lgamma(torch_sum(aa + xx + a0, dim=2)))
            term2 <- torch_sum(torch_where(xx > .tol,
                                           torch_lgamma(aa + xx + a0) -
                                           torch_lgamma(aa + a0),
                                           torch_zeros_like(xx)),
                               dim=2)
            term1 + term2
        },
        #' @param .mean mean vector
        #' @param .lnvar log(variance)
        kl.loss = function(.mean, .lnvar) {
            -0.5 * torch_sum(1. + .lnvar - torch_square(.mean) - torch_exp(.lnvar), dim=2)
        })

#' Train a model and save the results
#'
#' @param .data
#' @param kk
#' @param encoder.layers
#' @param model.file
#' @param stat.file
train.model <- function(.data, kk,
                        encoder.layers,
                        model.file,
                        stat.file,
                        max.epoch = 3000,
                        .seed = 377,
                        .lr = 1e-2,
                        batch.size = 128){

    if(file.exists(model.file)) { return(NULL) }
    if(file.exists(stat.file)) { return(NULL) }

    torch_manual_seed(.seed)

    d.prot <- length(.data$proteins)
    d.gene <- length(.data$genes)

    model <- build.ETM(d.gene, d.prot, kk, encoder.layers)$to(device=.data$DEV)
    opt <- optim_adam(model$parameters, lr = .lr)

    model$proteins <- .data$proteins
    model$genes <- .data$genes

    llik.g <- llik.p <- c()

    ncells <- .data$.length()
    n.batch <- ceiling(ncells/batch.size)
    nn <- n.batch * batch.size
    model.kl.scale <- 1/nn

    CPU <- torch_device("cpu")

    for(tt in seq(0, max.epoch)){

        .samples <- sample(ncells, nn, replace=TRUE)

        .llik.p <- 0
        .llik.g <- 0

        for(b in seq(0, n.batch-1)){
            .batch <- .samples[seq(b*batch.size+1, (b+1)*batch.size)]
            .batch.data <- .data$.getitem(.batch)
            opt$zero_grad() # calculate gradient fresh
            model$train()   # set the training mode
            out <- model(.batch.data) # take the forward pass
            llik <- out$llik.gene + out$llik.protein
            .model.kl <- out$model.kl * model.kl.scale
            loss <- torch_mean(out$z.kl - llik) + .model.kl
            loss$backward()
            opt$step()
            .llik.g <- .llik.g + out$llik.gene$sum()$item()
            .llik.p <- .llik.p + out$llik.protein$sum()$item()
        }
        llik.g <- c(llik.g, as.numeric(.llik.g)/ncells)
        llik.p <- c(llik.p, as.numeric(.llik.p)/ncells)

        hh <- as.matrix(out$h$to(device=CPU))
        h.hist <- table(apply(hh, 1, which.max))
        message(paste(str_c("", names(h.hist), ":", h.hist), collapse="|"))

        message("epoch: ", tt,
                ", gene: ", tail(llik.g, 1),
                ", protein: ", tail(llik.p, 1))
    }

    ## final latent encoding
    encoder <- model$encoder
    encoder$eval()

    .make.latent.mat <- function() {
        matrix(NA,
               nrow = .data$.length(),
               ncol = model$K)
    }

    latent.mean <- .make.latent.mat()
    latent.gene.mean <- .make.latent.mat()
    latent.prot.mean <- .make.latent.mat()

    .softmax <- model$decoder$soft.max

    for(lb in seq(0, ncells, batch.size)) {
        ub <- min(lb + batch.size, ncells)
        r <- seq(lb+1, ub)
        dr <- .data$.getitem(r)
        .out <- encoder(dr)
        latent.mean[r, ] <- as.matrix(.softmax(.out$z.mean)$to(device=CPU))
        latent.gene.mean[r, ] <- as.matrix(.softmax(.out$z.gene.mean)$to(device=CPU))
        latent.prot.mean[r, ] <- as.matrix(.softmax(.out$z.prot.mean)$to(device=CPU))
    }

    torch_save(model, path=model.file)
    saveRDS(list(llik = list(gene = llik.g,
                             protein = llik.p),
                 latent.mean = latent.mean,
                 latent.gene.mean = latent.gene.mean,
                 latent.prot.mean = latent.prot.mean),
            stat.file)
}
