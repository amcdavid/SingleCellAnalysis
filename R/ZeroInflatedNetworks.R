##' Fit Meinhousen-Buhlmann to a SingleCellAssay object 
##'
##' Regresses the dichotomous and continuous components of each gene in \code{sc}
##' on every other gene.
##' \code{additive.effects} from \code{cData(sc)} are included unpenalized.
##' Currently only \code{gene.predictors} as 'zero.inflated' is supported.
##' ... is passed along to cv.glmnet, see documentation there.
##' @param sc SingleCellAssay object on a thresholded layer
##' @param additive.effects character vector, possibly using formula syntax of columns from \code{cData(sc)} to be included as unpenalized terms.
##' @param min.freq genes below this frequency are excluded as predictors and dependent variables
##' @param gene.predictors one of 'zero.inflated' or 'hurdle'.  See details.
##' @param precenter How should centering/scaling be done with respect to continuous regressions.  TRUE if centering should be done with respect to all cells; FALSE if centering should be done only with respect to expressed cells
##' When precenter=TRUE, cv.glmnet will not standardize.
##' @param precenter.fun a function called to center the expression matrix prior to calling glmnet
##' @param response a character vector, one of 'zero.inflated', 'hurdle', or 'cg.regression'
##' @param ... passed to cv.glmnet
##' @return 2-D list of cv.glmnet objects with attributes
##' @importFrom glmnet glmnet cv.glmnet
##' @export
fitZifNetwork <- function(sc, additive.effects, min.freq=.05, gene.predictors='zero.inflated', precenter=TRUE, precenter.fun=scale, response='hurdle', ...){
    gene.predictors <- match.arg(gene.predictors, c('zero.inflated', 'hurdle'))
    response <- match.arg(response, c('hurdle', 'zero.inflated'))
    sub <- sc[, freq(sc)>min.freq]
    genes <- fData(sub)$primerid

    ## Additive (un-penalized) variables named from cData
    additive.mat <- model.matrix(as.formula(sprintf('~ %s', paste(additive.effects, collapse='+'))), cData(sub))[, -1, drop=FALSE] #no intercept
    additive.dim <- ncol(additive.mat)
        

    ## Untested code to decompose predictors into continuous/dichtomous
    if(gene.predictors == 'hurdle'){
        genes.appear <- 2
        ## Additive effects, then centered continuous, then dichotomous
        model.mat <- as.matrix(cbind(additive.mat,
                                     as.data.frame(xform(exprs(sub))), as.data.frame(exprs(sub)>0)))
    } else if(gene.predictors == 'zero.inflated'){
        ## Additive effects, then zero inflated genes
        model.mat <- as.matrix(cbind(additive.mat, as.data.frame(exprs(sub))))
        genes.appear <- 1
    } 
    
    if(precenter)
        model.mat <- precenter.fun(model.mat)

    ## Holds output from glmnet
    fits <- vector(mode='list', length=2*length(genes))
    sigma <- nobs <- lambda <- lambda0 <- rep(NA, length=2*length(genes))
    dim(sigma) <- dim(nobs) <- dim(lambda) <- dim(lambda0) <- dim(fits) <- c(length(genes), 2)
    dimnames(sigma) <- dimnames(nobs) <- dimnames(lambda) <- dimnames(lambda0) <- dimnames(fits) <- list(primerid=genes, type=c('dichotomous', 'continuous'))

    ## Begin loop thru genes
    for(i in seq_along(genes)){
        this.gene <- fData(sub)$primerid[i]
        y.zif <- exprs(sub)[,i]
        y.dichot <- y.zif>0
        y.real <- exprs(sub)[,i][y.dichot]
        genes.diff <- setdiff(genes, this.gene)
        ## remove response gene from design
        this.gene.idx <- seq(from=i, by=length(genes), length.out=genes.appear)
        this.model <- model.mat[,-this.gene.idx-additive.dim]
        
        penalty.factor<-rep(c(0, 1), times=c(additive.dim, genes.appear*length(genes.diff)))
        if(any(this.gene %in% colnames(this.model))) stop('ruhroh')
        tt <- try({
            if(response == 'hurdle'){
            fit.dichot <- cv.glmnet(this.model, y.dichot, family='binomial', penalty.factor=penalty.factor, standardize=!precenter, ...)
        } else if(response=='cg.regression'){

        }else if(response == 'zero.inflated'){
            fit.dichot <- cv.glmnet(this.model, y.zif, family='gaussian', penalty.factor=penalty.factor, standardize=!precenter, ...)
            fit.real <- fit.dichot
            }
            fits[[i, 'dichotomous']] <- if(fit.dichot$glmnet.fit$jerr==-1 || min(fit.dichot$glmnet.fit$lambda) > 1e2 ) NULL else fit.dichot
            lambda[i,'dichotomous'] <- fit.dichot$lambda.min[1]
            lambda0[i, 'dichotomous'] <- fit.dichot$glmnet.fit$lambda[1]
            nobs.d <- nrow(this.model)
            nobs[i, 'dichotomous'] <- nobs.d
            my <- mean(y.dichot)
            sigma[i, 'dichotomous'] <- my*(1-my)

            nobs.c <- length(y.real)
            sigma.y <- var(y.real)
            nobs[i, 'continuous'] <- nobs.c 
            sigma[i, 'continuous'] <- sigma.y
            wy <- nobs.c/(nobs.d*sigma.y)
            if(response=='hurdle') fit.real <- cv.glmnet(this.model[y.dichot,], y.real, family='gaussian', penalty.factor=penalty.factor, standardize=!precenter, ...)
            ## set things to null if we didn't converge rather than return empty fit
            fits[[i, 'continuous']] <- if(fit.real$glmnet.fit$jerr==-1 || min(fit.real$glmnet.fit$lambda) > 1e2 ) NULL else fit.real
            lambda[i,'continuous'] <- fit.real$lambda.min[1]
            lambda0[i, 'continuous'] <- fit.real$glmnet.fit$lambda[1]
        })
        if(class(tt) == 'try-error') warning(sprintf('There was an error with gene %s', this.gene))
        message(this.gene, '\n')
    }

    structure(fits, genes=genes, gene.predictors=gene.predictors, additive.dim=additive.dim, lambda=lambda, lambda0=lambda0, nobs=nobs, sigma=sigma, response=response)
}

## Alternately: factor out looping code as function
## needs to write to lambda, lambda0, nobs, etc



fortify.zifnetwork <- function(fits, lc.range, ld.range, nknots=20, ebic.lambda=0){
    sigma <- rename(cast(melt(attr(fits, 'sigma')), primerid ~ type), c('continuous'='sigma.c', 'dichotomous' = 'sigma.d'))
    null <- is.na(attr(fits, 'nobs')[,1])
    cv.fit <- fits[!null,]
    genes <- attr(fits, 'genes')[!null]

    ## out.nativepath: use solution path specific to the gene
    ## out: line things up using knots over the l1 norm on the betas
    ## grp.norm.list: not currently used
    out.nativepath <- grp.norm.list <- out <- vector(mode='list', length=nrow(cv.fit))
    names(grp.norm.list) <- names(out) <- genes


    if(missing(lc.range) || missing(ld.range)){
         L.Dmax <- max(attr(fits, 'lambda0')[,1], na.rm=TRUE)
         L.Dmin <- quantile(attr(fits, 'lambda')[,1], na.rm=TRUE, probs=.1)

         L.Cmax <- max(attr(fits, 'lambda0')[,2], na.rm=TRUE)
         L.Cmin <- quantile(attr(fits, 'lambda')[,2], na.rm=TRUE, probs=.1)
         ld.range <- c(L.Dmin, L.Dmax)
         lc.range <- c(L.Cmin, L.Cmax)
         message(sprintf('Guessing `lc.range`=[%f, %f] and `ld.range`=[%f, %f]', L.Cmin, L.Cmax, L.Dmin, L.Dmax))
     }
    
    if(length(lc.range)!=2 || length(ld.range) != 2) stop("'lc.range' and 'ld.range' must both be length 2")
    
    knots.c <- seq(from=lc.range[1], to=lc.range[2], length=nknots)
    knots.d <- seq(from=ld.range[1], to=ld.range[2], length=nknots)
    
    for(g in seq_len(nrow(cv.fit))){
        twofit <- list(cv.fit[[g,1]]$glmnet.fit, cv.fit[[g,2]]$glmnet.fit)
        nobs.d <- twofit[[1]]$nobs
        nobs.c <- twofit[[2]]$nobs
        ndev.d <- (1-twofit[[1]]$dev.ratio)*twofit[[1]]$nulldev
        ndev.c <- (1-twofit[[2]]$dev.ratio)*twofit[[2]]$nulldev
        fixed.d <-twofit[[1]]$df[1]
        fixed.c <- twofit[[2]]$df[1]
        if(fixed.d != fixed.c) warning(sprintf('mismatch between length of fixed predictors in %s', genes[g]))
        norm.d <- apply(twofit[[1]]$beta, 2, function(x) sum(abs(x[-seq_len(fixed.d) ])))
        norm.c <- apply(twofit[[2]]$beta, 2, function(x) sum(abs(x[-seq_len(fixed.c) ])))

        l.d <- twofit[[1]]$lambda
        l.c <- twofit[[2]]$lambda

        nnz.d <- twofit[[1]]$df-fixed.d
        nnz.c <- twofit[[2]]$df-fixed.c

        bic.d <- ndev.d+nnz.d*log(nobs.d) + 2*ebic.lambda*nnz.d*log(length(genes))
        bic.c <- ndev.c+nnz.c*log(nobs.c) +2*ebic.lambda*nnz.c*log(length(genes))

        ## Now line things up by knots
        out.nativepath[[g]] <- rbind(data.frame(lambda=l.d, nnz=nnz.d, ndev=ndev.d, primerid=genes[g], bic=bic.d,component='discrete'),
                                     data.frame(lambda=l.c, nnz=nnz.c, ndev=ndev.c, primerid=genes[g], bic=bic.c,component='continuous'))

        norm.d <- approx(l.d, norm.d, knots.d, rule=2)$y
        norm.c <- approx(l.c, norm.c, knots.c, rule=2)$y
        ndev.d <- approx(l.d, ndev.d, knots.d, rule=2)$y
        ndev.c <- approx(l.c, ndev.c, knots.c, rule=2)$y

        ## coef.d <- as.matrix(coef(twofit[[1]], s=l.d))
        ## coef.c <- as.matrix(coef(twofit[[2]], s=l.c))
        ## things.to.estimate <- c('grp.l1', 'comb.l1', 'ndev.comb')
        ## comb.norm <- array(NA, dim=c(ncol(coef.d), ncol(coef.c), length(things.to.estimate)), dimnames=list(l.d=l.d, l.c=l.c, estimand=things.to.estimate))
        
        ## fixed.plus.intercept <- fixed.c+1
        ## for(i in seq_len(ncol(coef.d))){
        ##     for(j in seq_len(ncol(coef.c))){
        ##         cd <- coef.d[-seq_len(fixed.plus.intercept),i]
        ##         cc <- coef.c[-seq_len(fixed.plus.intercept),j]
        ##         comb.norm[i,j, 'grp.l1'] <- sum(sqrt(cc^2 + cd^2))
        ##         comb.norm[i, j, 'comb.l1'] <- sum(abs(cc)+abs(cd))
        ##         comb.norm[i,j,'ndev.comb'] <- ndev.d[i]+ndev.c[j]
        ##     }
        ## }
        ## grp.norm.list[[g]] <-cbind(cast(melt(comb.norm), ...~estimand, fun.aggregate='[', x=1), primerid=genes[g], stringsAsFactors=FALSE) #fun.aggregate='[' because we might have duplicate lambda is we're on the boundary

        
        nnz.d <- approx(l.d, nnz.d, knots.d, method='constant', rule=2)$y
        nnz.c <- approx(l.c, nnz.c, knots.c, method='constant', rule=2)$y
        l.d <- approx(l.d, l.d, knots.d, rule=2)$y
        l.c <- approx(l.c, l.c, knots.c, rule=2)$y

       
        #data.frame(norm=knots, ndev1, ndev2, l1, l2)
        out[[g]] <- data.frame(norm.d, norm.c, nnz.d, nnz.c, knots.d, knots.c, l.d, l.c, ndev.d, ndev.c, nobs.d, nobs.c, primerid=genes[g])
    }
    fortified <- merge(sigma, do.call(rbind, out))
    native.path <- do.call(rbind, out.nativepath)
    #norm.grid <- do.call(rbind, grp.norm.list)
    
    
    list(fortified=fortified, norm.grid=NA, native.path=native.path)
        
        }


##' Put coefficients from a set of network regressions into a matrix/array
##'
##' The array is ngenes X ngenes X {2,4}, with the last dimension depending
##' on whether zero-inflated or hurdle predictors were used.
##' @param listOfFits output from fitZifNetwork
##' @param l.c continuous lambda value
##' @param l.d discrete lambda value
##' @param collapse should the network be collapsed between layers?
##' @param union currently ignored
##' @param layers upon which layers in the listOfFits should we operate (eg, discrete, continuous or both)
##' @return an array
getZifNetwork <- function(listOfFits, l.c, l.d, collapse=FALSE, union=TRUE, layers){
    genes <- dimnames(listOfFits)[['primerid']]
    gene.predictors <- attr(listOfFits, 'gene.predictors')
    additive.dim <- attr(listOfFits, 'additive.dim')
    if(missing(layers)) layers <- seq_len(dim(listOfFits)[2])

    if(length(l.c) == 1){
        message("Taking 'l.c' to be constant over genes")
        l.c <- rep(l.c, length(genes))
    }

    if(length(l.d) == 1){
        message("Taking 'l.d' to be constant over genes")
        l.d <- rep(l.d, length(genes))
    }


    if(gene.predictors=='hurdle'){
         out <- array(0, c(length(genes), length(genes), 4), dimnames=list(genes, genes, c('di.cont', 'di.di', 'cont.di', 'cont.cont')))
         genes.appear <- 2
    } else if(gene.predictors=='zero.inflated'){
         out <- array(0, c(length(genes), length(genes), 2), dimnames=list(genes, genes, c('di', 'cont')))
         genes.appear <- 1
    }

    ## Genes X component

    lambdaList <- attr(listOfFits, 'lambda')
    
    for( i in seq_along(genes)){
        this.gene <- genes[i]
        genes.diff <- setdiff(genes, this.gene)
        for(j in layers){
            this.comp <- dimnames(listOfFits)[[2]][j]
            this.lambda <- if(this.comp =='continuous') l.c[i] else l.d[i]
            if(!is.null(listOfFits[[i, j]]) && length(this.lambda)>0){

                ## kill intercept
                coeft <- as.numeric(coef(listOfFits[[i,j]], s=this.lambda))[-(1:(additive.dim+1))]
                coeft <- ifelse(abs(coeft)<1e-3, 0, coeft)
                check <- sum(abs(coeft))
                #browser(expr=check>1.5*constraint)
                out[i,match(genes.diff, genes),(j-1)*genes.appear+1] <- coeft[1:(length(genes.diff))]
                if(gene.predictors=='hurdle')  out[i,match(genes.diff, genes),(j-1)*genes.appear+2] <- coeft[(length(genes.diff)+1):(2*length(genes.diff))]
            }
        }
    }
                                 # penalty.factor=penalty.factor)
        out[is.na(out)] <- 0
    if(collapse){
        ## Take union of layers 
        adj.nonsym <- apply(out!=0, c(1,2), any)*1
        return(adj.nonsym)
    }
    return(out)
    
}

color <- function(attr, palette=brewer.pal, ...){
    attr.num <- as.numeric(as.factor(attr))
    attr.col <- palette(max(attr.num), ...)[attr.num]
    structure(attr.col, attr=attr)
}

##' Turn an object from fitZifNetwork into a igraph object
##'
##' Edges/Vertex may be given colors by setting Vattr or Ettr
##' @param zifFit fitZifNetwork output
##' @param Vattr character vector named with vertex names. 
##' @param Eattr data.frame with columns 'X1' 'X2' naming vertex pairs and 'value' giving the value for this set of edges
##' @param collapse should the discrete/continuous parts be collapsed?
##' @param union within each layer, should we take the union of a vertex's neighborhood, or should we take the intersection?
##' @param weight currently ignored
##' @param ... passed to getZifNetwork
##' @return igraph object, with attributes
##' @import igraph
##' @importFrom RColorBrewer brewer.pal
layoutZifNetwork <- function(zifFit, Vattr=NULL, Eattr=NULL, collapse=TRUE, weight=FALSE, union=TRUE, ...){
    if(weight && collapse) stop("Cannot provide both 'weight' and 'collapse'")
    if(!collapse && !is.null(Eattr)) stop("Cannot provide edge attributes when 'collapse = FALSE'")
    
    adj <- SingleCellAnalysis:::getZifNetwork(zifFit, collapse=collapse, ...)
    layerFun <- if(union) function(layer) (layer+ t(layer))/2 else function(layer) layer*t(layer)
    totalEdges <- sum(abs(adj)>0)
    if(!collapse){
        stopifnot(length(dim(adj))==3)
        ## Symmetrize, get support and assign a bitmask to the layer
        for(l in seq_len(dim(adj)[3])){
            adj[,,l] <- layerFun(adj[,,l])
            adj[,,l] <- (abs(adj[,,l])>0)*2^(l) 
        }
        
        collapse <- aaply(adj, c(1,2), sum)
        Eattr <- subset(melt(collapse), value>0)
        Eattr$value <- factor(Eattr$value, levels=1:(2^(dim(adj)[3]+1)))
    } else{
        collapse <- abs(layerFun(adj))>0
    }

    connected <- colSums(collapse)>0
    collapse.connect <- (collapse[connected,][,connected]>0)*1
    gadj <- graph.adjacency(collapse.connect, mode = 'upper')
    el <- setNames(as.data.frame(get.edgelist(gadj)), c('X1', 'X2'))
    Ecol <- NULL
    if(!is.null(Eattr)){
    el2 <- merge(el, Eattr)
    Ecol <- color(el2$value, name='Paired')
    E(gadj)$color <- Ecol
}
    Vcol <- NULL
    if(!is.null(Vattr)){
        Vattr <- Vattr[V(gadj)$name]
        Vcol <- color(Vattr, name='Set1')
        V(gadj)$color <- Vcol
    }
    
       structure(gadj,Ecol=Ecol, Vcol=Vcol, Vlab=colnames(collapse.connect), totalEdges=totalEdges, adjacencyMatrix=collapse)
}

##' @importFrom plyr ddply
##' @import reshape
plotNetworksBIC <- function(FITS, ebic.lambda, layers=1:2, ...){

fort.out <- SingleCellAnalysis:::fortify.zifnetwork(FITS, ebic.lambda=ebic.lambda)$native.path

min.bic <- ddply(fort.out, ~primerid + component, function(df){
    df[which.min(df$bic), ]
    })

min.bic.pid <- cast(min.bic, primerid  ~ component, value='lambda')
min.bic.pid <- min.bic.pid[match(attr(FITS, 'genes'), min.bic.pid$primerid),]


SingleCellAnalysis:::layoutZifNetwork(FITS, collapse=FALSE, l.c=min.bic.pid$continuous + .001, l.d=min.bic.pid$discrete + .001, union=TRUE, layers=layers, ...)

}

plotNetworks <- function(FITS, nedges, layers=1:2, printNNZperLambda=TRUE, ...){

   
fort.out <- SingleCellAnalysis:::fortify.zifnetwork(FITS, nknots=100)
    
lambda.cv <- as.data.frame(attr(FITS, 'lambda'))
lambda.cv$primerid <- row.names(lambda.cv)
fort <- merge(fort.out$fortified, lambda.cv, by='primerid')

nnz.per.lambda <- cast(melt(fort, measure.vars=c('nnz.d','nnz.c')), knots.d + knots.c  ~ variable, fun.aggregate=sum)
nnz.per.lambda <- nnz.per.lambda[order(-nnz.per.lambda$knots.d),]

if(printNNZperLambda){
    print(ggplot(nnz.per.lambda)+geom_line(aes(y=knots.d, x=nnz.d), col='blue')+geom_line(aes(y=knots.c, x=nnz.c), col='red') + xlab('Edges') + ylab('Lambda'))
}


nnz.per.lambda.fun <- list(
    dichot=with(nnz.per.lambda, approxfun(nnz.d, knots.d)),
    cont=with(nnz.per.lambda, approxfun(nnz.c, knots.c)))

    grList <- list()
for(i in seq_along(nedges)){
    ## Spread edges evenly over the selected layers
    edges <- nedges[i]/length(layers)
    
    l.c <- nnz.per.lambda.fun$cont(edges)
    l.d <- nnz.per.lambda.fun$dichot(edges)
    
        grList[[i]] <- SingleCellAnalysis:::layoutZifNetwork(FITS, collapse=FALSE, l.c=l.c, l.d=l.d, union=TRUE, layers=layers, ...)
}
    grList
   }
