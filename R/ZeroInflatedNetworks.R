without <- function(obj, idx){
    if(length(idx)==0) return(obj )
    if(!is.null(dim(obj)) && length(dim(obj))>1) return(obj[,-idx,drop=FALSE]) else return(obj[-idx])
}

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
##' @param gene.predictors ignored
##' @param precenter How should centering/scaling be done with respect to continuous regressions.  TRUE if centering should be done with respect to all cells; FALSE if centering should be done only with respect to expressed cells
##' When precenter=TRUE, cv.glmnet will not standardize.
##' @param precenter.fun a function called to center the expression matrix prior to calling glmnet
##' @param response a character vector, one of 'zero.inflated', 'hurdle', or 'cg.regression', 'cg.regression2', 'cg.mle'
##' @param modelSelector a function called gene and component-wise on each fit, that should return an index to the glmnet lambda sequence for that gene and component
##' @param onlyReturnFitter if TRUE, return an undocumented fitter function that is internally called on each gene/component.
##' @param ... passed to cv.glmnet
##' @return 2-D list of cv.glmnet objects with attributes
##' @importFrom glmnet glmnet cv.glmnet
##' @import nloptr
##' @importFrom numDeriv jacobian
##' @export
fitZifNetwork <- function(sc, additive.effects, min.freq=.05, gene.predictors='zero.inflated', precenter=TRUE, precenter.fun=scale, response='hurdle', modelSelector, onlyReturnFitter=FALSE, debug=FALSE, ...){
    ## gene.predictors <- match.arg(gene.predictors, c('zero.inflated', 'hurdle'))
    response <- match.arg(response, c('hurdle', 'zero.inflated', 'cg.regression',  'cg.regression2', 'cg.mle'))
    sub <- sc[, freq(sc)>min.freq]
    genes <- fData(sub)$primerid
    ngenes <- length(genes)

    ## Additive (un-penalized) variables named from cData
    additive.mat <- model.matrix(as.formula(sprintf('~ %s', paste(additive.effects, collapse='+')), env=parent.frame()), cData(sub))[, -1, drop=FALSE] #no intercept
    additive.dim <- ncol(additive.mat)

    ## transform expression and generate design
    expr <- exprs(sub)
    if(response %in% c('cg.regression', 'cg.regression2')){
        expr <- xform(expr)
    }

    model.mat <- as.matrix(cbind(additive.mat, as.data.frame(expr)))
    model.mat.zero <- as.matrix(cbind(additive.mat, as.data.frame(1*(exprs(sub)>0))))

    if(precenter){
        model.mat <- precenter.fun(model.mat)
        model.mat.zero <- precenter.fun(model.mat.zero)
    }

    ## Holds output from glmnet
    fits <- vector(mode='list', length=2*length(genes))
    sigma2 <- nobs <- lambda <- lambda0 <- rep(NA, length=2*length(genes))
    dim(sigma2) <- dim(nobs) <- dim(lambda) <- dim(lambda0) <- dim(fits) <- c(length(genes), 2)
    dimnames(sigma2) <- dimnames(nobs) <- dimnames(lambda) <- dimnames(lambda0) <- dimnames(fits) <- list(primerid=genes, type=c('dichotomous', 'continuous'))
    pf<-rep(c(0, 1), times=c(additive.dim, ngenes-1))

    ## called for each gene/component
    glmnetFit <- function(y.zif, this.gene, this.model, this.model.zero, component, fits, lambda, sigma2, ...){
        y.dichot <- y.zif>0
        y.real <- y.zif[y.dichot]
        if(component=='continuous'){
            family <- 'gaussian'
        } else{
            family <- 'binomial'
        }

        off <- NA
        ## Gaussian
        if(family=='gaussian' && response %in% 'hurdle'){
            fit <- glmnet(this.model[y.dichot,], y.real, family=family, penalty.factor=pf, standardize=!precenter, ...)
            nobs <- length(y.real)
        } else if(family=='gaussian' && response == 'cg.regression2'){
            pf <- c(pf, rep(1, ngenes))
            fit <- glmnet(cbind(this.model, without(this.model.zero, seq_len(additive.dim)))[y.dichot,], y.real, family=family, standardize=!precenter, penalty.factor=pf, ...)
            nobs <- length(y.real)
        } else if(family=='gaussian' && response == 'cg.regression'){
            fit <- glmnet(this.model[y.dichot,], y.real-mean(y.real), family=family, penalty.factor=pf, standardize=!precenter, ...)
            nobs <- length(y.real)
        }else if(family=='gaussian' && response =='zero.inflated'){
            fit <- glmnet(this.model, y.zif, family=family, penalty.factor=pf, standardize=!precenter, ...)
            nobs <- length(y.zif)
        } else if(family=='binomial' && response == 'hurdle'){
            fit <- glmnet(this.model, y.dichot, family=family, penalty.factor=pf, standardize=!precenter, ...)
            nobs <- length(y.dichot)
        } else if(family=='binomial' && response == 'zero.inflated'){
            fit <- fits[[this.gene, 'continuous']]
            nobs <- length(y.dichot)
        } else if(family=='binomial' && response %in% c('cg.regression', 'cg.regression2')){
            nobs <- length(y.dichot)
            fit.c <- fits[[this.gene, 'continuous']]
            if(is.null(fit.c)) stop('Empty continuous fit')
            l.cont <- lambda[this.gene, 'continuous']
            Kbb <- 1/sigma2[this.gene, 'continuous']

            ## calculate offset given lambda
            coef.c <- as.numeric(coef(fit.c, s=l.cont))
            newx <- if(response =='cg.regression') this.model else cbind(this.model, without(this.model.zero, seq_len(additive.dim)))
            ## H[b|a]/Kbb
            cont.fitted <- predict(fit.c, s=l.cont, newx=newx)
            #Hb <- coef.c[1]*Kbb    #Intercept times precision
            #Kba <- as.matrix(-coef.c[-1]*Kbb) #others times -precision

            ## H[b|a]^2*/Kbb^2
            off <- Kbb*cont.fitted^2/2
            #fitted2 <- this.model %*% Kba
            #off2 <- (fitted2^2/2 -fitted2)/(Kbb*Hb)
            offt <- off-mean(off)
            TOP <- 2
            offt[offt < -TOP] <- -TOP
            offt[offt > TOP] <- TOP
            offt <- offt-mean(offt)
            message(summary(offt))
            df <- data.frame(fit=drop(cont.fitted), pos=factor(y.dichot), off=drop(offt))
            aplot <- ggplot(df, aes(x=fit, col=pos))+geom_density(adjust=3)
            bplot <- ggplot(df,aes(x=off, y=fit, col=pos)) + geom_point()
            ## browser(expr=this.gene=='CXCL1')

            if(response == 'cg.regression'){
                thisx <- this.model.zero
            }else{
                thisx <- cbind(this.model.zero, without(this.model, seq_len(additive.dim)))
                pf <- c(pf, rep(1, ngenes))
            }
            fit <- glmnet(thisx, y.dichot, family=family, penalty.factor=pf, offset=offt, standardize=!precenter, ...)
        }
        l.idx <- modelSelector(fit, ngenes=ngenes)
        sigma2 <- (1-fit$dev.ratio)*fit$nulldev/nobs
        ## Need to wrap fit into object/method
        structure(fit, nobs=nobs, sigma2=sigma2[l.idx], selectedLambda=fit$lambda[l.idx], off=off)
    }

    if(response=='cg.mle'){
        stopifnot(additive.dim==0)
        th0 <- rep(0, ngenes*4-1)
        names(th0) <- parmap(ngenes)
        th0['hbb'] <- th0['kbb'] <- 1
        th0['gbb'] <- -10
        lb <- setNames(rep(-Inf, length(th0)), names(th0))
        lb['kbb'] <- .001
        ## th0['gbb'] <- -13
        ## th0['hbb'] <- 5

        glmnetFit <- function(y.zif, this.gene, this.model, this.model.zero, j, fits, lambda, sigma2, ...){
            if(j=='dichotomous'){
                fit <- fits[[this.gene, 'continuous']]
                return(fit)
            }
            ll <- generatelogLik(y.zif, this.model, debug=debug, ...)
            #oo <- optim(th0, ll, gr, method='BFGS',control=list(maxit=4000), hessian=TRUE)
            oo2 <- nloptr(th0, ll, opts=list(algorithm='NLOPT_LD_TNEWTON_PRECOND_RESTART', maxeval=200, check_derivatives=debug, check_derivatives_print='errors'), lb=lb)
            sol <- setNames(oo2$solution, names(lb))

            gr <- function(th) ll(th)$gradient
            hess <- jacobian(gr, oo2$sol)
            fit <- list(coefficients=sol, jerr=oo2$status, lambda=0, hess=hess)
            structure(fit, nobs=length(y.zif), sigma2=1/sol['kbb'], selectedLambda=0,
                      genes=c(this.gene, colnames(this.model), #gbb, gba
                          this.gene, colnames(this.model), #hbb, hab
                          colnames(this.model), colnames(this.model),#hba, kba   
                          this.gene ))  #kbb
        }
    }


    if(onlyReturnFitter) return(glmnetFit)
    
    ## Loop through components
    for(j in c('continuous', 'dichotomous')){
        ## Begin loop thru genes
        for(i in seq_along(genes)){
            this.gene <- fData(sub)$primerid[i]
            message(this.gene)
            y.zif <- exprs(sub)[,i]
            ## remove response gene from design
            this.gene.idx <- i
            this.model <- model.mat[,-this.gene.idx-additive.dim,drop=FALSE]
            this.model.zero <- model.mat.zero[,-this.gene.idx-additive.dim,drop=FALSE]
            if(any(this.gene %in% colnames(this.model))) stop('ruhroh')
            this.fit <- tryCatch({
                glmnetFit(y.zif, this.gene, this.model, this.model.zero, j, fits, lambda, sigma2, ...)
            }, error=function(e) SingleCellAssay:::reraise(e, convertToWarning=!debug))
            ## recover from glmnet errors or non-convergence
            if(inherits(this.fit, 'error') || this.fit$jerr==-1 ||  any(!is.finite(this.fit$lambda)) || min(this.fit$lambda) > 1e2){
                warning(sprintf('There was an error with gene %s', this.gene))
            }else{
                fits[[i,j]] <- this.fit
                lambda[i,j] <- attr(this.fit, 'selectedLambda')
                lambda0[i,j] <- this.fit$lambda[1]
                nobs[i, j] <- attr(this.fit, 'nobs')
                sigma2[i, j] <- attr(this.fit, 'sigma2')
            }
        }                               #end gene loop
    }                                   #end component loop
    
    structure(fits, genes=genes, gene.predictors=gene.predictors, response=response,additive.dim=additive.dim, lambda=lambda, lambda0=lambda0, nobs=nobs, sigma2=sigma2, response=response, class=c('FittedZifNetwork', 'class'))
}

nullSelector <- function(fit, ngenes){
    c(l.idx=1)
}

bicSelector <- function(fit, ngenes, ebic.lambda=1){
    l <- fit$lambda
    fixed <- fit$df[which.max(l)]
    nnz <- fit$df-fixed
    ndev <- (1-fit$dev.ratio)*fit$nulldev
    bic <- ndev+nnz*log(fit$nobs) + 2*ebic.lambda*nnz*log(ngenes)
    l.idx <- which.min(bic)
    c(l.idx=l.idx)
}


##' Derive global properties of network fit
##'
##' .. content for \details{} ..
##' @param fits \code{fitZifNetwork} object
##' @param lc.range optional
##' @param ld.range optional
##' @param nknots 
##' @param ebic.lambda lambda penalty for extended Bayesian Info Crit. (Rigel and Drton)
##' @return list with entries \code{fortified}, \code{norm.grid}, \code{native.path}
##' fortified and native.path are both data.frames with entries for each primerid, containing statistics of the fit as lambda varies.
##' fortified has cartesian product of continuous lambda values and discrete, over the same grid for each gene.
##' native.path has 
##' @import reshape
##' @import data.table
fortify.zifnetwork <- function(fits, lc.range, ld.range, nknots=20, ebic.lambda=0){
    #sigma <- rename(cast(melt(attr(fits, 'sigma2')), primerid ~ type), c('continuous'='sigma.c', 'dichotomous' = 'sigma.d'))
    null <- is.na(attr(fits, 'nobs')[,1]) & is.na(attr(fits, 'nobs')[,2])
    cv.fit <- fits[!null,]
    genes <- attr(fits, 'genes')[!null]

    ## out.nativepath: use solution path specific to the gene
    ## out: line things up using knots over the l1 norm on the betas
    ## grp.norm.list: not currently used
    out.nativepath <- grp.norm.list <- out <- vector(mode='list', length=nrow(cv.fit))
    names(grp.norm.list) <- names(out) <- genes



    if(missing(lc.range) || missing(ld.range)){
         L.Dmax <- max(attr(fits, 'lambda0')[,1], na.rm=TRUE)
         L.Dmin <- quantile(attr(fits, 'lambda')[,1], na.rm=TRUE, probs=.5)

         L.Cmax <- max(attr(fits, 'lambda0')[,2], na.rm=TRUE)
         L.Cmin <- quantile(attr(fits, 'lambda')[,2], na.rm=TRUE, probs=.5)
         ld.range <- c(L.Dmin, L.Dmax)
         lc.range <- c(L.Cmin, L.Cmax)
         message(sprintf('Guessing `lc.range`=[%f, %f] and `ld.range`=[%f, %f]', L.Cmin, L.Cmax, L.Dmin, L.Dmax))
     }
    
    if(length(lc.range)!=2 || length(ld.range) != 2) stop("'lc.range' and 'ld.range' must both be length 2")
    
    knots <- list(disc=seq(from=ld.range[1], to=ld.range[2], length=nknots),
                  cont=seq(from=lc.range[1], to=lc.range[2], length=nknots))
    comp <- c('disc', 'cont')

    fout <- CJ(primerid=genes, component=c('disc', 'cont'), norm=NA_real_, nnz=NA_real_, l=NA_real_, ndev=NA_real_,  nobs=NA_real_)
    setkey(fout, component)
    kout <- data.table(knots=c(knots[['disc']],knots[['cont']]),
                       component=rep(c('disc', 'cont'), each=nknots))
    setkey(kout, component)
    fout <- merge(fout, kout, allow.cartesian=TRUE)
    setkey(fout, component, primerid)
    
    np <- CJ(primerid=genes, component=c('disc', 'cont'), lseq=seq_along(fits[[1,2]]$lambda), lambda=NA_real_, nnz=NA_real_, ndev=NA_real_, bic=NA_real_)    
    setkey(np, primerid, component, lseq)
    
    for(g in seq_len(nrow(cv.fit))){
        twofit <- list(cv.fit[[g,1]], cv.fit[[g,2]])
        for(j in 1:2){
        nobs.d <- twofit[[j]]$nobs
        if(is.null(nobs.d)) next
        
          ndev.d <- (1-twofit[[j]]$dev.ratio)*twofit[[j]]$nulldev
          fixed.d <-twofit[[j]]$df[1]
          norm.d <- apply(twofit[[j]]$beta, 2, function(x) sum(abs(without(x, seq_len(fixed.d)))))
  
        l.d <- twofit[[j]]$lambda
  
        nnz.d <- twofit[[j]]$df-fixed.d
  
        bic.d <- ndev.d+nnz.d*log(nobs.d) + 2*ebic.lambda*nnz.d*log(length(genes))

        this.seq <- as.integer(seq_along(l.d))

        np[list(primerid=genes[g], component=comp[j], lseq=this.seq), c('lambda', 'nnz', 'ndev', 'bic'):= list(l.d, nnz.d, ndev.d, bic.d), by=NULL]


        #np[primerid==genes[g]& component==component[j] & lseq==this.seq, c('lambda', 'nnz', 'ndev', 'bic'):= list(l.d, nnz.d, ndev.d, bic.d)]

        ## Now line things up by knots

        this.knots <- knots[[j]]
        norm.d <- approx(l.d, norm.d, this.knots, rule=2)$y
          ndev.d <- approx(l.d, ndev.d, this.knots, rule=2)$y
          nnz.d <- approx(l.d, nnz.d, this.knots, method='constant', rule=2)$y
          l.d <- approx(l.d, l.d, this.knots, rule=2)$y
        #data.frame(norm=knots, ndev1, ndev2, l1, l2)

         fout[J(comp[j], genes[g]), c('norm', 'nnz', 'l', 'ndev', 'nobs'):= list(norm.d, nnz.d, l.d, ndev.d, nobs.d), by=NULL]
    }
    }
    #fortified <- merge(sigma, do.call(rbind, out))

    fortified <- fout[!is.na(norm),]
    native.path <- np[!is.na(bic),]
    list(fortified=fortified, norm.grid=NA, native.path=native.path)
        
        }

## lof = list of fits from fitZifNetwork
coefLayer <- function(lof, s, layer){
    genes <- attr(lof, 'genes')
    ngenes <- length(genes)
    add <- attr(lof, 'additive.dim') +1 #intercept
    stopifnot(length(s)==ngenes)
    stopifnot(length(layer)==1)
    out <- matrix(0, nrow=ngenes, ncol=ngenes, dimnames=list(genes, genes))
    if(is.integer(layer)){
        layer <- c('cont', 'disc', 'cont2', 'disc2')[layer]
        warning('Assuming layer is ', paste(layer, collapse=','))
    }

    ngenesNotSelf <- ngenes-1#because response gene is always omitted
    coefIdx <- 1:(ngenesNotSelf)+add  
    if(layer=='cont'){
        comp <- 'continuous'
        } else if(layer=='cont2'){
            comp <- 'continuous'
            stopifnot(attr(lof, 'response')=='cg.regression2')
            coefIdx <- coefIdx + ngenesNotSelf
        } else if(layer=='disc'){
            comp <- 'dichotomous'
        } else if(layer=='disc2'){
            comp <- 'dichotomous'
            stopifnot(attr(lof, 'response')=='cg.regression2')
            coefIdx <- coefIdx+ngenesNotSelf
        }

    for(i in seq_along(genes)){
        try({
        co <- coef(lof[[i,comp]], s=s[i])
        rn <- row.names(co)
        co <- setNames(as.numeric(co), rn)
        stopifnot(all(genes[-i] == names(co)[coefIdx]))
        out[,i][-i] <- co[coefIdx]
            })
    }
    out
}

##' Put coefficients from a set of network regressions into a matrix/array
##'
##' The array is ngenes X ngenes X {2,3}, with the last dimension depending
##' on whether zero-inflated or cg.regression2 predictors were used.
##' @param listOfFits output from fitZifNetwork
##' @param l.c continuous lambda value, if missing use lambda attribute from listofFits
##' @param l.d discrete lambda value, see l.c
##' @param collapse should the network be collapsed between layers?
##' @param union currently ignored
##' @param layers upon which layers in the listOfFits should we operate (eg, discrete, continuous or both)
##' @return an array
getZifNetwork <- function(listOfFits, l.c, l.d, collapse=FALSE, union=TRUE, layers){
    genes <- dimnames(listOfFits)[['primerid']]
    response <- attr(listOfFits, 'response')
    additive.dim <- attr(listOfFits, 'additive.dim')
    lambda <- attr(listOfFits, 'lambda')

    if(missing(l.c)){
        l.c <- lambda[,'continuous']
    } else if(length(l.c)==1){
        message("Taking 'l.c' to be constant over genes")
        l.c <- rep(l.c, length(genes))
     }

    if(missing(l.d)){
        l.d <- lambda[,'dichotomous']
    }else if(length(l.d)==1){
        message("Taking 'l.d' to be constant over genes")
        l.d <- rep(l.d, length(genes))
    }

    if(missing(layers)) layers <- c('cont', 'disc')
    out <- array(0, c(length(genes), length(genes), length(layers)), dimnames=list(genes, genes, layers))
    
        for(j in layers){
            this.lambda <- if(j %in% c('cont', 'cont2')) l.c else l.d
            out[,,j] <- coefLayer(listOfFits, s=this.lambda, layer=j)
        }

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
plotNetworksBIC <- function(FITS, layers=1:2, ...){
SingleCellAnalysis:::layoutZifNetwork(FITS, collapse=FALSE, union=TRUE, layers=layers, ...)
}

plotNetworks <- function(FITS, nedges, layers=1:2, printNNZperLambda=TRUE, ...){

   
fort.out <- SingleCellAnalysis:::fortify.zifnetwork(FITS, nknots=100)
    
lambda.cv <- as.data.frame(attr(FITS, 'lambda'))
lambda.cv$primerid <- row.names(lambda.cv)
#fort <- merge(fort.out$fortified, lambda.cv, by='primerid', all.y=TRUE)
fort <- fort.out$fortified
nnz.per.lambda <- fort[,list(nnz=sum(nnz)), by=list(component, knots)]

if(printNNZperLambda){
    print(ggplot(nnz.per.lambda, aes(y=knots, x=nnz, col=component))+geom_line() + xlab('Edges') + ylab('Lambda') + xlim(0, 500))
}


nnz.per.lambda.fun <- list(
    dichot=with(nnz.per.lambda[component=='disc'], approxfun(nnz, knots)),
    cont=with(nnz.per.lambda[component=='cont'], approxfun(nnz, knots)))


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
