##' Chi-square tests of independence for each pair of genes
##'
##' Conduct a chi-square tests of independence for each pair of genes.
##' The binary (expressed/not) values are used to build the contingency tables.
##' The chi-square test statistics for each pair of genes is returned, with zeroes placed along the main diagonal.
##' Significance could be assessed (in large samples) by comparison to a chi-square distribution with 1-dof.
##' @param sc SingleCellAssay, on a thresholded (et) layer
##' @return matrix of chisq test statistics
##' @author andrew
##' @export
chisq.pairwise <- function (sc) {
  if(layername(sc) != 'et') warning('Tests fail unless on a thresholded layer')
  scr3 <- exprs(sc)==0
  marg = apply(scr3, 2, mean)
  onames = names(marg)
  #onames = substr(onames, 5, 30)
  exp=array(dim=c(length(marg),length(marg), 4), dimnames=list(onames, onames, NULL))
  exp[,,1] = outer(marg, marg)*nrow(scr3)
  exp[,,2] = outer(1-marg, marg)*nrow(scr3)
  exp[,,3] = outer(marg, 1-marg)*nrow(scr3)
  exp[,,4] = outer(1-marg, 1-marg)*nrow(scr3)
  
  
  xy=exp
  xy[,,1] = t(as.matrix(scr3)) %*% as.matrix(scr3)
  xy[,,2] = t(as.matrix(!scr3)) %*% as.matrix(scr3)
  xy[,,3] = t(as.matrix(scr3)) %*% as.matrix(!scr3)
  xy[,,4] = t(as.matrix(!scr3)) %*% as.matrix(!scr3)
  
  chi = (xy-exp - .5*sign(xy-exp))^2/(exp+.0001) #continuity correction, and avoid division by 0
  chi2 = apply(chi, c(1, 2), sum)
  diag(chi2) = 0 #of course X_i is not independent of X_i!
  #smallexp = apply(exp<5, c(1, 2), any)
  #chi2[smallexp] = NA
  chi2
}

##' Calculate pairwise correlation coefficient
##'
##' Calculates the correlation coefficient between each pair of genes on the
##' dichotomous level.
##' The diagonal is set to zero to avoid having the fisher transform blow up.
##' @param sc SingleCellAssay object
##' @param diagonals currently ignored
##' @param partial adjust for ngeneson
##' @return a PairwisePearson object, which is a matrix with slots "nsamp" giving the number of cells the cor. coef. were calculated over, "npos" giving the number of expressed cells in each gene.
##' @import SingleCellAssay
##' @export
pearson.pairwise <- function(sc, diagonals='zero', partial=FALSE){
    if(layername(sc) != 'et') warning('Tests fail unless on a thresholded layer')
    ee <- (exprs(sc)>0)*1
    ngeneson <- apply(ee, 1, mean)
    if(partial){
        ee <- aaply(ee, 2, function(x){
            resid(glm(x~ngeneson, family='binomial'))
                  })
        ee <- t(ee)
    }
 r2 <- cor(ee)
 diag(r2) <- 0
    new('PairwisePearson', r2, nsamp=nrow(sc), npos=apply(exprs(sc)>0, 2, sum))
}

fisher.z <- function(r) .5 * (log(1+r)-log(1-r))
inv.fisher <- function(z) (exp(2*z)-1)/(1+exp(2*z))

##' Generate a matrix of z-statistics comparing two PairwisePearson objects
##'
##' 
##' @param pp1 First PairwisePearson object
##' @param pp2 Second PairwisePearson object
##' @return matrix of z-statistics
##' @importFrom abind abind
##' @export
test.pearson.pairwise <- function(pp1, pp2){
    z1 <- fisher.z(pp1)
    z2 <- fisher.z(pp2)
    zdiff <- (z2-z1)/sqrt(sum(1/(pp1@nsamp-3), 1/(pp2@nsamp-3)))
    diag(zdiff) <- 0
    zdiff
}

ci.pearson.pairwise <- function(pairwisePearson, ci=.95, method='asymptotic'){
    z <- fisher.z(pairwisePearson)
    norm.q <- qnorm((1-ci)/2, lower.tail=FALSE)/sqrt(pairwisePearson@nsamp-3)
    ci.lower <- inv.fisher(z-norm.q)
    ci.upper <- inv.fisher(z+norm.q)
    abind(lower=ci.lower, est=pairwisePearson@.Data, upper=ci.upper, along=0)
}


##' Calculate pairwise gene odds ratios
##'
##' Calculates odds ratios between each pair of genes on the
##' dichotomous level.
##' @param sc SingleCellAssay object
##' @param Formula a formula used to build the design matrix.  Probably should include the term 'logOR.primer2' (if the pairwise odds ratios are indeed desired).
##' @param surrenderFreq any pair of genes in which either is expressed above or below this threshold will be skipped
##' @param useFirth should biased-reduced regression be used via brglm.  Slow.
##' @param lm.hook a function to be called on the fitted regression
##' @param MAX_IT if useFirth==TRUE, how many iterations should be permitted?
##' @return PairwiseOddsRatio class, with slots \code{fits} (2-D list of fit) and \code{genes} (character vector of gene names).
##' @import SingleCellAssay
##' @importFrom brglm brglm brglm.control
##' @export
oddsratio.pairwise <- function(sc, Formula, surrenderFreq=.05, useFirth=FALSE, lm.hook, MAX_IT=30){
    if(layername(sc) != 'et') warning('Tests fail unless on a thresholded layer')
    ee <- (exprs(sc)>0)*1
    ## if(!missing(Formula)){
    ##     cd <- cData(sc)
    ##     if(str_detect(as.character(Formula)[2], 'ngeneson') && !('ngeneson'%in% names(cData(sc))) ){
    ##         ngeneson <- apply(ee, 1, mean)
    ##         cd$ngeneson <- ngeneson
    ##     }
    ##     mf <- model.frame(Formula, cd)
    ## } else{
    ##     mf <- data.frame('(Intercept)'=rep(1, nrow(ee)))
    ## }

    mf <- cData(sc)

   genes <- fData(sc)$primerid
   freqs <- freq(sc)
   fits <- vector(mode='list', length=length(genes)^2)
   dim(fits) <- c(length(genes), length(genes))
   dimnames(fits) <- list(primerid1=genes, primerid2=genes)

    if(!missing(lm.hook)){
       hook <- fits
       dim(hook) <- dim(fits)
       dimnames(hook) <- dimnames(fits)
   }


   EPS <- 10000 * .Machine$double.eps
   for(i in seq_along(genes)){
       if(freqs[i] < surrenderFreq || freqs[i] > 1-surrenderFreq) next

       mf.new <- cbind(mf, logOR.primer2=ee[,i])
       Xnew <- model.matrix(Formula, mf.new)
       message(genes[i])
       for(j in setdiff(seq_along(genes), i)){
           if(freqs[j]<surrenderFreq || freqs[j]>1-surrenderFreq) next
           tt <- try({
               if(useFirth){
                   if(genes[i] %in% c('GZMA', 'IL-21'))
                       message(genes[j])
                   #browser(expr=genes[i]=='GZMA')
                   #ll <- logistf(ee[,j] ~ .+0, data=as.data.frame(Xnew), pl=FALSE)
                   ll <- brglm(ee[,j] ~ .+0, data=as.data.frame(Xnew), control.brglm=brglm.control(br.maxit=MAX_IT))
                   ll$converged <- ll$nIter<MAX_IT
                   vcov <- vcov(ll)
                   rownames(vcov) <- colnames(vcov) <- names(coef(ll))
                   ll$boundary <- FALSE
               } else{
                   ll <- glm.fit(Xnew, ee[,j], family=binomial())
                   ll$boundary <- ll$boundary | any((ll$fitted < EPS) | (ll$fitted > 1-EPS))
                   Ahat <- crossprod(Xnew*ll$weights, Xnew)
                   trA <- sum(diag(Ahat))
                   p <- ncol(Ahat)
                   vcov <- solve(Ahat)
               }
               
           })
           if(!inherits(tt, 'try-error')){
               fits[[i,j]] <- list(coef=coef(ll), vcov=vcov, converged=ll$converged, boundary=ll$boundary, primer1=genes[i], primer2=genes[j])
               if(!missing(lm.hook))
                   hook[[i, j]] <- lm.hook(ll)
           }
       }
   }

    if(!missing(lm.hook)){
        return(new('PairwiseOddsWithHook', fits, genes=genes, hook=hook))
    } else{
    return(new('PairwiseOddsratio', fits, genes=genes))
}
}

safeSubset <- function(list, index){
    !is.null(list[[index]]) && list[[index]]
}

setMethod(coef, 'PairwiseOddsratio', function(object, onlyConverged=FALSE, se=FALSE){
    tmp <- lapply(object, function(x) x[['coef']])
    tmp.names <- expand.grid(primer1=object@genes, primer2=object@genes, stringsAsFactors=FALSE)
    converged <- rep(TRUE, length(tmp))
    if(onlyConverged){
        converged <- sapply(object, safeSubset, 'converged') & !sapply(object, safeSubset, 'boundary')
    }
    empty <- sapply(tmp, is.null) | !converged
    m <- do.call(rbind, tmp[!empty])
    tmp.names <- tmp.names[!empty,]
        if(se){
        se.list <- lapply(object, function(x) sqrt(diag(x[['vcov']])))
        se.list <- se.list[!empty]
        se.bind <- do.call(rbind, se.list)
        res <- cbind(m, se=se.bind, tmp.names)
    } else{
    res <- cbind(m, tmp.names)
}
    return(res)

 })
