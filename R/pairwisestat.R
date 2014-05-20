gencor <- function(ct, maxmissingness = .99){
  #trim names
  colnames(ct) = substr(colnames(ct), 5, 20)
  #remove genes that were observed with lower frequency than 1-missingness
  numobs = apply(ct>0, 2, sum, na.rm=T)
  nullgenes = numobs < nrow(ct)*(1-maxmissingness)
  ct2 = ct[,!nullgenes]
  corbygene = cor(ct2, use="pairwise.complete.obs")
  corbygene = submatrix(corbygene) #remove genes without significant correlations
  heatmap.2(corbygene, trace="none", keysize=1)
  invisible(corbygene)
}

rmse.sim <- function(N=100){
p = c(.005, .01, .1)
out = expand.grid(p1=p, p2=p)
out = cbind(out,matrix(0, nrow(out), ncol=N))
reppp = vector()
for(i in 1:nrow(out)){
  for(n in 1:N){
    x<- rbinom(500, 1, p=out[i,1])
    y<-rbinom(500, 1, p=out[i,2])
    reppp[n] = rmse.test(x,y, P=5e3)$p
  }
out[i,3:(N+2)] =reppp
}
out
}


cap = function(covarmats, Nsamp){
  #Bartlet's modified likelihood ratio test
  p = nrow(covarmats[[1]])
  stopifnot(p==ncol(covarmats[[1]]))
  stopifnot(length(covarmats)==length(Nsamp))
  k = length(covarmats)
  covsum = do.call(sum, covarmats)
  Nsum = sum(Nsamp)
  dets = unlist(lapply(det, covarmats))
  detsum = det(covsum)
  prod(abs(dets)^((Nsamp-1)/2))*(Nsum - k)^(p*(Nsum-k)/2)/
    (abs(det(covsum))^((Nsum-k)/2)*prod((Nsamp-1)^(p*(Nsamp-1))) )
}

submatrix = function(mat, rowcrit=function(arow){sort(abs(arow), decreasing=T)[2]>.2 }){
  val = apply(mat, 2, rowcrit)
  mat1 = mat[val,]
  mat1[,val]
}

rmse.test <- function(x, y, P=5e4, nulldist=F){
  stopifnot(length(x)==length(y))
  N = length(x)
  v0 = matrix(c(sum(x*y), sum(x*(1-y)), sum((1-x)*y), sum((1-x)*(1-y))), nrow=4)
  out = .rmse(v0, N, expected=T)
  expHat = out[[2]]
  obs = out[[1]]
  v = rmultinom(P, N, expHat)
  h0 = .rmse(v, N)
  p = mean(obs <= h0)
  if(nulldist){
    return(list(p.value=p, obs=obs, null=h0))
  }
  return(list(p.value=p, obs=obs))
}

.rmse <- function(v, N, expected=F){
  #v is a 4xN matrix of simulated obs
  A = (v[1,]+v[2,])/N
  B = (v[1,]+v[3,])/N
  exp = rbind(A*B, A*(1-B), (1-A)*B, (1-A)*(1-B))*N
  out = sqrt(apply((v - exp)^2, 2, sum))
  if(!expected){
    return(out)
  }
    
  return(list(out, exp))
}

rmsepairws <- function(scr3, assays){
  pval = matrix(0, ncol=length(assays), nrow=length(assays), dimnames=list(assays, assays))
  oval = pval
  for(i in 1:(length(assays)-1)){
    for(j in (i+1):length(assays)){
      test = rmse.test(scr3[,assays[i]], scr3[,assays[j]])
      pval[i,j] = test$p
      oval[i,j] = test$obs
    }
  }
  
  pval = t(pval) + pval
  oval = t(oval) + oval
  
  return(list(pval, oval))
}

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

setClass('PairwisePearson', contains='matrix', representation=representation(nsamp='numeric', npos='numeric'), validity=function(object){
    ncol(object@.Data)==nrow(object@.Data) && ncol(object@.Data) == length(object@npos)
})

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

setClass('PairwiseOddsratio', contains='matrix', representation=representation(genes='character'))

setClass('PairwiseOddsWithHook', contains='PairwiseOddsratio', representation=representation(hook='matrix'))

##' Calculate pairwise gene odds ratios
##'
##' Calculates odds ratios between each pair of genes on the
##' dichotomous level.
##' @param sc SingleCellAssay object
##' @param partial a character vector of terms to adjust for.  \code{cData(sc)} is used as a model frame.  If contains 'ngeneson' and it is not present in \code{cData} then it will be calculated.
##' @return PairwiseOddsRatio class, with slots \code{fits} (2-D list of fit) and \code{genes} (character vector of gene names).
##' @import SingleCellAssay
##' @importFrom brglm brglm
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

setGeneric('coef')
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
