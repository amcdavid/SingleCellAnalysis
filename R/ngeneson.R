calcNgeneson <- function(sca){
    ee <- exprs(sca)
    ngeneson <- apply(ee>1, 2, mean)
    zmean <- apply(ee, 2, mean)
    ee[ee==0] <- NA
    cmean <- apply(ee, 2, mean, na.rm)
    cd <- cData(sca)
    alreadySet <- intersect(names(cd), c('ngeneson', 'condmean', 'zmean'))
    if(length(alreadySet)>0) message('Overwriting field(s) ', paste(alreadySet, collapse = ','), ' in cData.')
    cd$ngeneson <- ngeneson
    cd$condmean <- cmean
    cd$zmean <- zmean
    cData(sca) <- cd
    sca
}

##' Get coefficients and standard errors from a zlm fit
##'
##' For each coefficient name in cname, return discrete and continuous values and standard errors
##' @param zlm zlm object
##' @param cname character vector of coefficient names
##' @return 3D array: (disc, cont) X (coefficient, standard error) X coefficient
getCoef <- function(zlm, cname){
    res <- array(NA, dim=c(2, 2, length(cname)), dimnames=list(component=c('disc', 'cont'), arg=c('coef', 'se'), coef=cname))
    for(component in c('disc', 'cont')){
        res[component,'coef',] <- coef(zlm[[component]])[cname]
        res[component, 'se',] <- sqrt(diag(vcov(zlm[[component]])))[cname]
    }
        res
    }

examineNgeneson <- function(zlm, sca){
    model.info <- laply(zlm$models, getCoef, cname='ngeneson')
    dimnames(model.info)[[1]] <- names(zlm$models)
    names(dimnames(model.info))[1] <- 'primerid'
    this.freq <- freq(sca)
    this.condmean <- condmean(sca)
    gene.info <- data.frame(primerid=names(this.freq), freq=this.freq, condmean=this.condmean)
    merge(melt(model.info), gene.info)
}
