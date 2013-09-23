library(glmnet)
function(sc, additive.effects, min.freq=.05){
    sub <- sc[, freq(sc)>min.freq]
    genes <- fData(sub)$primerid
    out <- array(0, c(length(genes), length(genes), 4), dimnames=list(genes, genes, c('di.cont', 'di.di', 'cont.di', 'cont.cont')))
    additive.mat <- model.matrix(as.formula(sprintf('~ %s', paste(additive.effects, sep='+'))), cData(sub))[, -1] #no intercept
    additive.dim <- ncol(additive.mat)
    model.mat <- as.matrix(cbind(additive.mat,
                                 as.data.frame(xform(exprs(sub))), as.data.frame(exprs(sub)>0)))
    for(i in seq_along(genes)){
        this.gene <- fData(sub)$primerid[i]
        y.dichot <- exprs(sub)[,i]>0
        y.real <- exprs(sub)[,i][y.dichot]
        genes.diff <- setdiff(genes, this.gene)
        this.model <- model.mat[,-c(i, length(genes)+i)-additive.dim]
        penalty.factor<-rep(c(0, 1), times=c(additive.dim, 2*length(genes.diff)))
        if(any(this.gene %in% colnames(this.model))) stop('ruhroh')
        tt <- try({
            fit.dichot <- cv.glmnet(this.model, y.dichot, family='binomial', alpha=.9)
            coeft <- as.numeric(coef(fit.dichot))[-(1:(additive.dim+1))]}, silent=TRUE)    #kill intercept)
                                        # penalty.factor=penalty.factor)
        if(class(tt) == 'try-error') {warning(sprintf('There was an error with gene %s', this.gene));coeft<- rep(0, length(coeft))}
        fit.real <- cv.glmnet(this.model[y.dichot,], y.real, family='gaussian', alpha=.9)
        out[i,match(genes.diff, genes),1] <- coeft[1:(length(genes.diff))]
        out[i,match(genes.diff, genes),2] <- coeft[(length(genes.diff)+1):(2*length(genes.diff))]
        coeft <- as.numeric(coef(fit.real))[-(1:(additive.dim+1))]
        out[i,match(genes.diff, genes),3] <- coeft[1:(length(genes.diff))]
        out[i,match(genes.diff, genes),4] <- coeft[(length(genes.diff)+1):(2*length(genes.diff))]
    }
    out[is.na(out)] <- 0
}
