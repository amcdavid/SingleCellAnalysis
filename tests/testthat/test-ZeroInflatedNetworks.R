library(SingleCellAssay)
data(vbetaFA)
vbetaFA <- subset(vbetaFA, ncells==1)
ngeneson <- data.frame(ngeneson=apply(exprs(vbetaFA)>0, 1, mean))
ntest <- 4
vbetaFA <- vbetaFA[,freq(vbetaFA)>.1]
vbetaFA <- combine(vbetaFA, ngeneson)
vbetaT <- vbetaFA[,seq_len(ntest)]

getNullDF <- function(fit){
unlist(lapply(fit, function(fit) fit$df[1]))
}


nullout <- fitZifNetwork(vbetaT, additive.effects='0', precenter.fun=function(x) scale(x, scale=FALSE), response='zero.inflated', modelSelector=nullSelector)
context('Testing network fits')
test_that('Can fit', {
expect_is(nullout, 'FittedZifNetwork')
expect_equivalent(nullout[,1], nullout[,2])
expect_equal(attr(nullout, 'additive.dim'),0)
expect_true(all(getNullDF(nullout)==0))

outHurdle <- fitZifNetwork(vbetaT, '0', precenter.fun=function(x) scale(x, scale=FALSE), response='hurdle', modelSelector=nullSelector)
expect_is(all.equal(outHurdle[,1], outHurdle[,2], check.attributes=FALSE, use.names=FALSE), 'character')

})

test_that('Can include unpenalized covariate', {
out <- fitZifNetwork(vbetaT, 'ngeneson', precenter.fun=function(x) scale(x, scale=FALSE), response='zero.inflated', modelSelector=nullSelector)
expect_is(out, 'FittedZifNetwork')
expect_equal(attr(out, 'additive.dim'),1)
expect_true(all(getNullDF(out)==1))

out <- fitZifNetwork(vbetaT, c('factor(Stim.Condition)', 'ngeneson'), precenter.fun=function(x) scale(x, scale=FALSE), response='zero.inflated', modelSelector=nullSelector)
expect_equal(attr(out, 'additive.dim'),2)
expect_true(all(getNullDF(out)==2))

out <- fitZifNetwork(vbetaT, c('factor(Stim.Condition)', 'ngeneson'), precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression', modelSelector=nullSelector)
expect_equal(attr(out, 'additive.dim'),2)
expect_true(all(getNullDF(out)==2))

out <- fitZifNetwork(vbetaT, c('factor(Stim.Condition)', 'ngeneson'), precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression2', modelSelector=nullSelector)
expect_equal(attr(out, 'additive.dim'),2)
#expect_true(all(getNullDF(out)==2))


})

test_that('Throw error on bad unpenalized covariate', {
    expect_error(fitZifNetwork(vbetaT, 'BADDDD', precenter.fun=function(x) scale(x, scale=FALSE), response='zero.inflated', modelSelector=nullSelector))
    expect_error(fitZifNetwork(vbetaT, 'ngenes', precenter.fun=function(x) scale(x, scale=FALSE), response='zero.inflated', modelSelector=nullSelector))
})


context('Testing propriety of estimates')
test_that('Variance is calculated properly', {
    unitout <- fitZifNetwork(vbetaT, '0', precenter.fun=function(x) scale(x, scale=FALSE), response='hurdle', modelSelector=nullSelector)
    s <- attr(unitout, 'sigma2')
    ee <- exprs(vbetaT)
    exprs(vbetaT) <- ee*10
    scaleout <- fitZifNetwork(vbetaT, '0', precenter.fun=function(x) scale(x, scale=FALSE), response='hurdle', modelSelector=nullSelector)
    sscale <- attr(scaleout, 'sigma2')
    expect_equal(s[,'dichotomous'], sscale[,'dichotomous'])
    expect_equal(s[,'continuous']*100, sscale[,'continuous'])
})

context('Testing derived design matrices')
test_that('Matrices are scaled properly',{
    hurdle <- fitZifNetwork(vbetaT, 'ngeneson', precenter.fun=function(x) scale(x, scale=FALSE), response='hurdle', modelSelector=nullSelector, onlyReturnFitter=TRUE)
    model.mat <- get('model.mat', environment(hurdle))
    model.mat.zero <- get('model.mat.zero', environment(hurdle))
    expect_equivalent(apply(model.mat, 2, mean), rep(0, ntest+1))
    expect_equivalent(apply(model.mat.zero, 2, mean), rep(0, ntest+1))

    cg <- fitZifNetwork(vbetaT, 'ngeneson', precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression2', modelSelector=nullSelector, onlyReturnFitter=TRUE)
    model.mat <- get('model.mat', environment(cg))
    model.mat.zero <- get('model.mat.zero', environment(cg))
    expect_equivalent(apply(model.mat, 2, mean), rep(0, ntest+1))
    expect_equivalent(apply(model.mat.zero, 2, mean), rep(0, ntest+1))

    expressed1 <- abs(model.mat[,-1])>1e-7
    expressed2 <- model.mat.zero[,-1]>0
    expect_equal(expressed1, expressed2)
})


context('Testing cgRegression code')
test_that('Can fit', {
    cg <- fitZifNetwork(vbetaT, '0', precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression2', modelSelector=bicSelector)
    expect_is(cg, 'FittedZifNetwork')
    cg <- fitZifNetwork(vbetaT, '0', precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression', modelSelector=bicSelector)


    ## Add tests for translation invariance
    
})


context('Testing error handling')
test_that('Not too zealous in squelching errors', {
    expect_error(cg <- fitZifNetwork(vbetaFA, 'ngeneson', precenter.fun=function(x) x <- matrix('abcd', nrow=nrow(vbetaFA), ncol=ncol(vbetaFA)), min.freq=0, response='hurdle', modelSelector=bicSelector, lambda.min.ratio=0, debug=TRUE), regexp='NA/NaN/Inf')
    expect_warning(cg <- fitZifNetwork(vbetaFA, 'ngeneson', precenter.fun=function(x) x <- matrix('abcd', nrow=nrow(vbetaFA), ncol=ncol(vbetaFA)), min.freq=0, response='hurdle', modelSelector=bicSelector, lambda.min.ratio=0, debug=FALSE))
    
})

context('Testing network fortification')
cg <- fitZifNetwork(vbetaFA, 'ngeneson', precenter.fun=function(x) scale(x, scale=FALSE), response='cg.regression', modelSelector=bicSelector, lambda.min.ratio=.05)
hurdle <-fitZifNetwork(vbetaFA, 'ngeneson', precenter.fun=function(x) scale(x, scale=FALSE), response='hurdle', modelSelector=bicSelector, lambda.min.ratio=.05)
test_that('Can fortify', {
    f <- fortify.zifnetwork(cg, ebic.lambda=1)
    expect_is(f$fortified, 'data.frame')
    expect_is(f$native.path, 'data.frame')
    fmax <- subset(f$fortified, knots>=max(knots) & component=='cont')
    ngenes <- length(attr(cg, 'genes'))
    expect_equal(nrow(fmax), ngenes)
    expect_equal(fmax$nnz, rep(0, ngenes))
})

context('Testing network getting')
test_that('Can get layer', {
    g <- coefLayer(cg, s=rep(.1, ncol(vbetaFA)), layer='cont')
    expect_is(g, 'matrix')
    expect_false(all(g==0))
    expect_equal(dim(g),c(ncol(vbetaFA), ncol(vbetaFA)))
                   
    ## g <- coefLayer(cg, s=rep(.1, ncol(vbetaFA)), layer='cont2')
    ## expect_is(g, 'matrix')
    ## expect_false(all(g==0))

    g <- coefLayer(cg, s=rep(.05, ncol(vbetaFA)), layer='disc')
    expect_is(g, 'matrix')
    expect_false(all(g==0))

    g <- coefLayer(hurdle, s=rep(c(0,3), c(1,ncol(vbetaFA)-1)), layer='disc')
    expect_is(g, 'matrix')
    expect_true(all(is.na(g[,-1]) | g[,-1]==0))
    expect_false(all(g[,1]==0))
})

test_that('Can get array', {
    nc <- ncol(vbetaFA)
    g <- getZifNetwork(cg, layer=c('cont', 'disc'))
    expect_is(g, 'array')
    expect_false(all(g==0))
    expect_equal(dim(g),c(nc, nc, 2))
                   

    g <- getZifNetwork(hurdle, layer=c('cont', 'disc'))
    expect_is(g, 'array')
    expect_false(all(g==0))
    expect_equal(dim(g),c(nc, nc, 2))
    })



library(numDeriv)
context('MV hurdle MLE')
test_that('Finite difference gradients approximate analytic', {
    ee <- exprs(vbetaT)    
    gl <- generatelogLik(ee[,1], ee[,-1, drop=FALSE])
    dr <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], returnGrad=TRUE)
    par <- get('par', environment(gl))
    th <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, th), dr(th), tolerance=1e-4)

    gl <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], lambda=1)
    dr <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], returnGrad=TRUE, lambda=1)
    par <- get('par', environment(gl))
    th <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, th), dr(th),tolerance=1e-4)
})
