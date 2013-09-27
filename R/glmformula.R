##' @importFrom glmnet cv.glmnet
setGeneric('cv.glmnet', function(x, y, ...) standardGeneric('cv.glmnet'))

contr.dummy <- function(n, base=1, contrasts=FALSE, sparse=FALSE){
  contr.treatment(n, base, contrasts, sparse)
}

##' Generic for cv.glmnet
##'
##' Accepts formula arguments
##' @param x formula
##' @param y data.frame or environment in which \code{x} is evaluated
##' @param ... arguments passed to cv.glmnet
##' @return see cv.glmnet
##' @export
setMethod('cv.glmnet', signature='formula', function(x, y, ...){
  family <- list(...)$family
  subset <- list(...)$subset
  if(is.null(family)) family <- 'gaussian'
  stopifnot(is.data.frame(y))
  if(is.null(subset)){
    mf <- model.frame(x, y, subset=NULL)
  } else {
    mf <- model.frame(x, y, subset=subset)
  }
  
  response <- model.response(mf)
  opar <- options('contrasts')
  npar <- opar
  npar[[1]][1] <- 'contr.dummy'
  options(npar)
  mm <- model.matrix(x, mf)
  options(opar)
  if(family=='binomial')
    response <- as.numeric(as.factor(response))
  x <- mm
  y <- response
  callNextMethod()
}
          )

## Run lasso on singlecellassay object to predict comparison
## Omit genes above/below min.freq
## Use either continuous, dichotomous, or both predictors
## If using both, take orthogonal xform of continuous
## Return crossvalidated glm, model matrix, response and singlecellassay object
glmSingleCellAssay <- function(sca, comparison, min.freq, predictor=c('continuous', 'dichotomous'), addn, addn.penalty, user.mm, alpha=.9, only.mm=FALSE, ...){
  predOpts <-  c('dichotomous', 'continuous', 'interaction', 'user')
  predictor <- match.arg(predictor, predOpts, several.ok=TRUE)
  sel <- freq(sca) > min.freq & (1-freq(sca)) > min.freq
  ngenes <- sum(sel)
  #setkeyv(melt(sca),getMapping(sca,"idvars"))
  ee <- exprs(sca)
  ee <- ee[,sel]
  ee.dichot <- ee>0
  if(all(c('dichotomous', 'continuous') %in% predictor)){
  ee.real <- xform(ee)
} else{
  ee.real <- ee
}
  df <- as.data.frame(cbind(ee.dichot, ee.real))
  names(df) <- make.names(c(paste(colnames(ee), 'd', sep='.'), paste(colnames(ee), 'c', sep='.')))
  
  ## select indices of desired predictors
  idx <- list(1:ngenes, (ngenes+1):(2*ngenes) )[match(predictor, predOpts)]
  idx <- do.call(c, idx)
  df <- df[,idx]
 if('interaction' %in% predictor){
  form <- sprintf('~ (%s)^2 + 0', paste(names(df), collapse="+"))
} else{
 form <- sprintf('~ %s + 0', paste(names(df), collapse="+"))
} 
  resp <- as.factor(cData(sca)[,comparison])
  fam <- 'binomial'
  if(length(levels(resp))>2)
    fam <- 'multinomial'

  opar <- options('contrasts')
  npar <- opar
  npar[[1]][1] <- 'contr.dummy'
  options(npar)
  if('user' %in% predictor){
      mm <- user.mm(df)
  } else{
  mm <- model.matrix(formula(form), df)
}
  if(!missing(addn)){
  mm <- cbind(mm, addn)
}
  pf <- rep(1, ncol(mm))
  if(!missing(addn.penalty)){
  pf <- c(rep(1, ncol(mm)), addn.penalty)
}
  options(opar)

  mm <- mm[!is.na(resp),]
  resp <- resp[!is.na(resp)]
  
  if(only.mm) return(list(mm=mm, resp=resp))
  
  fit <- cv.glmnet(mm, resp,  family=fam, alpha=alpha, penalty.factor=pf, ...)
  list(cv.fit=fit, mm=mm, response=resp, sca=sca)
}

## predict response in glmsca using model model matrix in glmsca
## unless alt.fit (a cv.glmnet object) is provided
## in which case, response is predicted via the fit in cv.glmnet
glmMisclass <- function(glmsca, alt.fit, groups, s='lambda.min'){
  stopifnot(missing(groups)|| length(groups)==1)
  if(!missing(alt.fit)){
  pred <- predict(alt.fit$cv.fit, glmsca$mm, type='response', s=s)[,,1]
  nresp <- levels(alt.fit$response)
} else{
   pred <- predict(glmsca$cv.fit, glmsca$mm, type='response', s=s)[,,1]
   nresp <- levels(glmsca$response)
}
  resp <- glmsca$response
  w.max <- apply(pred, 1, which.max)
  if(missing(groups)){fact <- rep('(all)', length(glmsca$resp))}
  else{ fact <- cData(glmsca$sca)[,groups,drop=FALSE]}
  names(fact) <- NULL
  truthvec <- data.frame(pred=nresp[w.max], truth=resp, fact=fact)
  tab <- with(truthvec, table(pred, truth, fact))
  ## tab <- tapply(seq_along(nrow(pred)), nrow(fact), function(i){
  ##   x <- table(pred[i,], levels(resp)[w.max[i]])
  ##               print(x)
  ##               x
  ##             })
  print(tab)
  correct.class <- apply(tab, 3, function(x){1-sum(diag(prop.table(x)))})
  message('Error rate ', paste(names(correct.class), round(correct.class, 3), sep=':', collapse=' '))
}

.sparseGlmToMat <- function(sparseList, ...){
  mat <- lapply(sparseList, as.matrix)
  mm <- melt(mat)
  coefmat <- cast(mm,  X1 ~ L1, ...)
  coefmat <- rename(coefmat, c('X1'='predictor'))
  coefmat[is.na(coefmat)] <- 0
  coefmat
}

summarizeCoef <- function(glmsca, s='lambda.min'){
  co <- coef(glmsca$cv.fit, s=s)
  names(co) <- levels(glmsca$response)
  coefmat <- .sparseGlmToMat(co, subset=value>0)
  traj <- .sparseGlmToMat(glmsca$cv.fit$glmnet.fit$beta, fun.aggregate=function(x){sum(x>0)})
  class(coefmat) <- c('CoefficientMatrix', class(coefmat))
  origcandidate <- str_split_fixed(coefmat$predictor, '[.][dc]$', 2)[,1]
  hasCandidate <- origcandidate %in% colnames(exprs(glmsca$sca))
  origcandidate <- origcandidate[hasCandidate]
  path.enter <- subset(traj, predictor %in% coefmat$predictor)
  path.enter2 <- apply(path.enter[hasCandidate,][,-1], 1, max)
  path.order <- order(tapply(path.enter2, origcandidate, max))
  attr(coefmat, 'path.order') <- path.order
  attr(coefmat, 'originalvar') <- unique(origcandidate)
  coefmat
}

plot.CoefficientMatrix <- function(x, sca){
  ee <- exprs(sca)[,attr(x)$originalvar]
  heat2(ee, Colv=attr(x)$path.order)
}
