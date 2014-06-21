## setGeneric('cv.glmnet',  function(x, y, ...) standardGeneric('cv.glmnet'))

## contr.dummy <- function(n, base=1, contrasts=FALSE, sparse=FALSE){
##   contr.treatment(n, base, contrasts, sparse)
## }


## setMethod('cv.glmnet', signature='ANY', function(x, y, ...){
##     callNextMethod()
## })

## ##' Generic for cv.glmnet
## ##'
## ##' Accepts formula arguments
## ##' @param x formula
## ##' @param y data.frame or environment in which \code{x} is evaluated
## ##' @param ... arguments passed to cv.glmnet
## ##' @return see cv.glmnet
## ## @export
## setMethod('cv.glmnet', signature='formula', function(x, y, ...){
##   family <- list(...)$family
##   subset <- list(...)$subset
##   if(is.null(family)) family <- 'gaussian'
##   stopifnot(is.data.frame(y))
##   if(is.null(subset)){
##     mf <- model.frame(x, y, subset=NULL)
##   } else {
##     mf <- model.frame(x, y, subset=subset)
##   }
  
##   response <- model.response(mf)
##   mm <- model.matrix(x, mf)
##   if(family=='binomial')
##     response <- as.numeric(as.factor(response))
##   x <- mm
##   y <- response
##   callNextMethod()
## }
##           )

##' Run a multinomial lasso on a SingleCellAssay object to predict group membership
##'
##' This function generates a design matrix based on the expression values in \code{sca}
##' and calls \code{cv.glmnet} to try to classify a group named by \code{comparison}, which keys a column in the \code{cData} of \code{sca}
##'
##' The design matrix is generated according to the option \code{predictor}. If \code{predictor} vector includes the term 'dichotomous', then each gene is treated as binary indicators.  If the term 'continuous' is included, then the zero-inflated (continuous) value for the gene is used.  If both 'continuous' and 'dichotomous' are included, then the both values for the gene are used, however the continuous values are centered about their conditional mean using the function \code{xform}.  If 'interaction' is included, then all the terms are crossed with each other to generate pairwise interactions.
##' 
##' @param sca SingleCellAssay object
##' @param comparison character naming a column in \code{cData(sca)}
##' @param min.freq minimum frequency for a gene to be considered in the classifier
##' @param predictor character vector naming some combination of 'continuous', 'dichotomous' or 'interaction'.  See details.
##' @param pen.scale.interaction multiply the l1 penalty by this factor if interactions are included
##' @param precenter should the gene predictors be centered? Recommended if there are interactions present to reduce co-linearity of the interaction with the marginal term.
##' @param prescale should the gene predictors be scaled to have unit variance?
##' @param addn character vector, giving additional columns of design, interpreted in the context of cData(sca)
##' @param addn.penalty an optional numeric giving the relative scale of the penalty for add
##' @param user.mm a function to be applied to exprs(sca) instead of the defaults given by \code{predictor}
##' @param alpha elasticnet penalty parameter. Default =.9.
##' @param only.mm Should only the model matrix be returned, rather than actually calling cv.glmnet?
##' @param ... additional arguments to cv.glmnet.
##' @return list with components 'cv.fit' giving the output from cv.glmnet, 'mm' giving the model matrix, 'response' giving the response vector and 'sca' containing the 'sca' passed as input to the function
##' @seealso glmMisclass, getNZdesign, doGLMnet, cv.glmnet
##' @export
##' @importFrom Matrix cBind
##' @import glmnet
##' @import stringr
glmSingleCellAssay <- function(sca, comparison, min.freq, predictor=c('continuous', 'dichotomous'), pen.scale.interaction=2, precenter=FALSE, prescale=FALSE, addn=NULL, addn.penalty=0, user.mm, alpha=.9, only.mm=FALSE, ...){
  predOpts <-  c('dichotomous', 'continuous', 'interaction', 'user')
  predictor <- match.arg(predictor, predOpts, several.ok=TRUE)
  sel <- freq(sca) > min.freq & (1-freq(sca)) > min.freq
  ngenes <- sum(sel)
  #setkeyv(melt(sca),getMapping(sca,"idvars"))
  ee <- exprs(sca)
  ee <- ee[,sel]
  ee.dichot <- ee>0
  if(all(c('dichotomous', 'continuous') %in% predictor)){
  ee.real <- SingleCellAssay:::xform(ee)
} else{
    ee.real <- scale(ee, scale=prescale, center=precenter)
        
}
  df <- as.data.frame(cbind(ee.dichot, ee.real))
  names(df) <- make.names(c(paste(colnames(ee), 'd', sep='.'), paste(colnames(ee), 'c', sep='.')))
  
  ## select indices of desired predictors
  idx <- list(1:ngenes, (ngenes+1):(2*ngenes) )[match(predictor, c('dichotomous', 'continuous'))]
  idx <- do.call(c, idx)
  df <- df[,idx]
 if('interaction' %in% predictor){
  form <- sprintf('~ (%s)^2', paste(names(df), collapse="+"))
} else{
 form <- sprintf('~ %s', paste(names(df), collapse="+"))
} 
  resp <- as.factor(cData(sca)[,comparison])
  fam <- 'binomial'
  if(length(levels(resp))>2)
    fam <- 'multinomial'


  ## df <- cbind(df, cData(sca))
  
  if('user' %in% predictor){
      mm <- user.mm(df)[,-1]
  } else{
  mm <- model.matrix(formula(form), df)[, -1]
}
  
  pf <- rep(1, ncol(mm))
  if('interaction' %in% predictor){
      interIdx <- str_detect(colnames(mm), fixed(':'))
      pf[interIdx] <- pf[interIdx]*pen.scale.interaction
  }


    if(!is.null(addn)){
      addn <- model.matrix(formula(paste0('~', addn)), cData(sca))
      pf <- c(pf, rep(addn.penalty, ncol(addn)))
      additive.dim <- seq(from=ncol(mm)+1, length.out=ncol(addn))
      mm <- cbind(mm, addn)
} else{
    additive.dim <- 0
}



  mm <- mm[!is.na(resp),]
  resp <- resp[!is.na(resp)]
  
  if(only.mm) return(list(mm=mm, resp=resp))
  message('Calling cv.glmnet on ',  ncol(mm), ' predictors and ', nrow(mm), ' cells')
  fit <- cv.glmnet(mm, resp,  family=fam, alpha=alpha, penalty.factor=pf, ...)
  list(cv.fit=fit, mm=mm, response=resp, sca=sca, additive.dim=additive.dim)
}

## predict response in glmsca using model model matrix in glmsca
## unless alt.fit (a cv.glmnet object) is provided
## in which case, response is predicted via the fit in cv.glmnet
glmMisclass <- function(glmsca, alt.fit, groups, s='lambda.min'){
  stopifnot(missing(groups)|| length(groups)==1)
  if(missing(alt.fit))
      alt.fit <- glmsca
  
  nresp <- levels(glmsca$response)
  resp <- glmsca$response
  if(missing(groups)){fact <- rep('(all)', length(glmsca$resp))}
  else{ fact <- cData(glmsca$sca)[,groups,drop=FALSE]}
  names(fact) <- NULL
  pred <- factor(as.vector(predict(glmsca$cv.fit, alt.fit$mm, type='class')), levels=nresp)
  truthvec <- data.frame(pred=pred, truth=resp, fact=fact)
  tab <- with(truthvec, table(pred, truth, fact))
  ## tab <- tapply(seq_along(nrow(pred)), nrow(fact), function(i){
  ##   x <- table(pred[i,], levels(resp)[w.max[i]])
  ##               print(x)
  ##               x
  ##             })

  correct.class <- apply(tab, 3, function(x){1-sum(diag(prop.table(x)))})
  message('Error rate ', paste(names(correct.class), round(correct.class, 3), sep=':', collapse=' '))
  invisible(list(confusion=tab, errorRate=correct.class))
}

.sparseGlmToMat <- function(glmsca, s='lambda.1se', additive=TRUE){
    scores <- coef(glmsca$cv.fit, s=s)
    Scores <- do.call(cBind, scores)[-1,] #No intercept
    if(!additive && all(glmsca$additive.dim>0)){
        Scores <- Scores[-glmsca$additive.dim,]
}
    nz <- apply(abs(Scores)>1e-3, 1, any)
    Scores <- as.matrix(Scores[nz,])
    colnames(Scores) <- names(scores)
    Scores
}

getNZDesign <- function(glmsca, s='lambda.1se', additive=TRUE){
    scores <- .sparseGlmToMat(glmsca, s=s, additive)
    nz <- row.names(scores)
    glmsca$mm[,nz]
}

summarizeCoef <- function(glmsca, s='lambda.min'){
  co <- coef(glmsca$cv.fit, s=s)
  names(co) <- levels(glmsca$response)
  coefmat <- .sparseGlmToMat(co)
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
