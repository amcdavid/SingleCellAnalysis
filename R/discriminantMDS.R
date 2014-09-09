##' Draw a confidence ellipse
##'
##' @param mapping 
##' @param data 
##' @param geom 
##' @param position 
##' @param ... 
##' @return ggplot2 layer
##' @importFrom MASS cov.trob
##' @import proto
stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
  StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
}


StatEllipse <- proto(ggplot2:::Stat,
	{
		required_aes <- c("x", "y")
		default_geom <- function(.) GeomPath
		objname <- "ellipse"

		calculate_groups <- function(., data, scales, ...){
			.super$calculate_groups(., data, scales,...)
		}
		calculate <- function(., data, scales, level = 0.75, segments = 51,...){
      dfn <- 2
      dfd <- length(data$x) - 1
      if (dfd < 3){
      	ellipse <- rbind(c(NA,NA))	
      } else {
          v <- cov.trob(cbind(data$x, data$y))
          shape <- v$cov
          center <- v$center
          radius <- sqrt(dfn * qf(level, dfn, dfd))
          angles <- (0:segments) * 2 * pi/segments
          unit.circle <- cbind(cos(angles), sin(angles))
          ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
      }
    
      ellipse <- as.data.frame(ellipse)
      colnames(ellipse) <- c("x","y")
      return(ellipse)
		}
	}
)


##' Run Linear Discriminant Analysis
##'
##' Runs Fisher's LDA to discriminate 'grouping' via 'train',
##' optionally attempting out-of-sample prediction using sample 'predict'
##' @param train \code{SingleCellAssay}
##' @param predict \code{SingleCellAssay}
##' @param grouping \code{character} vector of length 1 naming a column in the \code{cData} of 'train' and 'predict'
##' @param Layer which layer should be used to train and predict?
##' @return object of class 'LinearClassifier' and 'matrix'
##' @export
##' @seealso \link{doGLMnet}, \link{annotateBiPlot}
##' @importFrom MASS lda
doLDA <- function(train, predict, grouping, Layer='lCount'){
        layer(train) <- Layer
        ll <- lda(exprs(train), grouping=cData(train)[,grouping])
        pX <- predict(ll, exprs(predict))$x
        scores <- coef(ll)
        colnames(scores) <- colnames(pX) <- paste('Dim', seq_len(ncol((pX))), sep='')
        cls <- structure(pX, cData=cData(predict), scores=scores)
        class(cls) <- c('LinearClassifier', 'matrix')
        prepClassifier(cls)
}
##' Evaluate a glmSingleCellAssay classifier
##'
##' Wraps a glmnet classifier (derived from glmSingleCellAssay) to get the scores (coefficients) and predictions.
##' @param train output from \code{glmSingleCellAssay}, used for its fitted coefficients
##' @param predict output from \code{glmSingleCellAssay}, used for its model.matrix to predict (possibly) out of sample
##' @param s glmnet penalty parameter, default 'lambda.1se'
##' @return object of class 'Linear Classifier' and 'matrix'
##' @seealso \link{doLDA}, \link{annotateBiPlot}, \link{glmSingleCellAssay}
##' @export
doGLMnet <- function(train, predict, s='lambda.1se'){
    pX <- predict(train$cv.fit, predict$mm, s=s)
    Scores <- .sparseGlmToMat(train, s=s, additive=FALSE)
    cls <- structure(drop(pX), cData=cData(predict$sca), scores=Scores)
    class(cls) <- c('LinearClassifier', 'matrix')
    prepClassifier(cls)
}

##' Generate biplot vectors from a 'LinearClassifier'
##'
##' If there are G discrminant functions and P genes,
##' Generate  'biNorm', an array of G x G x P, giving the 2-dimensional 2-norm of each gene along each 2-combination of groups
##' 'biGenes', also G x G x P giving the ordering of the genes along the 2-norms,
##' 'norm', the G-dimensional 2-norm of each gene.
##' @param lc 'LinearClassifier'
##' @return 'LinearClassifer' with attributes
prepClassifier <- function(lc){
    scores <- attr(lc, 'scores')
    features <- row.names(scores)
    ranks <- apply(scores, 2, rank)
    rankScores <- setNames(merge(melt(scores), melt(ranks), by=c('X1', 'X2')), c('primerid', 'dimension', 'score', 'rank'))
    
    scoreNorm <- cast(rankScores,primerid ~ ., value='score', fun.aggregate=function(x) sum(x^2))
    scoreNorm$primerid <- as.character(    scoreNorm$primerid)

    ## Get 2-norm for each predictor along each pair of dimensions
    biDist <- biNorm <- array(NA, dim=c(ncol(scores), ncol(scores), nrow(scores)), dimnames=list(colnames(scores), colnames(scores), features))
    ## features ordered by 2-norm for this pair of dimensions
    biDistGenes <- biGenes <- array('', dim=dim(biNorm))
    dimnames(biDistGenes)[1:2] <- dimnames(biGenes)[1:2] <- dimnames(biNorm)[1:2]
    for(d1 in seq_len(ncol(scores)-1)){
        for(d2 in seq(d1+1, ncol(scores))){
            biNorm[d1, d2,] <- biNorm[d2,d1,] <- sqrt(scores[,d1]^2+scores[,d2]^2)
            biDist[d1, d2,] <- biDist[d2,d1,] <- (scores[,d1]-scores[,d2])^2
            biGenes[d1, d2,] <- features[order(-biNorm[d1,d2,])]
            biDistGenes[d1,d2,] <- features[order(-biDist[d1,d2,])]
        }
    }

    structure(lc, biNorm=biNorm, biGenes=biGenes, norm=scoreNorm, biDist=biDist, biDistGenes=biDistGenes)
}

##' Combine predictions from a LinearClassifier with cData from the object used to generate it
##'
##' @param lc LinearClassifier
##' @param optional ignored
##' @param row.names ignored
##' @param ... ignored
##' @return data.frame
##' @export
as.data.frame.LinearClassifier <- function(lc, row.names, optional, ...){
    cbind(as.data.frame.matrix(lc), attr(lc, 'cData'))
}

getCoord <- function(lc, d1, d2, genesToShow, metric='norm'){
    metric <-  match.arg(metric, c('norm', 'distance'))
    bn <- if(metric=='norm') attr(lc, 'biNorm') else attr(lc, 'biDist')
    bg <- if(metric=='norm') attr(lc, 'biGenes') else attr(lc, 'biDistGenes')
    genesSelected <- dim(bg)[3]
    thisGenes <- bg[d1,d2,1:min(genesToShow,genesSelected)]
    as.data.frame(attr(lc, 'scores')[thisGenes,c(d1, d2)])
}

##' Add an confidence ellipses or other features to a ggpairs object
##'
##' This function is necessary because there is no way to set aesthetics in ggpairs for only a portion of the plot matrix.
##' If ellipseArgList is empty, then no ellipses will be drawn.
##' If globalGGobjs is non-empty, then it will be added (+) to each panel.  You might use this to alter the color scale 
##' @param ggpairsObj output from ggpairs
##' @param panels numeric vector giving rows/columns to modify from the _lower_ triangle of the plot matrix
##' @param ellipseArgList list of ggplot aesthetics to be set to fixed values for the ellipses, eg, \code{lwd} or \code{alpha}
##' @param pointsArgList list of ggplot aesthetics to be rewritten for the purposes of calling geom point, eg, \code{size} or \code{alpha}
##' @param globalGGobjs vector of ggplot2 theme or scale elements 
##' @return ggpairs object with modified panels
##' @export
##' @import ggplot2
addEllipse <- function(ggpairsObj, panels, ellipseArgList=list(lwd=1, alpha=1), pointsArgList=list(), globalGGobjs=list()){
    if(length(pointsArgList)>0)         update_geom_defaults('point', pointsArgList)
    ncol <- length(ggpairsObj$columns)
    for(p1 in seq_len(ncol)){
        for( p2 in seq_len(ncol)){
            panel <- getPlot(ggpairsObj, p1, p2)
            if(length(ellipseArgList)>0 && p2<p1 && p1<=max(panels)) panel <- panel + do.call(stat_ellipse, ellipseArgList)
            if(length(globalGGobjs)>0)
            ggpairsObj <- putPlot(ggpairsObj, panel+ globalGGobjs, p1, p2)
        }
    }
    ggpairsObj
}

##' Add bi-plot vectors onto a ggpairs scatter plot of classifiers
##'
##' Draws vectors representing the "angle" and "length" the top 'genesToShow' predictors have when trying to discminant objects.
##' Vectors are drawn on the lower triangle of a ggpairs object (probably plotted by calling ggpairs(as.data.frame(lc), ...)) by matching names in the ggpairs object to the names of the classifier functions in lc.
##' So this function probably won't work unless ggpairs was called on 'lc'.
##' The 'top' genes can be selected either by the 2-norm (x^2+y^2)^(1/2), or the their 2-distance (x-y)^2 in each panel.
##' There is an arbitrary scaling between the scatter plot and the length of these vectors.  By default the longest vector is scaled to have length of the narrowest plot limit, but this can be adjusted via 'expand'.
##' @param ggpairsObj output from a call to 'ggpairs'
##' @param lc object of class 'LinearClassifier'
##' @param metric character, one of 'distance' or 'norm'.  See details.
##' @param genesToShow number of bi-plot vectors to show?
##' @param expand scaling factor of bi-plot vectors
##' @param where character vector, one or more of "upper" or "lower"
##' @param debug 
##' @param ... additional arguments passed to ggplot
##' @return modified ggpairs object, which can be plotted by evaluating it.
##' @importFrom GGally getPlot putPlot
##' @importFrom grid unit arrow
##' @export
##' @seealso \link{doLDA}, \link{doGLMnet}
annotateBiPlot <- function(ggpairsObj, lc, metric='norm', genesToShow=5, expand=1, where='lower', debug=FALSE, ...){
    ## Precondition: lower triangle contains scatter plots
    ## and possibly some condition on the order of the scatter plots and lc
    dims <- match(colnames(lc), names(ggpairsObj$data)[ggpairsObj$columns])
    where <- match.arg(where, c('upper', 'lower'), several.ok=TRUE)
    metric <- match.arg(metric, c('norm', 'distance'))
    for(d1Idx in seq_along(dims)[-length(dims)]){ #gives index in terms of lc
        d1 <- dims[d1Idx]               #gives index in terms of ggpairs
        for(d2Idx in seq(d1Idx+1, length(dims))){
            d2 <- dims[d2Idx]
            gp <- getPlot(ggpairsObj, d2, d1)
            panel <- ggplot_build(gp)$panel
            rx <- min(abs(panel$ranges[[1]]$x.range))
            ry <- min(abs(panel$ranges[[1]]$y.range))
            rangeMin <- min(rx, ry)
            gc <- getCoord(lc, d1Idx, d2Idx, genesToShow, metric)
            gcMax <- max(abs(gc))
            scale <- rangeMin/gcMax

            gc <- gc*scale*expand
            gc$id <- row.names(gc)
            gcNames <- names(gc)
            gcTrans <- c(gcNames[2], gcNames[1])
           
            if('lower' %in% where){
                 segs <- list(geom_segment(data=gc, aes_string(x=0, y=0, xend=gcNames[1], yend=gcNames[2], col=NULL, shape=NULL), arrow=grid::arrow(length=unit(0.2,"cm")), color="red", ...), geom_text(data=gc, aes_string(x=gcNames[1], y=gcNames[2], col=NULL, shape=NULL, label="id"), col='black', size=2.5, ...))
                ggpairsObj <- putPlot(ggpairsObj, gp + segs, d2, d1)
                 if(debug){
                     cat('row=',d2, 'col=', d1)
                     print(gp + segs)
                 }
            }
            if('upper' %in% where){
                 segs <- list(geom_segment(data=gc, aes_string(x=0, y=0, xend=gcNames[2], yend=gcNames[1], col=NULL, shape=NULL), arrow=grid::arrow(length=unit(0.2,"cm")), color="red"), geom_text(data=gc, aes_string(x=gcNames[2], y=gcNames[1], col=NULL, shape=NULL, label="id"), col='black', size=2.5))
                gp <- getPlot(ggpairsObj, d1, d2)
                ggpairsObj <- putPlot(ggpairsObj, gp + segs, d1, d2)
                 if(debug){
                     cat('row=',d1, 'col=', d2)
                     print(gp + segs)
                 }

            }
        }
    }
    ggpairsObj
}

##' @importFrom scales hue_pal
twoPal <- hue_pal()(2)

annotateUniplot <- function(plot, lc, genesToShow=5, expand.x=1, expand.y=1, negColor=twoPal[2], posColor=twoPal[1], ...){
    panel <- ggplot_build(plot)$panel
    rangeMin <- min(abs(panel$ranges[[1]]$x.range))
    ytop <- max(panel$ranges[[1]]$y.range)
    gc <- getCoord(lc, 1, 2, genesToShow)
    gcMax <- max(abs(gc))
    scale <- rangeMin/gcMax

    gc <- gc*scale*expand.x
    gc$id <- row.names(gc)
    gcNames <- names(gc)
    gc$y <- seq(ytop*.9, to=ytop*.7/expand.y, length.out=nrow(gc))
    gc$col <- ifelse(gc[,1]>0, posColor, negColor)
   
    segs <- list(geom_segment(data=gc, aes_string(x=0, y='y', xend=gcNames[1], yend='y', color=NULL, shape=NULL), arrow=grid::arrow(length=unit(0.2,"cm")), color=gc$col,...), geom_text(data=gc, aes_string(x=gcNames[1], y='y', col=NULL, shape=NULL, label="id"), col='black', size=2.5, ..., vjust=0))
    plot + segs
}
