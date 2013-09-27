heat2 <- function(mat, Rowv=NA, Colv=NA, symbreaks=TRUE, col=if(symbreaks) redblue else redgreen, scale='none', ...){
  mat <- as.matrix(mat)
  library(gplots, pos=length(search()))
  heatmap.2(mat, trace='none', scale='none', Rowv=Rowv, Colv=Colv, symbreaks=symbreaks, col=col, ...)
}

panel.pos <- function(x, y){
  sel <- x>0 & y>0 & !is.na(x) & !is.na(y)
 x.pos <- x[sel]
 y.pos <- y[sel]
 par(new=TRUE)
 plot(x.pos, y.pos, col='gray', cex=.4)
 lines(lowess(x.pos, y.pos))
}
panel.discrete <- function(x, y){
 tt <- table(x>0, y>0, deparse.level=2)
  par(new=TRUE)
 if(any(dim(tt)==1)) plot(0, 0, type='n')
 else {
   #dimnames(tt) <- list(c('-', '+'))[rep(1, length(dimnames(tt)))]
   plot(prop.table(tt), main='', cex=2)
 }
}

panel.mixed <- function(x, y){
  par(new=TRUE)
  y.pos <- y[y>0]
  x.sub <- x[y>0]
  boxplot(y.pos ~ x.sub>0)
}

uniParams <- function(sc.sub){
  ee.pos <- exprs(sc.sub)
  ee.pos[ee.pos==0] <- NA
  m <- condmean(sc.sub)
  pi <- freq(sc.sub)
  s <- apply(ee.pos, 2, sd, na.rm=TRUE)
  nu <- exp(m+.5*s^2)
  tau <- nu*sqrt(exp(s^2)-1)
  mu <- nu*pi
  sigma <- sqrt(pi*(nu^2*(1-pi)+tau^2))
  kappa <- sqrt(log(1+tau^2/nu^2))
  data.frame(m=m, s=s, pi=pi, nu=nu, tau=tau, mu=mu, sigma=sigma, kappa=kappa)
}



plotScDiscrete <- function(sc.sub, d){
  ee.sc <- exprs(sc.sub)
  suppressWarnings(cor.sc <- cor0(ee.sc>0))
  if(!missing(d) && !is.null(d)){
    heat2(cor.sc, Rowv=d$rowDendrogram, Colv=d$colDendrogram, breaks=d$breaks)
  } else{
    d <- heat2(cor.sc, Rowv=TRUE, Colv=TRUE)
  }
  d
}

panel.xyreg <- function(...){panel.xyplot(...); panel.smoother(..., method='rlm')}
panel.r2 <- function(...){
  args <- list(...)
  x <- args$x; y <- args$y
  d <- args$heatmapargs
  r <- cor(x, y, use='pairwise.complete.obs')
  if(!is.null(d)){
  r.col <- d$col[findInterval(r, d$breaks, all.inside=TRUE)]
  panel.fill(col=r.col)
}
  panel.text(x=mean(range(x, na.rm=TRUE)), y=mean(range(y, na.rm=TRUE)), labels=round(r, 2))
}

corByStrata <- function(sc.sub, groups, dendro, ...){
patxstim <- unique(cData(sc.sub)[,groups, drop=FALSE])
sc.bypat <- split(sc.sub, groups)
suppressWarnings(cor.bypat <- lapply(sc.bypat, function(x){
  ee <- exprs(x)
  n <- nrow(ee)
 corout <- cor0(ee, ...)
}))

cc <- opts_chunk$get()
opts_chunk$set(fig.width=6, fig.height=6)
  if(missing(dendro) || is.null(dendro)){
    dendro <-  heat2(cor.bypat[[1]], Rowv=TRUE, Colv=TRUE, symm=TRUE, sub=paste(patxstim[1,], collapse=':'))
 } else{
     heat2(cor.bypat[[1]], Rowv=dendro$rowDendrogram, Colv=dendro$colDendrogram, symm=TRUE, breaks=dendro$breaks, sub=paste(patxstim[1,], collapse=':'))
   }


opts_chunk$set(cc)
for(i in seq(2, nrow(patxstim))){
  ptid <- patxstim[i,]
  x <- cor.bypat[[i]]
 heat2(x, Rowv=dendro$rowDendrogram, Colv=dendro$colDendrogram, symm=TRUE, breaks=dendro$breaks, key=FALSE, dendrogram='none', labRow=NA, labCol=NA, sub=paste(patxstim[i,], collapse=':'))
}
invisible(list(listofcor=cor.bypat, dendro=dendro))
}

randomForestClass <-  function(sc, groups){
library(randomForest)
class <- factor(do.call(paste, cData(sc)[,groups, drop=FALSE]))
rpartcont <- data.frame(class, exprs(sc))
rf <- randomForest(class ~., data=rpartcont,type='class', proximity=TRUE)
print(rf)
rf
}

predictRFOnSCA <- function(rf, sc, groups){
if(!missing(sc)){
rpartcont <- data.frame(exprs(sc))
class <- factor(do.call(paste, cData(sc)[,groups, drop=FALSE]), levels=rf$classes)
tab <- table(true=class, predicted=predict(rf, rpartcont))
} else{
  tab <- rf$confusion[,-ncol(rf$confusion)]
}
err <- 1-sum(diag(tab))/sum(tab)
list(tab=tab, err=err)
}

randomizeMatrix <- function(mat, grouping){
for(i in seq_along(unique(grouping))){
   for(j in seq_len(ncol(mat))){
     subs <- mat[grouping %in% unique(grouping)[i],j]
     mat[grouping %in% unique(grouping)[i],j] <- sample(subs)
   }
 }
mat
}

randomizeStrata <- function(sca, groups){
 pat <- cData(sca)[[groups]]
 ee <- exprs(sca)
 randomizeMatrix(ee, pat)
 scNew <- copy(sca)
 exprs(scNew) <- ee
 scNew
}


emNorm <- function(sca){
require(norm)
ee.na <- exprs(sca)
ee.na[ee.na==0] <- NA
pp <- SingleCellAnalysis::prelim.norm(ee.na)
mu0 <- condmean(sca)
cov0 <- cov(ee.na, use='pairwise.complete.obs')
lambda.min <- abs(min(eigen(cov0, only.values=TRUE)$values))+.4
cov0 <- cov0+diag(nrow(cov0))*lambda.min
emout <- em.norm(pp, makeparam.norm(pp, list(mu0, cov0)), criterion=1e-4, maxits=3000)
r <- getparam.norm(pp, emout, cor=TRUE)$r
dimnames(r) <- list(dimnames(ee.na)[[2]], dimnames(ee.na)[[2]])
r
}

compareCorRank <- function(cor1, cor2, alpha=TRUE, ncomp=100){
  diag(cor1) <- 0
  diag(cor2) <- 0
  m1 <- melt(cor1)[melt(upper.tri(cor1))$value,]
  m2 <- melt(cor2)[melt(upper.tri(cor1))$value,]
  m1 <- m1[order(-abs(m1$value)), ]
  m2 <- m2[order(-abs(m2$value)), ]
  m1 <- cbind(m1, idx=1:nrow(m1))
  m2 <- cbind(m2, idx=1:nrow(m2))
  mm <- merge(m1, m2, by=c('X1', 'X2'))
  mm <- mm[order(-abs(mm$value.x)), ]
  print(parallelplot(~mm[1:ncomp,c('idx.x', 'idx.y')], alpha=if(alpha) abs(mm$value.y[1:ncomp])^2 else 1, common.scale=TRUE))
  xyplot(idx.x ~ idx.y, alpha=1, mm)#mm$value.y^2, mm)
}

pairwiseCoexp <- function(ee.subset){
pairs(ee.subset, upper.panel=panel.mixed, lower.panel=panel.mixed)
pairs(ee.subset, upper.panel=panel.pos, lower.panel=panel.discrete)
}

## mat: exprs matrix (dichotomized or otherwise)
## test: NULL return raw correlations
##       permute truncate to zero if permuted pvalue is less than .1
##       confint return min(lowerint, zero) of bootstrapped conf intervals
cor0 <- function(mat, test=NULL, confint=.025, pvalue=.1){
  r <- suppressWarnings(cor(mat, use='pairwise.complete.obs'))
  if(!is.null(test)){
    test <- match.arg(test, c('permute', 'confint'))
    if(test == 'confint'){
     stopifnot(is.numeric(confint))
    rb <- replicate(50, suppressWarnings(cor(mat[sample(nrow(mat), replace=TRUE),], use='pairwise.complete.obs')))
    ub <- apply(rb, c(1, 2), quantile, na.rm=TRUE, p=1-confint)
    lb <- apply(rb, c(1, 2), quantile, na.rm=TRUE, p=confint)
    r <- apply(rb, c(1, 2), mean, na.rm=TRUE)
    return(list(lower=lb, mean=r, upper=ub))
  } else if(test=='permute'){
    stopifnot(is.numeric(pvalue))
    rb <- replicate(100, suppressWarnings(cor(randomizeMatrix(mat, rep(1, nrow(mat))), use='pairwise.complete.obs')))
    crit <- apply(abs(rb), c(1, 2), quantile, na.rm=TRUE, p=1-pvalue)
    r[abs(r)<crit] <- 0
  }
  }
diag(r) <- 0
r[is.na(r)] <- 0
r
}

compareCor <- function(listofcor){
upper <- lapply(listofcor, function(x){x[lower.tri(x)] <- NA; x})
m <- melt(upper)
m <- subset(m, !is.na(value))
vectorized <- cast(m, X1 + X2 ~L1)
dropidx <- 1:2
mat <- as.matrix(vectorized[, -dropidx])
dimnames(mat)[[2]] <- names(vectorized)[-dropidx]
heat2(mat, Rowv=TRUE, Colv=TRUE, scale='row')
}

compareCorScatter <- function(listoftwocor){
  cor1 <- listoftwocor[[1]]
  cor2 <- listoftwocor[[2]]
   diag(cor1) <- 0
  diag(cor2) <- 0
  m1 <- melt(cor1)[melt(upper.tri(cor1))$value,]
  m2 <- melt(cor2)[melt(upper.tri(cor1))$value,]
  mm <- rbind(cbind(id=names(listoftwocor)[1], m1), cbind(id=names(listoftwocor)[2], m2))
  
}



heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }
 
    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors) && !is.null(ColSideColors)) {
           #if (!is.matrix(ColSideColors))
           #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei <- c(lhei[1], 0.2, lhei[2])
             lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], 0.2, lwid[2])
            lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
            }
        }
    }
 
    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
            }
        }
    }
 
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
 
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

##' Plot a heatmap of a SingleCellAssay expression matrix
##'
##' Plot a heatmap of a SingleCellAssay expression matrix (eg, from calling exprs).
##' Optionally cluster and label rows by rowFactor, columns by colFactor, and/or sort the expression matrix by rowFactor/colFactor.
##' @param ee.sc expression matrix 
##' @param rowFactor named list of character vectors to be used as a labelling color schema for the rows.
##' @param cluster boolean. Apply hclust to rows and columns.
##' @param order.by.factor Sort expression matrix by rowFactors and colFactors.  If cluster is TRUE than cluster within rowFactor/colFactor.
##' @param colFactor named list of character vectors to be used as a labelling color schema for the columns.
##' @param ... passed to heatmap.2
##' @return prints plot to current graphic device
##' @export
singleCellHeat <- function(ee.sc, rowFactor, cluster=TRUE, order.by.factor=!cluster, colFactor, ...){
  if(cluster){
rowD <- as.dendrogram(hclust(dist(ee.sc)))
colD <- as.dendrogram(hclust(dist(t(ee.sc))))
}
ee.na <- ee.sc
ee.na[ee.na==0] <- NA
row.names(ee.na) <- NULL
if(!missing(rowFactor)){
  if(!is.list(rowFactor)) rowFactor <- list(rowFactor)
  rowFactor <- lapply(rowFactor, function(rf) as.factor(rf)[TRUE, drop=TRUE])
  RowSideColors <- as.matrix(do.call(rbind, lapply(rowFactor, function(rf){nfact <- length(levels(rf)); rainbow(nfact, s=.7)[as.numeric(rf)] })))
  rownames(RowSideColors) <- names(rowFactor)
} else{
RowSideColors <- rep('clear', nrow(ee.sc))
}
  
  if(!missing(colFactor)){
  if(!is.list(colFactor)) colFactor <- list(colFactor)
  colFactor <- lapply(colFactor, function(rf) as.factor(rf)[TRUE, drop=TRUE])
  ColSideColors <- as.matrix(do.call(cbind, lapply(colFactor, function(rf){nfact <- length(levels(rf)); rainbow(nfact, s=.7)[as.numeric(rf)] })))
  colnames(ColSideColors) <- names(colFactor)
} else{
  ColSideColors=matrix('black', nrow=ncol(ee.sc), ncol=1)
}

  
  if(order.by.factor){
    rf.inter <- as.ordered(interaction(rowFactor, drop=TRUE))
    ee.sc <- ee.sc[order(as.numeric(rf.inter)),]
    RowSideColors <- RowSideColors[,order(rf.inter), drop=FALSE]
    rf.inter <- sort(rf.inter)
    rn <- as.character(rf.inter)
    rn <- abbreviate(rn)
    rn[duplicated(rn)] <- ''
    row.names(ee.sc) <- rn

    if(!missing(colFactor)){
    col.inter <- as.ordered(interaction(colFactor, drop=TRUE))
    ee.sc <- ee.sc[,order(as.numeric(col.inter))]
    ColSideColors <- ColSideColors[order(col.inter),, drop=FALSE]
    col.inter <- sort(col.inter)
  } 
    if(cluster){
    ee.reorder <- lapply(levels(rf.inter), function(rf){
      ee.sub <- ee.sc[rf.inter==rf,]
      d <- as.dendrogram(hclust(dist(ee.sub)))
      ee.sub[order.dendrogram(d),]
    })
    ee.sc <- do.call(rbind, ee.reorder)
    rowD <- FALSE

    if(!missing(colFactor)){
     ee.reorder <- lapply(levels(col.inter), function(rf){
      ee.sub <- ee.sc[,col.inter==rf]
      d <- as.dendrogram(hclust(dist(t(ee.sub))))
      ee.sub[,order.dendrogram(d)]
    })
    ee.sc <- do.call(cbind, ee.reorder)
    colD <- FALSE
   }
  }
  }

  suppressWarnings(
  if(cluster){
heatmap.3(ee.sc, Colv=colD, col=topo.colors, trace='none', RowSideColors=RowSideColors, Rowv=rowD, ColSideColors=ColSideColors, ...)
} else{
heatmap.3(ee.sc, col=topo.colors, trace='none', RowSideColors=RowSideColors, ColSideColors=ColSideColors, ...)
}
    )
}
