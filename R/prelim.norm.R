prelim.norm <- 
function (x) 
{
    if (is.vector(x)) 
        x <- matrix(x, length(x), 1)
    n <- nrow(x)
    p <- ncol(x)
    storage.mode(x) <- "double"
    r <- 1 * is.na(x)
    nmis <- as.integer(apply(r, 2, sum))
    names(nmis) <- dimnames(x)[[2]]
    mdp <- (r %*% (2^((1:ncol(x)) - 1))) + 1
    ro <- order(mdp)
    x <- matrix(x[ro, ], n, p)
    mdp <- mdp[ro]
    r <- matrix(r[ro, ], n, p)
    ro <- order(ro)
    mdpst <- as.integer(seq(along = mdp)[!duplicated(mdp)])
    mdp <- unique(mdp)
    npatt <- length(mdpst)
    r <- 1 - r
    r <- matrix(r[mdpst, ], npatt, p)
    if (npatt == 1) 
        tmp <- format(n)
    if (npatt > 1) 
        tmp <- format(c(mdpst[2:npatt], n + 1) - mdpst)
    dimnames(r) <- list(tmp, dimnames(x)[[2]])
    storage.mode(r) <- "integer"
    if (sum(is.na(x)) < length(x)) {
        mvcode <- as.double(max(x[!is.na(x)]) + 1000)
        x <- .na.to.snglcode(x, mvcode)
        tmp <- .Fortran("ctrsc", x, n, p, numeric(p), numeric(p), 
            mvcode, PACKAGE = "norm")
        x <- tmp[[1]]
        xbar <- tmp[[4]]
        sdv <- tmp[[5]]
        x <- .code.to.na(x, mvcode)
    }
    if (sum(is.na(x)) == length(x)) {
        xbar <- rep(0, p)
        sdv <- rep(1, p)
    }
    d <- as.integer((2 + 3 * p + p^2)/2)
    psi <- .Fortran("mkpsi", p, matrix(as.integer(0), p + 1, 
        p + 1), PACKAGE = "norm")[[2]]
    if (npatt > 1) 
        nmdp <- as.integer(c(mdpst[-1], n + 1) - mdpst)
    if (npatt == 1) 
        nmdp <- n
    sj <- .Fortran("sjn", p, npatt, r, integer(p), PACKAGE = "norm")[[4]]
    nmon <- .Fortran("nmons", p, npatt, r, nmdp, sj, integer(p), 
        PACKAGE = "norm")[[6]]
    last <- .Fortran("lasts", p, npatt, sj, integer(npatt), PACKAGE = "norm")[[4]]
    tmp <- .Fortran("layers", p, sj, integer(p), integer(1), 
        PACKAGE = "norm")
    layer <- tmp[[3]]
    nlayer <- tmp[[4]]
    list(x = x, n = n, p = p, r = r, nmis = nmis, ro = ro, mdpst = mdpst, 
        nmdp = nmdp, npatt = npatt, xbar = xbar, sdv = sdv, d = d, 
        psi = psi, sj = sj, nmon = nmon, last = last, layer = layer, 
        nlayer = nlayer)
}
