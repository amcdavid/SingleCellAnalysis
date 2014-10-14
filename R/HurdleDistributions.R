expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

.checkArgs <- function(x, G, H, K){
    stopifnot(is.matrix(H))
    stopifnot(ncol(H)==nrow(H))
    stopifnot(dim(G)==dim(H))
    stopifnot(dim(K)==dim(H))
    stopifnot(ncol(x)==ncol(H))
}

## joint density function
## x is row vector!!!
## But should be a column vector...
dHurdle210 <- function(x, G, H, K, tol=5e-2){
    if(length(dim(x))<2) x <- t(as.matrix(x))
    .checkArgs(x, G, H, K)
    xI <- (abs(x)>tol)*1
    logDens <- xI %*% G %*% t(xI) + (xI*x)%*%H%*%t(xI) - .5*(xI*x)%*%K%*%t(xI*x)
    exp(logDens)
}

## simulate from the conditional distribution of a hurdle model
## where i gives the indices being provided in x
rCondHurdle210 <- function(x, i, G, H, K, tol=5e-2){
    stopifnot(length(i) == ncol(x), length(i) == ncol(G)-1)
    .checkArgs(cbind(1, x), G, H, K)
    xI <- (abs(x)>tol)*1
    midx <- seq_len(ncol(G))
    noti <- setdiff(midx, i)
    
    Gba <- G[noti,noti,drop=TRUE]+2*xI%*%G[i,noti,drop=FALSE] + (xI*x)%*%t(H[noti, i,drop=FALSE])
    Hba <- H[noti, noti]+xI%*%H[i,noti,drop=FALSE] - (xI*x)%*%K[i,noti,drop=FALSE]
    Kba <- K[noti, noti,drop=TRUE]

    logitP <- Gba-.5*log(Kba/(2*pi))+Hba^2/(2*Kba)
    mu <- Hba/Kba
    yI <- runif(nrow(x))<expit(logitP)
    y <- rnorm(nrow(x), sd=1/sqrt(Kba))+mu
    y*yI
}

rv <- function(x){
    dim(x) <- c(1, length(x))
    x
}

rGibbsHurdle <- function(G, H, K, Nt, burnin=floor(Nt/2), tol=5e-2){
    p <- ncol(G)
    .checkArgs(matrix(ncol=p), G, H, K)
    ## samp is internally arrayed backwards from what rCondHurdle uses,
    ## but it is much easier to have it in column-major order for the purposes of gibbs sampling
    samp <- rep(NA, p*Nt)
    samp[seq_len(p)] <-0
    ## pointer to previous p-1 positions in vectorized sample
    notp <- seq_len(p-1)+1
    for(i in seq(p+1, Nt*p)){
            ## and coordinates of these samples
            notpidx <- ((notp-1) %% p)+1
            samp[i] <- rCondHurdle210(rv(samp[notp]), i=notpidx, G=G, H=H, K=K, tol=tol)
            notp <- notp+1
        }
    ## we are in column-major order, but will want to be in row-major for consistency elsewhere
    t(matrix(samp, nrow=p))[-seq_len(burnin),]
}
