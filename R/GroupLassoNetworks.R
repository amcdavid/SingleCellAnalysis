## Code to be used with grplasso
HurdleReg <- function(){
    grpl.model(
        invlink = function(eta, cont){
            mu.1 <- LogReg()@invlink(eta[!cont])
            mu.2 <- LinReg()@invlink(eta[cont])
            c(mu.1, mu.2)
        },
        link = function(mu, cont){
            eta.1 <- LogReg()@link(mu[!cont])
            eta.2 <- LinReg()@invlink(mu[cont])
            c(eta.1, eta.2)
        },
        nloglik= function(y, eta, weights, cont){
            nloglik.1 <- LogReg()@nloglik(y[!cont], eta[!cont], weights[!cont])
            nloglik.2 <- LinReg()@nloglik(y[cont], eta[cont], weights[cont])
            nloglik.1 + nloglik.2
        },
        ngradient=function(x, y, mu, weights, cont){
             is.cont <- substr(colnames(x), 1, 4)=='cont'
            x.1 <- x[!cont,!is.cont, drop=FALSE]
            x.2 <- x[cont,is.cont, drop=FALSE]
            grad.1 <- grad.2 <- matrix(NA, nrow=0, ncol=0)
            if(ncol(x.1)>0)
                grad.1 <- LogReg()@ngradient(x.1, y[!cont], mu[!cont], weights[!cont])
             if(ncol(x.2)>0)
                 grad.2 <- LinReg()@ngradient(x.2, y[cont], mu[cont], weights[cont])
            cbind(grad.1, grad.2)
        },
        nhessian=function(x, mu, weights, cont){
            is.cont <- substr(colnames(x), 1, 4)=='cont'
            x.1 <- x[!cont,!is.cont, drop=FALSE]
            x.2 <- x[cont,is.cont, drop=FALSE]
            h.1 <- h.2 <- matrix(NA, nrow=0, ncol=0)
            if(ncol(x.1)>0)
                h.1 <- LogReg()@nhessian(x.1, mu[!cont], weights[!cont])
            if(ncol(x.2)>0)
                h.2 <- LinReg()@nhessian(x.2, mu[cont], weights[cont])
            as.matrix(bdiag(h.1, h.2))
        },
        check=function(y) TRUE)
}

adjustLambda <- function(nobs.d, nobs.c, sigma.c){
    nobs.c/(nobs.d/sqrt(sigma[,2]))

}
