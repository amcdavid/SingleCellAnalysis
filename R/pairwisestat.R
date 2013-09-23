gencor <- function(ct, maxmissingness = .99){
  #trim names
  colnames(ct) = substr(colnames(ct), 5, 20)
  #remove genes that were observed with lower frequency than 1-missingness
  numobs = apply(ct>0, 2, sum, na.rm=T)
  nullgenes = numobs < nrow(ct)*(1-maxmissingness)
  ct2 = ct[,!nullgenes]
  corbygene = cor(ct2, use="pairwise.complete.obs")
  corbygene = submatrix(corbygene) #remove genes without significant correlations
  heatmap.2(corbygene, trace="none", keysize=1)
  invisible(corbygene)
}

rmse.sim <- function(N=100){
p = c(.005, .01, .1)
out = expand.grid(p1=p, p2=p)
out = cbind(out,matrix(0, nrow(out), ncol=N))
reppp = vector()
for(i in 1:nrow(out)){
  for(n in 1:N){
    x<- rbinom(500, 1, p=out[i,1])
    y<-rbinom(500, 1, p=out[i,2])
    reppp[n] = rmse.test(x,y, P=5e3)$p
  }
out[i,3:(N+2)] =reppp
}
out
}


cap = function(covarmats, Nsamp){
  #Bartlet's modified likelihood ratio test
  p = nrow(covarmats[[1]])
  stopifnot(p==ncol(covarmats[[1]]))
  stopifnot(length(covarmats)==length(Nsamp))
  k = length(covarmats)
  covsum = do.call(sum, covarmats)
  Nsum = sum(Nsamp)
  dets = unlist(lapply(det, covarmats))
  detsum = det(covsum)
  prod(abs(dets)^((Nsamp-1)/2))*(Nsum - k)^(p*(Nsum-k)/2)/
    (abs(det(covsum))^((Nsum-k)/2)*prod((Nsamp-1)^(p*(Nsamp-1))) )
}

submatrix = function(mat, rowcrit=function(arow){sort(abs(arow), decreasing=T)[2]>.2 }){
  val = apply(mat, 2, rowcrit)
  mat1 = mat[val,]
  mat1[,val]
}

rmse.test <- function(x, y, P=5e4, nulldist=F){
  stopifnot(length(x)==length(y))
  N = length(x)
  v0 = matrix(c(sum(x*y), sum(x*(1-y)), sum((1-x)*y), sum((1-x)*(1-y))), nrow=4)
  out = .rmse(v0, N, expected=T)
  expHat = out[[2]]
  obs = out[[1]]
  v = rmultinom(P, N, expHat)
  h0 = .rmse(v, N)
  p = mean(obs <= h0)
  if(nulldist){
    return(list(p.value=p, obs=obs, null=h0))
  }
  return(list(p.value=p, obs=obs))
}

.rmse <- function(v, N, expected=F){
  #v is a 4xN matrix of simulated obs
  A = (v[1,]+v[2,])/N
  B = (v[1,]+v[3,])/N
  exp = rbind(A*B, A*(1-B), (1-A)*B, (1-A)*(1-B))*N
  out = sqrt(apply((v - exp)^2, 2, sum))
  if(!expected){
    return(out)
  }
    
  return(list(out, exp))
}

rmsepairws <- function(scr3, assays){
  pval = matrix(0, ncol=length(assays), nrow=length(assays), dimnames=list(assays, assays))
  oval = pval
  for(i in 1:(length(assays)-1)){
    for(j in (i+1):length(assays)){
      test = rmse.test(scr3[,assays[i]], scr3[,assays[j]])
      pval[i,j] = test$p
      oval[i,j] = test$obs
    }
  }
  
  pval = t(pval) + pval
  oval = t(oval) + oval
  
  return(list(pval, oval))
}

##' Chi-square tests of independence for each pair of genes
##'
##' Conduct a chi-square tests of independence for each pair of genes.
##' The binary (expressed/not) values are used to build the contingency tables.
##' The chi-square test statistics for each pair of genes is returned, with zeroes placed along the main diagonal.
##' Significance could be assessed (in large samples) by comparison to a chi-square distribution with 1-dof.
##' @param sc SingleCellAssay, on a thresholded (et) layer
##' @return matrix of chisq test statistics
##' @author andrew
##' @export
chisq.pairwise <- function (sc) {
  if(layername(sc) != 'et') warning('Tests fail unless on a thresholded layer')
  scr3 <- exprs(sc)==0
  marg = apply(scr3, 2, mean)
  onames = names(marg)
  #onames = substr(onames, 5, 30)
  exp=array(dim=c(length(marg),length(marg), 4), dimnames=list(onames, onames, NULL))
  exp[,,1] = outer(marg, marg)*nrow(scr3)
  exp[,,2] = outer(1-marg, marg)*nrow(scr3)
  exp[,,3] = outer(marg, 1-marg)*nrow(scr3)
  exp[,,4] = outer(1-marg, 1-marg)*nrow(scr3)
  
  
  xy=exp
  xy[,,1] = t(as.matrix(scr3)) %*% as.matrix(scr3)
  xy[,,2] = t(as.matrix(!scr3)) %*% as.matrix(scr3)
  xy[,,3] = t(as.matrix(scr3)) %*% as.matrix(!scr3)
  xy[,,4] = t(as.matrix(!scr3)) %*% as.matrix(!scr3)
  
  chi = (xy-exp - .5*sign(xy-exp))^2/(exp+.0001) #continuity correction, and avoid division by 0
  chi2 = apply(chi, c(1, 2), sum)
  diag(chi2) = 0 #of course X_i is not independent of X_i!
  #smallexp = apply(exp<5, c(1, 2), any)
  #chi2[smallexp] = NA
  chi2
}


chisq.stim <- function (scr3, stim) {
  ng = ncol(scr3)
  gs = vector(length=ng)
  tt = list(exp=0)
  for(g in 1:ng){
    tab = table(scr3[,g], stim)
    if(all(dim(tab)==c(2, 2))){
      tt = fisher.test(scr3[,g], stim) 
      gs[g] = tt$p  
    }
     else{
       gs[g] = 1
     }
    }
  return(gs)
}

