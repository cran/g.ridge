
#' @title g.ridge (generalized ridge regression)
#' @description Generalized ridge regression with the optimal shrinkage parameter.
#'  Ridge regression (Hoerl and Kennard, 1970) and
#'  generalized ridge regression (Yang and Emura 2017) are implemented.
#'  Tuning parameters are optimized by minimizing the CGV function (by the function CGV(.)):
#'   See Golub et al. (1979), and Sections 2.3 and 3.3 of Yang and Emura (2017).
#' @param X design matrix of explanatory variables (regressors)
#' @param Y vector of response variables
#' @param method "HK" or "YE" for Hoerl and Kennard (1970) or Yang and Emura (2017)
#' @param kmax maximum possible value for the shrinkage parameter (the "lambda" parameter),
#'  where the parameter is optimized in the interval (0, kmax).
#' @return lambda: optimized shrinkage parameter
#' @return delta: the optimized thresholding parameter
#' @return estimate: regression coefficients (beta)
#' @return SE: Standard Error
#' @return Z: Z-value for testing beta=0
#' @return SE: P-value for testing beta=0
#' @return Sigma: variance estimate of the error distribution
#'  (the square of the standard deviation)
#' @return delta: thresholding parameter
#' @examples
#' n=100 # no. of observations
#' p=100 # no. of dimensions
#' q=r=10 # no. of nonzero coefficients
#' beta=c(rep(0.5,q),rep(0.5,r),rep(0,p-q-r))
#' X=X.mat(n,p,q,r)
#' Y=X%*%beta+rnorm(n,0,1)
#' g.ridge(X,Y-mean(Y),method="HK",kmax=200)
#' g.ridge(X,Y-mean(Y),method="YE",kmax=200)
#' @references Yang SP, Emura T (2017) A Bayesian approach with generalized ridge estimation
#'  for high-dimensional regression and testing, Commun Stat-Simul 46(8): 6083-105.
#' @references Hoerl AE, Kennard RW (1970) Ridge regression: Biased estimation for
#'  nonorthogonal problems. Technometrics 12:55–67.
#'
g.ridge=function(X,Y,method="HK",kmax=500){

X=as.matrix(X)
delta.hat=NULL
n=nrow(X)
p=ncol(X)
XtX=t(X)%*%X
XtY=t(X)%*%Y
varname=colnames(X)

if(method=="HK"){# Hoerl and Kennard (1970)
  gcv=function(k){ GCV(X,Y,k) }
  k.hat=optimize(gcv,interval=c(0,kmax))$minimum
  beta.hat=solve(XtX+k.hat*diag(p))%*%XtY
  gcv.vec=Vectorize(gcv)
  ylab=expression("V( "*lambda*" )")
  V=solve(XtX+k.hat*diag(p))
}

if(method=="YE"){# Yang and Emura (2017)
  ini.est=as.vector(diag(diag(XtX)^-1)%*%XtY)
  min=arg.min=numeric(101)
  Delta=seq(0,3,length=101)
  for(j in 1:101){
    delta = Delta[j]
    W=diag(((abs(ini.est/sd(ini.est))>delta)+1)^-1)
    gcv=function(k){ GCV(X,Y,k,W) }
    temp = optimize(gcv,interval = c(0,kmax))
    arg.min[j] = temp$minimum
    min[j] = temp$objective
  }
  k.hat=arg.min[which.min(min)]
  delta.hat= Delta[which.min(min)]
  W=diag(((abs(ini.est/sd(ini.est))>delta.hat)+1)^-1)
  beta.hat=solve(XtX+k.hat*W)%*%XtY
  W.hat=diag(((abs(ini.est/sd(ini.est))>delta.hat)+1)^-1)
  gcv=function(k){ GCV(X,Y,k,W.hat) }
  gcv.vec=Vectorize(gcv)
  ylab=expression("V( "*lambda*","*hat(Delta)*" )")
  V=solve(XtX+k.hat*W.hat)
}

A=X%*%V%*%t(X)
DF=n-sum(diag(2*A-A%*%A))
sigma.hat=sum((Y-X%*%beta.hat)^2)/DF
Cov=sigma.hat*t(V)%*%XtX%*%V
SE=sqrt(diag(Cov))
Z=beta.hat/SE
P=2*(1-pnorm(abs(Z)))

oldpar=par(no.readonly=TRUE)
on.exit(par(oldpar))

par(mar=c(5,5,3,2))
curve(gcv.vec,from=kmax/100,to=kmax,
      xlab=expression(lambda),ylab=ylab)
points(k.hat,gcv(k.hat),col="red")
mtext(expression(hat(lambda)),at=k.hat,side=1,col="red")

Beta=data.frame(estimate=beta.hat,SE=SE,Z=Z,P=P)

list(lambda=k.hat,delta=delta.hat,
     beta=Beta,sigma=sigma.hat)
}
