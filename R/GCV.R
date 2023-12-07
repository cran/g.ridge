
#' @title GCV (generalized cross-validation)
#' @description The CGV function gives the sum of cross-validated squared errors
#' that can be used to optimize tuning parameters in ridge regression and
#'  generalized ridge regression.
#'   See Golub et al. (1979), and Sections 2.3 and 3.3 of Yang and Emura (2017) for details.
#' @param X matrix of explanatory variables (design matrix)
#' @param Y vector of response variables
#' @param k shrinkage parameter (>0); it is the "lambda" parameter
#' @param W matrix of weights (default is the identity matrix)
#' @return The value of GCV
#' @examples
#' n=100 # no. of observations
#' p=100 # no. of dimensions
#' q=r=10 # no. of nonzero coefficients
#' beta=c(rep(0.5,q),rep(0.5,r),rep(0,p-q-r))
#' X=X.mat(n,p,q,r)
#' Y=X%*%beta+rnorm(n,0,1)
#' GCV(X,Y,k=1)
#' @references Yang SP, Emura T (2017) A Bayesian approach with generalized ridge estimation
#'  for high-dimensional regression and testing, Commun Stat-Simul 46(8): 6083-105.
#' @references Golub GH, Heath M, Wahba G (1979) Generalized cross-validation as
#'  a method for choosing a good ridge parameter. Technometrics 21:215–223.

GCV=function(X,Y,k,W=diag(ncol(X))){
  X=as.matrix(X)
  n=nrow(X)
  A=X%*%solve(t(X)%*%X+k*W)%*%t(X)
  I_A=diag(n)-A
  temp=I_A%*%as.matrix(Y)
  D=(sum(diag(I_A))/n)^2
  N=t(temp)%*%temp/n
  return(N/D)
}

