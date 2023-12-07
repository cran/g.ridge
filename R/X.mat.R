
#' @title X.mat (generating a design matrix)
#' @description A design matrix (X; nrow(X)=n, ncol(X)=p) is generated
#'  by random numbers as previously used in our simulation studies
#'  (Section 5 of Yang and Emura (2017); p.6093).
#'  The design matrix has two blocks of correlated regressors
#'   (Pearson correlation=0.5): the first q regressors and the second r regressors.
#'   Other p-q-r regressors are independent.
#'   If regressors are gene expressions, the correlated blocks may be regarded as
#'    "gene pathways" (Emura et al. 2012).
#' @param n the number of rows (samples)
#' @param p the number of columns (regressors)
#' @param q the number of correlated regressors in the first block (1<=q<p, q+r<p)
#' @param r the number of correlated regressors in the second block (1<=r<p, q+r<p)
#' @return a matrix X (nrow(X)=n, ncol(X)=p)
#' @examples
#' X.mat(n=10,p=5,q=2,r=2)
#' X.mat(n=100,p=50,q=10,r=10) # Case I in Section 5 of Yang and Emura (2017)
#' @references Yang SP, Emura T (2017) A Bayesian approach with generalized ridge estimation
#'  for high-dimensional regression and testing, Commun Stat-Simul 46(8): 6083-105
#' @references Emura T, Chen YH, Chen HY (2012) Survival prediction based on compound covariate method
#'   under Cox proportional hazard models PLoS ONE 7(10) doi:10.1371/journal.pone.0047627

X.mat=function(n,p,q,r){
  X=matrix(c(rnorm(n*(q+r),0,sqrt(2)^-1),rnorm(n*(p-(q+r)))),n,p)
  u.q=rnorm(n,0,sqrt(2)^-1)
  v.r=rnorm(n,0,sqrt(2)^-1)
  X[,1:q]=X[,1:q]+u.q
  X[,(q+1):(q+r)]=X[,(q+1):(q+r)]+v.r
  return(X)
}
