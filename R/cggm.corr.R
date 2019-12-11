#' Function to estimate the correlation matrix
#' @param X A sample by variable data matrix
#' @return The estimated correlation matrix
#' @description \code{cggm.corr} estimates the correlation between (potentially censored) continuous variables using the framework
#' of censored Gaussian graphical models.
#' @details When no variable is censored, \code{cggm.corr} estimates the correlation using Pearson's correlation coefficient.
#' When there are censored variables, a univariate Tobit model is first used to estimate the marginal distribution of
#' each variable. Pairwise correlations are estimated using maximum likelihood after plugging in the marginal mean
#' and standard deviation estimates.
#'
#' The correlation matrix based on pairwise estimates is not necessarily positive definite. When this is the case,
#' the function \code{\link[Matrix]{nearPD}} is used to compute the nearest positive definite correlation matrix with a
#' tolerance 1e-04 for enforcing positive definiteness.
#' @import censReg
#' @references
#' Ma, Jing. Joint Microbial and Metabolomic Network Estimation with the Censored Gaussian Graphical Model. Technical Report. 2019.
#'@examples
#' library(MASS)
#' p <- 20
#' S <- diag(1,p)
#' for (j in 1:(p-1)){
#'    S[j,j+1] <- S[j+1,j] <- 0.5
#' }
#' X <- mvrnorm(n=100, mu=rep(0,p), Sigma=S)
#' X_cens <- apply(X,2,function(a) ifelse(a>1,a,0))
#' S_hat <- cggm.corr(X_cens)
#'
#' @export
cggm.corr <- function(X){
  obj <- cgm.marginal(X)
  if (sum(obj$y_L==-Inf)==ncol(X)){
    message('All variables are continuous; use pearson correlation')
    hat_sigma <- cor(X)
  } else if (sum(obj$y_L==-Inf)==0){
    message('All variables are truncated.')
    hat_sigma <- cgm.covariance.mixed(X1=NULL, X2=X)$corr
  } else {
    new.order <- order(c(which(obj$y_L==-Inf),which(obj$y_L>-Inf)))
    X1 <- matrix(X[,which(obj$y_L==-Inf)],nrow=nrow(X))
    X2 <- matrix(X[,which(obj$y_L>-Inf)], nrow=nrow(X))
    b <- cgm.covariance.mixed(X1, X2)
    hat_sigma <- b$corr[new.order,new.order]
  }
  hat_sigma
}
