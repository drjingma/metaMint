#' Function to estimate the partial correlation matrix
#' @param X A sample by variable (\code{n} by \code{p}) data matrix
#' @param lambda A sequence of decreasing positive numbers to control the regularization.
#' @param method Method for estimating the precision matrix. Could be either \code{clime} or \code{glasso}.
#' @description
#' \code{cggm.pcorr} estimates the partial correlation matrix based on the censored Gaussian graphical model.
#' @return A list with items
#' \item{path}{A list of \code{p} by \code{p} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}}
#' \item{icov}{A list of \code{p} by \code{p} inverse covariance matrices corresponding to \code{lambda}}
#' \item{cov}{A list of \code{p} by \code{p} estimated covariance matrices corresponding to \code{lambda}}
#' \item{lambda}{The sequencing of decreasing regulatization parameters used in estimating the precision matrices}
#' \item{data}{The \code{n} by \code{p} data matrix from the input}
#' @details A correlation matrix is first estimated data \code{X} using \code{cggm.corr}. Given the estimated correlation matrix, the
#' \href{http://www-stat.wharton.upenn.edu/~tcai/paper/Precision-Matrix.pdf}{clime} method or the \href{https://cran.r-project.org/web/packages/glasso/index.html}{graphical lasso} method can be
#' used to estimate the precision matrix along a sequence of regularization parameters \code{lambda}.
#' @references
#' Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "Sparse inverse covariance estimation with the graphical lasso." Biostatistics 9.3 (2008): 432-441.
#'
#' Cai, Tony, Weidong Liu, and Xi Luo. "A constrained l1 minimization approach to sparse precision matrix estimation." Journal of the American Statistical Association 106.494 (2011): 594-607.
#' @import glasso
#' @examples
#' library(MASS)
#' library(glasso)
#' p <- 20
#' S <- diag(1,p)
#' for (j in 1:(p-1)){
#'    S[j,j+1] <- S[j+1,j] <- 0.5
#' }
#' X <- mvrnorm(n=100, mu=rep(0,p), Sigma=S)
#' X_cens <- apply(X,2,function(a) ifelse(a>1,a,0))
#' Sinv_hat <- cggm.pcorr(X_cens,c(0.2,0.1,0.05),'glasso')
#' @export
cggm.pcorr <- function(X,lambda,method=c('clime','glasso')){
  lambda.sorted <- sort(lambda,decreasing = TRUE)

  method <- match.arg(method)
  hat_sigma <- cggm.corr(X)

  if (method=='clime'){
    ## Fit clime
    pathList <- clime(hat_sigma, sigma=TRUE, lambda=lambda.sorted)
  }
  if (method=='glasso'){
    ## Fit the graphical lasso
    pathList <- lapply(lambda.sorted, function(m) glasso(hat_sigma, rho=m, penalize.diagonal = FALSE, approx=FALSE)$wi)
  }

  pathList <- lapply(pathList, function(a) (a+t(a))/2)
  adjList <- lapply(pathList, function(a) 1*(abs(a)>1e-08) - diag(1,ncol(X)))
  return(list(path=adjList, icov=pathList,cov=hat_sigma,lambda=lambda.sorted, data=X))
}
