#' Function to estimate the precision matrix based on the censored Gaussian graphical model
#' @param X A sample by variable (\code{n} by \code{p}) data matrix
#' @param lambda A sequence of decreasing positive numbers to control the regularization.
#' @param method Method for estimating the precision matrix. Could be either \code{clime} or \code{glasso}.
#' @description
#' \code{cggm.precision} estimates the precision, or inverse covariance, matrix based on the censording Gaussian graphical model.
#' A correlation matrix is first estimated using \code{cggm.corr}
#' @return A list with items
#' \item{path}{A list of \code{p} by \code{p} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}}
#' \item{icov}{A list of \code{p} by \code{p} inverse covariance matrices corresponding to \code{lambda}}
#' \item{cov}{A list of \code{p} by \code{p} estimated covariance matrices corresponding to \code{lambda}}
#' \item{lambda}{The sequencing of decreasing regulatization parameters used in estimating the precision matrices}
#' \item{data}{The \code{n} by \code{p} data matrix from the input}
#' @examples
cggm.precision <- function(X,lambda,method=c('clime','glasso')){
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
