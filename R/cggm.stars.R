#' Model selection for the censored Gaussian graphical model
#' @param est An object from fitting \code{cggm.pcorr}.
#' @param stars.thresh The variability threshold in stars. The default value is \code{0.1}. An alternative value is \code{0.05}.
#' @param stars.subsample.ratio The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size.
#' @param rep.num The number of subsamplings. The default value is \code{20}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @description This function implements the regularization parameter selection for high-dimensional censored Gaussian graphical models based on stability.
#' @details
#' Stability approach to regularization selection (\href{https://papers.nips.cc/paper/3966-stability-approach-to-regularization-selection-stars-for-high-dimensional-graphical-models.pdf}{stars}) is a way
#' to select optimal regularization parameter for high-dimensional graphical models. It selects the optimal graph that has the smallest variability of subsamplings.
#' It also provides an additional estimated graph which merges the subsampled graphs using the frequency counts.
#'
#' The subsampling procedure in stars may not be very efficient.
#'
#' @return A list with items
#' \item{refit}{The optimal graph selected from the graph path}
#' \item{opt.icov}{The optimal precision matrix from the path}
#' \item{merge}{The graph path estimated by merging the subsampling paths}
#' \item{variability}{The variability along the subsampling paths}
#' \item{opt.index}{The index of the selected regularization parameter}
#' \item{opt.lambda}{The selected regularization parameter}
#' \item{opt.sparsity}{The sparsity level of \code{refit}}
#' @references
#' Liu, Han, Kathryn Roeder, and Larry Wasserman. "Stability approach to regularization selection (stars) for high dimensional graphical models." Advances in neural information processing systems. 2010.
#'
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
#' sel_models <- cggm.stars(Sinv_hat, rep.num=5)
#' @export
cggm.stars <- function(est, stars.thresh=0.1,stars.subsample.ratio=NULL, rep.num=20, verbose=TRUE){
  n = nrow(est$data)
  d = ncol(est$data)
  nlambda = length(est$lambda)

  if(is.null(stars.subsample.ratio))
  {
    if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
    if(n<=144) stars.subsample.ratio = 0.8
  }

  est$merge = list()
  for(i in 1:nlambda) est$merge[[i]] = matrix(0,d,d)

  for(i in 1:rep.num)
  {
    if(verbose)
    {
      mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/rep.num), "%"), collapse="")
      cat(mes, "\r")
      flush.console()
    }
    ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)

    tmp = cggm.pcorr(est$data[ind.sample,], lambda = est$lambda, method='glasso')$path

    for(i in 1:nlambda)
      est$merge[[i]] = est$merge[[i]] + tmp[[i]]

    rm(ind.sample,tmp)
    gc()
  }

  if(verbose){
    mes = "Conducting Subsampling....done.                 "
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }

  est$variability = rep(0,nlambda)
  for(i in 1:nlambda){
    est$merge[[i]] = est$merge[[i]]/rep.num
    est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
  }

  est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
  est$refit = est$path[[est$opt.index]]
  est$opt.lambda = est$lambda[est$opt.index]
  est$opt.sparsity = est$sparsity[est$opt.index]
  est$opt.icov = est$icov[[est$opt.index]]

  class(est) = "select"
  return(est)
}
