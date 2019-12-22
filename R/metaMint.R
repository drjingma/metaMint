#' Joint estimation of metabolite and microbial interaction networks.
#'
#' The metaMint package provides functions for estimating the correlation and
#' partial correlations between microbial species and metabolites.
#' @details
#' \tabular{ll}{
#'  Package: \tab metaMint\cr
#'  Type: \tab Package\cr
#'  Version: \tab 0.1.0\cr
#'  Date: \tab 2019-12-21\cr
#'  License: \tab GPL (>=2)\cr
#' }
#' @author Jing Ma <jingma@fredhutch.org>
#' @examples
#' p <- 20
#' S <- diag(1,p)
#' for (j in 1:(p-1)){
#'    S[j,j+1] <- S[j+1,j] <- 0.5
#' }
#' X <- mvrnorm(n=100, mu=rep(0,p), Sigma=S)
#' X_cens <- apply(X,2,function(a) ifelse(a>1,a,0))
#' S_hat <- cggm.corr(X_cens)
#' Sinv_hat <- cggm.pcorr(X_cens,c(0.2,0.1,0.05),'glasso')
#'
#' @docType package
#' @name metaMint
#' @import MASS
#' @importFrom stats cor dnorm optimize pnorm sd
#' @importFrom utils flush.console
#' @importFrom mvtnorm dmvnorm pmvnorm
#' @importFrom clime clime
#' @importFrom compositions clr
#' @importFrom Matrix nearPD
NULL
