% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cggm.stars.R
\name{cggm.stars}
\alias{cggm.stars}
\title{Model selection for the censored Gaussian graphical model}
\usage{
cggm.stars(
  est,
  stars.thresh = 0.1,
  stars.subsample.ratio = NULL,
  rep.num = 20,
  verbose = TRUE
)
}
\arguments{
\item{est}{An object from fitting \code{cggm.pcorr}.}

\item{stars.thresh}{The variability threshold in stars. The default value is \code{0.1}. An alternative value is \code{0.05}.}

\item{stars.subsample.ratio}{The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size.}

\item{rep.num}{The number of subsamplings. The default value is \code{20}.}

\item{verbose}{If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.}
}
\value{
A list with items
\item{refit}{The optimal graph selected from the graph path}
\item{opt.icov}{The optimal precision matrix from the path}
\item{merge}{The graph path estimated by merging the subsampling paths}
\item{variability}{The variability along the subsampling paths}
\item{opt.index}{The index of the selected regularization parameter}
\item{opt.lambda}{The selected regularization parameter}
\item{opt.sparsity}{The sparsity level of \code{refit}}
}
\description{
This function implements the regularization parameter selection for high-dimensional censored Gaussian graphical models based on stability.
}
\details{
Stability approach to regularization selection (\href{https://papers.nips.cc/paper/3966-stability-approach-to-regularization-selection-stars-for-high-dimensional-graphical-models.pdf}{stars}) is a way
to select optimal regularization parameter for high-dimensional graphical models. It selects the optimal graph that has the smallest variability of subsamplings.
It also provides an additional estimated graph which merges the subsampled graphs using the frequency counts.

The subsampling procedure in stars may not be very efficient.
}
\examples{
library(MASS)
library(glasso)
p <- 20
S <- diag(1,p)
for (j in 1:(p-1)){
   S[j,j+1] <- S[j+1,j] <- 0.5
}
X <- mvrnorm(n=100, mu=rep(0,p), Sigma=S)
X_cens <- apply(X,2,function(a) ifelse(a>1,a,0))
Sinv_hat <- cggm.pcorr(X_cens,c(0.2,0.1,0.05),'glasso')
sel_models <- cggm.stars(Sinv_hat, rep.num=5)
}
\references{
Liu, Han, Kathryn Roeder, and Larry Wasserman. "Stability approach to regularization selection (stars) for high dimensional graphical models." Advances in neural information processing systems. 2010.
}
