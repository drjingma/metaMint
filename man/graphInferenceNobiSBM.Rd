% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testProcedure_biSBM.R
\name{graphInferenceNobiSBM}
\alias{graphInferenceNobiSBM}
\title{new graph inference procedure}
\usage{
graphInferenceNobiSBM(
  dataMatrix,
  nodeClusteringRow,
  nodeClusteringCol,
  theta,
  alpha = 0.05,
  modelFamily = "Gauss"
)
}
\arguments{
\item{dataMatrix}{observed adjacency matrix, n1 x n2 matrix}

\item{nodeClusteringRow}{n1-vector of hard node Clustering}

\item{nodeClusteringCol}{n2-vector of hard node Clustering}

\item{theta}{parameter of the noisy bipartite stochastic block model}

\item{alpha}{confidence level}

\item{modelFamily}{probability distribution for the edges. Possible values:
\code{Gauss}}
}
\value{
a list with:
 \describe{
    \item{\code{A}}{resulting binary adjacency matrix}
    \item{\code{qvalues}}{matrix with conditional q-values in the noisy bipartite stochastic block model}
 }
}
\description{
new graph inference procedure
}
\details{
graph inference procedure based on conditional q-values in the noisy bipartite stochastic block model.
It works in the Gaussian model
}
\examples{
set.seed(1)
theta <- list(pi=c(.5,.5), w=c(.8,.1,.2), nu0=c(0,1), nu=matrix(c(-1,5,10, 1,1,1), ncol=2))
obs <- rnsbm(n=30, theta)
# res_gauss <- fitNobiSBM(obs$dataMatrix, nbCores=1)
resGraph <- graphInference(obs$dataMatrix, res_gauss[[2]]$clustering, theta, alpha=0.05)
sum((resGraph$A))/2 # nb of derived edges
sum(obs$latentAdj)/2 # correct nb of edges
}
