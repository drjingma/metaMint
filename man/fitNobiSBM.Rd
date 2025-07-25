% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMalgorithm_biSBM.R
\name{fitNobiSBM}
\alias{fitNobiSBM}
\title{VEM algorithm to adjust the noisy bipartite stochastic block model to an observed dense adjacency matrix
a number of passes can be used to deal with potential sensitivity to initialization.}
\usage{
fitNobiSBM(
  dataMatrix,
  model = "Gauss0",
  sbmSize = list(Q1 = 1:5, Q2 = 1:5, explor = 1.5),
  exclusive_Rows_Cols = TRUE,
  filename = NULL,
  initParam = list(nbOfbeta = NULL, nbOfPointsPerbeta = NULL, maxNbOfPasses = NULL,
    minNbOfPasses = 1),
  nbCores = parallel::detectCores()
)
}
\arguments{
\item{dataMatrix}{observed dense adjacency matrix}

\item{model}{Implemented models:
\describe{
  \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
  \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
  \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
  \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
  \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
  \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
  \item{\code{Gauss02distr}}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
  \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
  \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
  \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
  \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
}}

\item{sbmSize}{list of parameters determining the size of SBM (the number of latent blocks) to be explored. Currently restrict the number of blocks on both sides to be equal.
\describe{
  \item{\code{Q1}}{range for the number of latent blocks in Rows}
  \item{\code{Q2}}{range for the number of latent blocks in Columns}
  \item{\code{explor}}{if \code{Qmax} is not provided, then \code{Qmax} is automatically determined as \code{explor} times the number of blocks where the ICL is maximal}
}}

\item{exclusive_Rows_Cols}{whether to set the number of latent blocks in Rows to be the same as the number of latent blocks in Columns.}

\item{filename}{results are saved in a file with this name (if provided)}

\item{initParam}{list of parameters that fix the number of initialization
\describe{
  \item{\code{nbOfbeta}}{number of initial points for the node clustering (i. e. for the variational parameters \code{beta})}
  \item{\code{nbOfPointsPerbeta}}{number of initial points of the latent binary graph}
  \item{\code{maxNbOfPasses}}{maximum number of passes through the SBM models, that is, passes from \code{Qmin} to \code{Qmax} or inversely}
  \item{\code{minNbOfPasses}}{minimum number of passes through the SBM models}
}}

\item{nbCores}{number of cores used for parallelization}
}
\value{
Returns a list of estimation results for all numbers of latent blocks considered by the algorithm.
Every element is a list composed of:
\describe{
  \item{\code{theta}}{estimated parameters of the noisy bipartite stochastic block model; a list with the following elements:
  \describe{
    \item{\code{alpha1}}{parameter estimate of alpha1}
    \item{\code{alpha2}}{parameter estimate of alpha2}
    \item{\code{pi}}{parameter estimate of pi}
    \item{\code{nu0}}{parameter estimate of nu0}
    \item{\code{nu}}{parameter estimate of nu}
   }}
  \item{\code{clustering_row}}{node clustering for rows obtained by the noisy bipartite stochastic block model, more precisely, a hard clustering given by the
  maximum a posterior estimate of the variational parameters \code{sbmParam$edgeProba[[1]]}}
  \item{\code{clustering_col}}{node clustering for columns obtained by the noisy bipartite stochastic block model, more precisely, a hard clustering given by the
  maximum a posterior estimate of the variational parameters \code{sbmParam$edgeProba[[2]]}}
  \item{\code{sbmParam}}{further results concerning the latent binary stochastic block model. A list with the following elements:
  \describe{
    \item{\code{Q1,Q2}}{number of latent blocks in the noisy bipartite stochastic block model}
    \item{\code{clusterProba}}{a list of length two, consisting of soft clustering given by the conditional probabilities of a node to belong to a given latent block.
    In other words, these are the variational parameters \code{beta1}; (Q1 x n1)-matrix and \code{beta2}; (Q2 x n2)-matrix}
    \item{\code{edgeProba}}{conditional probabilities \code{rho} of an edges given the node memberships of the interacting nodes; (N_Q x N)-matrix}
    \item{\code{ICL}}{value of the ICL criterion at the end of the algorithm}
    }}
  \item{\code{convergence}}{a list of convergence indicators:
  \describe{
    \item{\code{J}}{value of the lower bound of the log-likelihood function at the end of the algorithm}
    \item{\code{complLogLik}}{value of the complete log-likelihood function at the end of the algorithm}
    \item{\code{converged}}{indicates if algorithm has converged}
    \item{\code{nbIter}}{number of iterations performed}
 }}
}
}
\description{
\code{fitNobiSBM()} estimates model parameters of the noisy bipartite stochastic block model and provides a clustering of the nodes
}
\details{
\code{fitNobiSBM()} supports different probability distributions for the edges and can estimate the number of node blocks
}
\examples{
n1 <- 10
n2 <- 15
Q1 <- 2
Q2 <- 2
theta <- list(alpha1 = c(0.5, 0.5), 
              alpha2 = c(0.4, 0.6),
              nu0 = c(0,.1),
              nu = list(mean = matrix(2, Q1, Q2), 
                        sd = matrix(0.1, Q1, Q2)),  
              pi=matrix(0.01,Q1,Q2) )
diag(theta$pi) <- 0.3
obs <- rnbisbm(n1, n2, theta, modelFamily='Gauss')
res <- fitNobiSBM(obs$dataMatrix, sbmSize = list(Q1 =1:3, Q2=1:3),
       initParam=list(nbOfbeta=1, nbOfPointsPerbeta=1), nbCores=1)
}
