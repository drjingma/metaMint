% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMalgorithm_biSBM.R
\name{JEvalMstep}
\alias{JEvalMstep}
\title{Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution}
\usage{
JEvalMstep(VE, mstep, dataMatrix, modelFamily)
}
\arguments{
\item{VE}{list with variational parameters beta and rho}

\item{mstep}{list with current model parameters and additional auxiliary terms}

\item{dataMatrix}{observed dense adjacency matrix}

\item{modelFamily}{probability distribution for the edges. Possible values:
\code{Gauss}, \code{Gamma}}

\item{param}{parameters of the Gamma distribution}

\item{L}{weighted mean of log(data)}

\item{M}{weighted mean of the data}

\item{param.old}{parameters of the Gamma distribution}

\item{epsilon}{threshold for the stopping criterion}

\item{nb.iter.max}{maximum number of iterations}
}
\value{
updated parameters of the Gamma distribution
evaluate the objective in the Gamma model

value of the lower bound of the log-likelihood function
compute the MLE in the Gamma model using the Newton-Raphson method

updated parameters of the Gamma distribution
evaluation of the objective in the Gauss model
this evaluates both the ELBO and the complete data log-likelihood

value of the ELBO and the complete log likelihood function
}
\description{
Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution
}
