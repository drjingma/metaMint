% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxiliaryFunctions_biSBM.R
\name{listNodePairs}
\alias{listNodePairs}
\title{transform a pair of nodes (i,j) into an identifying integer}
\usage{
listNodePairs(n1, n2, directed = FALSE)
}
\arguments{
\item{n1}{number of nodes on side one}

\item{n2}{number of nodes on side two}

\item{directed}{indicates if the graph is directed}

\item{i}{scalar or vector}

\item{j}{scalar or vector, same length as i}

\item{n}{number of vertices}
}
\value{
a 2-column matrix with all possible node pairs (i,j)
}
\description{
Associates an identifying integer with a pair of nodes (i,j)
}
\details{
returns the row number of the matrix build by listNodePairs(n)
    containing the pair (i,j)
}
