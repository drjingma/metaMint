% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VEMinitialization_biSBM.R
\name{betaDown}
\alias{betaDown}
\title{Create new initial values by merging pairs of groups of provided beta}
\usage{
betaDown(beta, nbOfMerges, by = "row")
}
\arguments{
\item{beta}{soft node clustering}

\item{nbOfMerges}{number of required merges of blocks}
}
\value{
a list of length nbOfMerges (at most) of initial points for beta
}
\description{
Create nbOfMerges new initial values by merging nbOfMerges (or all possible) pairs of groups of provided beta
}
