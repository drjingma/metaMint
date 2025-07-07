#' transform a pair of nodes (i,j) into an identifying integer
#'
#' Associates an identifying integer with a pair of nodes (i,j)
#'
#' @details returns the row number of the matrix build by listNodePairs(n)
#'     containing the pair (i,j)
#'
#' @param i scalar or vector
#' @param j scalar or vector, same length as i
#' @param n number of vertices
#' @param directed booelan to indicate whether the model is directed or undirected
# convertNodePair <- function(i,j,n, directed){
#   if (sum((i>n) | (j>n))>0){
#     stop("Your index is out of range")
#   }
#   if (directed){#directed case
#     dyads <- (i-1)*(n-1)+j-(i<j)
#   } else {#undirected case
#     dyads <- c(0,cumsum((n-1):1))[pmin(i,j)] + abs(j-i)
#   }
#   return(dyads)
# }



#' returns a list of all possible node pairs (i,j)
#'
#' @param n1 number of nodes on side one
#' @param n2 number of nodes on side two
#' @param directed indicates if the graph is directed
#'
#' @return a 2-column matrix with all possible node pairs (i,j)
listNodePairs <- function(n1, n2, directed=FALSE){
  return(expand.grid(1:n1,1:n2))
}


#' transform a pair of block identifiers (q,l) into an identifying integer
#'
#' this is the inverse function of convertGroupPairIdentifierBipartite()
#'
#' @param q indicator of a latent block
#' @param l indicator of a latent block; can be a vector
#' @param Q1 number of latent blocks on side one
#' @param Q2 number of latent blocks on side two
#' @param directed indicates if the graph is directed
convertGroupPair <- function(q, l, Q1, Q2){
  all <- listNodePairs(Q1,Q2)
  index <- match(sapply(l, function(a) paste(q, a, collapse =" ")),apply(all, 1, paste, collapse =" "))
  return(index)
}


#' takes a scalar index of a group pair (q,l) and returns the values q and l
#'
#' this is the inverse function of convertGroupPair()
#'
#' @param ind_ql indicator for a pair of latent blocks
#' @param Q1 number of latent blocks on side one
#' @param Q2 number of latent blocks on side two
# convertGroupPairIdentifierBipartite <- function(ind_ql, Q1, Q2){
#   q <- ind_ql %% Q1
#   if (q==0){
#     q <- Q1
#   }
#   l <- ceiling(ind_ql/Q1)
#   return(c(q,l))
# }

#' takes a scalar indice of a group pair (q,l) and returns the values q and l
#'
#' this is the inverse function of convertGroupPair()
#'
#' @param ind_ql indicator for a pair of latent blocks
#' @param Q number of latent blocks
convertGroupPairIdentifier <- function(ind_ql, Q){
  w <- cumsum((Q-1):1)
  q <- which.max(ind_ql<=w)
  w <- c(0, w)
  l <- ind_ql - w[q] + q
  return(c(q,l))
}


#' corrects values of the variational parameters tau that are too close to the 0 or 1
#'
#' @param tau variational parameters
correctBeta <- function(tau){
  tau <- pmin(tau,.Machine$double.xmax)
  tau <- pmax(tau,.Machine$double.xmin)
  tau <- tau/sum(tau)
  tau <- pmin(tau,1-1e-7)
  tau <- pmax(tau,1e-7)
  tau <- tau/sum(tau)
  
  return(tau)
}


#' evaluate the density in the current model
#'
#' @param x vector with points where to evaluate the density
#' @param nu distribution parameter
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}, \code{Poisson}
modelDensity <- function(x, nu, modelFamily='Gauss'){
  if (modelFamily=='Gauss')
    res <- stats::dnorm(x, nu[1], nu[2])
  if (modelFamily=='Gamma')
    res <- stats::dgamma(x, nu[1], nu[2])
  res[res<=.Machine$double.eps] <- .Machine$double.eps
  return(res)
}


#' Evaluate beta_q_1*beta_l_2 in the noisy bipartite stochastic block model
#'
#' @param q indicator of a latent block on side one
#' @param l indicator of a latent block on side two
#' @param beta variational parameters in a list of length 2
getBetaql <- function(q, l, beta){
  beta1 <- beta[[1]]
  beta2 <- beta[[2]]
  n1 <- ncol(beta1)
  n2 <- ncol(beta2)
  Q1 <- nrow(beta1)
  Q2 <- nrow(beta2)
  
  # would like to get a matrix output
  if (Q1 == 1 && Q2 == 1) # one block on both sides
    betaql <- matrix(1,n1,n2)
  else{
    betaql <- outer(beta1[q, ],beta2[l, ], '*')
  }
  return(betaql)
}

#' Function to calculate test statistics
#' @param X,Y feature by sample matrices
teststat_general <- function(X,Y=NULL){
  n1 <- ncol(X)
  if (is.null(Y)){
    Y <- X
  }
  empcor <- cor(t(X),t(Y),method = 'pearson')
  X.scale <- t(scale(t(X), center = T, scale = T))
  Y.scale <- t(scale(t(Y), center = T, scale = T))
  Theta <- matrix(0, nrow(X), nrow(Y))
  for (i in 1:nrow(Theta)){
    for (j in 1:ncol(Theta)){
      Theta[i,j] <- mean((2 * X.scale[i,] * Y.scale[j,] - empcor[i,j] * X.scale[i,]^2 - empcor[i,j] * Y.scale[j,]^2)^2)
    }
  }
  dataMatrix <- 2*empcor*sqrt(n1/Theta)
  
  return(list(mean = empcor, sd = Theta, statistic = dataMatrix))
}
#' Function to calculate two-sample test statistics 
teststat_general_2sample <- function(sample1, sample2){
  n1 <- ncol(sample1$X)
  n2 <- ncol(sample2$X)
  
  out1 <- teststat_general(sample1$X,sample1$Y)
  out2 <- teststat_general(sample2$X,sample2$Y)
  
  dataMatrix <- 2*(out1$mean - out2$mean)/sqrt(out1$sd/n1 + out2$sd/n2)
  
  dataMatrix
}









