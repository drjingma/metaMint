#' compute conditional l-values in the noisy bipartite stochastic block model
#'
#' @param dataMatrix observation matrix
#' @param Z1 a node clustering for rows
#' @param Z2 a node clustering for columns
#' @param theta list of parameters for a noisy bipartite stochastic block model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} 
#'
#' @return conditional l-values in the noisy stochastic block model
lvaluesNobiSBM <- function(dataMatrix, Z1, Z2, theta, modelFamily='Gauss'){
  n1 <- length(Z1)
  n2 <- length(Z2)
  Q1 <- length(theta$alpha1)
  Q2 <- length(theta$alpha2)

  Z1_matrix <- matrix(0, nrow=Q1, ncol=n1)
  ind <- matrix(c(Z1,1:n1),ncol=2)
  Z1_matrix[ind] <- 1
  Z2_matrix <- matrix(0, nrow=Q2, ncol=n2)
  ind <- matrix(c(Z2,1:n2),ncol=2)
  Z2_matrix[ind] <- 1

  lval <- matrix(0, n1, n2)
  fnu0data <- modelDensity(dataMatrix, theta$nu0, modelFamily)
  for (q in 1:Q1){
    for (l in 1:Q2){
      mask <-  outer(Z1_matrix[q, ], Z2_matrix[l, ], '*')
      part_f0 <- fnu0data*(1-theta$pi[q,l])
      lval <- lval + part_f0/ (modelDensity(dataMatrix, 
                                            c(theta$nu$mean[q,l],theta$nu$sd[q,l]), 
                                            modelFamily) * theta$pi[q,l] + part_f0)*mask
    }
  }

  return(lval)
}




#' auxiliary function for the computation of q-values
#'
#' @param theta list of parameters for a noisy bipartite stochastic block model
#' @param ind indicators for a pair of latent blocks
#' @param t l-values
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} 
q_delta_ql_biSBM <- function(theta, ind, t, modelFamily='Gauss'){
  # for a given (q,l)
  q <- ind[1]
  l <- ind[2]
  if (modelFamily=='Gauss'){
    mu0 <- theta$nu0[1]
    sigma0 <- theta$nu0[2]
    mu <- theta$nu$mean[q,l]
    sigma <- theta$nu$sd[q,l]

    nb_t <- length(t)
    res <- matrix(0, nb_t, 2)
    res[t==1,] <- 1
    ind_t <- (t>0)&(t<1)
    if (length(ind_t)>0){
      a <- sigma^(-2)-sigma0^(-2) # scalar
      b <- -2*(mu/sigma^2 -mu0/sigma0^2)    # scalar
      cVec <- mu^2/sigma^2 -mu0^2/sigma0^2 + 2*log(sigma/sigma0*(1/theta$pi[q,l]-1)*(1/t-1))    # vector
      if(a!=0){
        res[ind_t,] <- res[ind_t,] + (a<0)
        ind_2 <- (b^2>(4*a*cVec))&ind_t
        if(sum(ind_2)>0){
          z <- (-b+matrix(sqrt(pmax(0,b^2-4*a*cVec[ind_2])),ncol=1)%*%c(1,-1))/(2*a) # matrix sum(ind_2) x 2
          res[ind_2,2] <- res[ind_2,2] + stats::pnorm(z[,1], mu, sigma) - stats::pnorm(z[,2], mu, sigma)
          res[ind_2,1] <- res[ind_2,1] + stats::pnorm(z[,1], mu0, sigma0) - stats::pnorm(z[,2], mu0, sigma0)
        }
      }else{
        if (b!=0){
          res[ind_t,2] <- if (b<0) 1-stats::pnorm(-cVec[ind_t]/b, mu, sigma) else stats::pnorm(-cVec[ind_t]/b, mu, sigma)
          res[ind_t,1] <- if (b<0) 1-stats::pnorm(-cVec[ind_t]/b, mu0, sigma0) else stats::pnorm(-cVec[ind_t]/b, mu0, sigma0)
        }else{
          res[ind_t,] <- 1*(t[ind_t] >= 1- theta$pi[q,l])
        }
      }
    }
  }
  # if (modelFamily=='Gamma')  {
  #   nu0 <- theta$nu0
  #   nu <- theta$nu[ind, ]
  #   ## warning si nu0[1] ! = nu[1]
  #   if (nu0[1] != nu[1]){
  #     stop('shapes between null and alternative distribution are different')
  #     } else {
  #     nb_t <- length(t)
  #     res <- matrix(0, nb_t, 2)
  #     res[t==1,] <- 1
  #     ind_t <- (t>0)&(t<1)
  #     if (length(ind_t)>0){
  #       if (nu0[2]>nu[2]){
  #         res[ind_t,1] <- 1- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$pi[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu0[1] , nu0[2])
  #         res[ind_t,2] <- 1- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$pi[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu[1], nu[2])
  #       }
  # 
  #       if (nu0[2]< nu[2]){
  #         res[ind_t,1] <- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$pi[ind]-1)*(nu0[2]/nu[2])^nu0[1]),  nu0[1] , nu0[2])
  #         res[ind_t,2] <- stats::pgamma( (1/(nu0[2]-nu[2]))*log((1/t[ind_t]-1)*(1/theta$pi[ind]-1)*(nu0[2]/nu[2])^nu0[1]), nu[1] , nu[2])
  #       }
  # 
  #       if (nu0[2]==nu[2]){
  #         res[ind_t,1] <- 1*(1-theta$pi[ind] <= t[ind_t])
  #         res[ind_t,2] <- 1*(1-theta$pi[ind] <= t[ind_t])
  #       }
  #     }
  #   }
  # }
  return(res)
}




#' compute q-values in the noisy bipartite stochastic block model
#'
#' @param Z1 a node clustering for rows
#' @param Z2 a node clustering for columns
#' @param theta list of parameters for a noisy bipartite stochastic block model
#' @param lvalues conditional l-values in the noisy bipartite stochastic block model
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}
#'
#' @return q-values in the noisy bipartite stochastic block model
qvaluesNobiSBM <- function(Z1, Z2, theta, lvalues, modelFamily='Gauss'){
  Q1 <- length(theta$alpha1)
  Q2 <- length(theta$alpha2)
  num <- den <- result <- matrix(0, length(Z1), length(Z2))
  for (q in 1:Q1){
    for (l in 1:Q2){
      q01_lvalues_seq <- theta$alpha1[q]*theta$alpha2[l]*q_delta_ql_biSBM(theta, c(q,l), as.vector(lvalues), modelFamily)
      f <- 1 
      num <- num + f*(1-theta$pi[q,l]) * q01_lvalues_seq[,1]
      den <- den + f*theta$pi[q,l]*q01_lvalues_seq[,2]
    }
  }
  den <- den + num
  ind <- (den!=0)
  result[ind] <- num[ind]/den[ind]
  return(result)
}

#' new graph inference procedure
#'
#' @details graph inference procedure based on conditional q-values in the noisy bipartite stochastic block model.
#' It works in the Gaussian model
#'
#' @param dataMatrix observed adjacency matrix, n1 x n2 matrix
#' @param nodeClusteringRow n1-vector of hard node Clustering
#' @param nodeClusteringCol n2-vector of hard node Clustering
#' @param theta parameter of the noisy bipartite stochastic block model
#' @param alpha confidence level
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss} 
#'
#' @return a list with:
#'  \describe{
#'     \item{\code{A}}{resulting binary adjacency matrix}
#'     \item{\code{qvalues}}{matrix with conditional q-values in the noisy bipartite stochastic block model}
#'  }
#' @export
#'
#' @examples
#' set.seed(1)
#' theta <- list(pi=c(.5,.5), w=c(.8,.1,.2), nu0=c(0,1), nu=matrix(c(-1,5,10, 1,1,1), ncol=2))
#' obs <- rnsbm(n=30, theta)
#' # res_gauss <- fitNobiSBM(obs$dataMatrix, nbCores=1)
#' resGraph <- graphInference(obs$dataMatrix, res_gauss[[2]]$clustering, theta, alpha=0.05)
#' sum((resGraph$A))/2 # nb of derived edges
#' sum(obs$latentAdj)/2 # correct nb of edges
graphInferenceNobiSBM <- function(dataMatrix, nodeClusteringRow, nodeClusteringCol, theta, alpha=0.05, modelFamily='Gauss'){

  lval_results <- lvaluesNobiSBM(dataMatrix, nodeClusteringRow, nodeClusteringCol, theta, modelFamily)
  qval_results <- qvaluesNobiSBM(nodeClusteringRow, nodeClusteringCol, theta, lval_results, modelFamily)

  A <- 1*(qval_results<alpha)
  return(list(A=A, qvalues=qval_results))
}
#' graph inference via the Benjamini-Hochberg multiple testing procedure
graphInferenceBH <- function(dataMatrix,alpha=0.05){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  pvals <- 2*pnorm(abs(as.vector(dataMatrix)),lower.tail = FALSE)
  qval_results <- p.adjust(pvals,"BH")
  A <- matrix(1*(qval_results<alpha), n1, n2)
  return(list(A=A,qvalues=qval_results))
}

graphInferenceStorey <- function(dataMatrix,alpha=0.05){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  pvals <- 2*pnorm(abs(as.vector(dataMatrix)),lower.tail = FALSE)
  qval_results <- qvalue::qvalue(pvals,fdr.level=alpha)$qvalues
  A <- matrix(1*(qval_results<alpha), n1, n2)
  return(list(A=A,qvalues=qval_results))
}


#' graph inference via the multiple testing procedure in Sun and Cai (07)
graphInferenceSC <- function(dataMatrix,alpha=0.05){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  A <- matrix(0, n1, n2)
  adaptiveZ <- adaptZ.func(as.vector(dataMatrix), alpha, 0.1, model = 'Gauss01')
  A[adaptiveZ$re] <- 1
  return(A)
}

#' graph inference via the multiple testing procedure in Sun and Cai (09)
#' which accounts for dependency via a latent HMM model
#' LIS is short for local index of significance
#' @param ncomp number of component in the mixture model. Default is 2.
graphInferenceLIS <- function(dataMatrix,ncomp=2){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  A <- matrix(0, n1, n2)
  res <- em.hmm(as.vector(dataMatrix), L=ncomp, maxiter=500)
  return(res)
}



#' plot the data matrix, the inferred graph and/or the true binary graph
#'
#' @param dataMatrix observed data matrix
#' @param inferredGraph graph inferred by the multiple testing procedure via graphInference()
#' @param binaryTruth true binary graph
#'
#' @return a list of FDR and TDR values, if possible
#' @export
plotGraphsNobiSBM <- function(dataMatrix=NULL, inferredGraph=NULL, binaryTruth=NULL){
  res <- NULL
  if (!is.null(dataMatrix))
    plt_data <- pheatmap::pheatmap(dataMatrix,cluster_rows = F,cluster_cols = F, color = RColorBrewer::brewer.pal(9,'YlGn'),
            main='data matrix',silent = TRUE)[[4]]
  xlab2 <- 'inferred graph'
  xlab3 <- 'true binary graph'
  if ((!is.null(binaryTruth))&(!is.null(inferredGraph))){
    truthVec <- as.vector(binaryTruth)
    inferredGraphVec <- as.vector(inferredGraph)
    FDR <- sum(inferredGraphVec[truthVec==0])/ sum(inferredGraphVec)
    TDR <- sum(inferredGraphVec[truthVec==1])/ sum(truthVec==1)
    xlab2 <- paste('inferred graph, FDR=', round(FDR, digits=3), sep='')
    xlab3 <- paste('true binary graph, TDR=', round(TDR, digits=3), sep='')
    res <- list(FDR=FDR, TDR=TDR)
  }
  if (!is.null(inferredGraph))
    plt_inf <- pheatmap::pheatmap(inferredGraph, cluster_rows = F,cluster_cols = F, col = RColorBrewer::brewer.pal(9,'YlGn'),
            main=xlab2,silent = TRUE)[[4]]
  if (!is.null(binaryTruth))
    plt_truth <- pheatmap::pheatmap(binaryTruth, cluster_rows = F,cluster_cols = F, col = RColorBrewer::brewer.pal(9,'YlGn'),
            main=xlab3,silent = TRUE)[[4]]

  if (!is.null(inferredGraph)){
    print(plot_grid(plt_data,plt_inf,plt_truth))
  } else{
    print(plot_grid(plt_data,plt_truth))
  }

  return(res)
}

# StructDiff <- function(inferredGraph, binaryTruth){
#   truthVec <- as.vector(binaryTruth)
#   inferredGraphVec <- as.vector(inferredGraph)
#   FDR <- sum(inferredGraphVec[truthVec==0])/ sum(inferredGraphVec)
#   TDR <- sum(inferredGraphVec[truthVec==1])/ sum(truthVec==1)
#   return(list(FDR=FDR,TDR=TDR))
# }

