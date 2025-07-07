
#' compute a list of initial points for the VEM algorithm
#'
#' compute a list of initial points of beta and rho for the VEM algorithm
#' for a given number of blocks; returns nbOfBeta*nbOfPointsPerBeta inital points
#'
#' @param Q1 number of latent blocks in the Rows
#' @param Q2 number of latent blocks in the Cols
#' @param dataMatrix observed dense adjacency matrix
#' @param nbOfBeta number of initializations for the latent block memberships
#' @param nbOfPointsPerBeta number of initializations of the latent binary graph
#' associated with each initial latent block memberships
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss02distr}}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param method type of method used when initializing the row/column block memberships
#'
#' @return list of inital points of beta and rho of length nbOfBeta*nbOfPointsPerBeta
initialPoints <- function(Q1, Q2, dataMatrix, modelFamily, model, method='kmeans'){
  beta <- initialBeta(Q1, Q2, dataMatrix, method=method)
  init <- initialRho(beta, dataMatrix, modelFamily, model)
  return(init)
}


#' compute initial values for beta
#'
#' returns a list of length nbOfBeta of initial points for beta using spectral clustering with absolute values, kmeans and random perturbations of these points
#'
#' @param Q1 number of latent blocks on side one in the noisy bipartite stochastic block model
#' @param Q2 number of latent blocks on side two in the noisy bipartite stochastic block model
#' @param dataMatrix observed dense adjacency matrix
#' @param method method used for initializing the latent block memberships
#'
#' @return a list of initial points for beta1 and beta2
initialBeta <- function(Q1, Q2, dataMatrix, method='spectral'){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  
  if (Q1==1){
    beta1 <- matrix(1,1,n1)
  }
  if (Q2==1){
    beta2 <- matrix(1,1,n2)
  }
  if (Q1>1){
    # cluster the rows
    
    if (method=='spectral'){
      trySpect <- try(spectralClustering(dataMatrix, Q1, Q2), silent=TRUE)
      if (is.list(trySpect)){
        cl <- trySpect$clustering_row
      }
    }
    
    if (method=='kmeans'){
      tryKmeans_row <- try(stats::kmeans(dataMatrix,Q1,nstart=50),silent=TRUE) 
      
      if (is.numeric(tryKmeans_row$cluster)){
        cl <- tryKmeans_row$cluster
      }
    }
    # cl_re <- 
    #   factor(cl,
    #          levels = 1:Q1,
    #          labels = unique(cl)) |> 
    #   as.character() |> 
    #   as.numeric()
    
    beta1 <- apply(t(classInd(cl, Q1)),2,correctBeta)
    
  }
  if (Q2 > 1){
    # cluster the columns
    if (method=='spectral'){
      trySpect <- try(spectralClustering(dataMatrix, Q1, Q2), silent=TRUE)
      if (is.list(trySpect)){
        cl <- trySpect$clustering_col
      }
    }
    
    if (method=='kmeans'){
      tryKmeans_col <- try(stats::kmeans(t(dataMatrix),Q2,nstart=50),silent=TRUE) 
      if (is.numeric(tryKmeans_col$cluster)){
        cl <- tryKmeans_col$cluster
      }
    }
    # cl_re <- 
    #   factor(cl,
    #          levels = 1:Q2,
    #          labels = unique(cl)) |> 
    #   as.character() |> 
    #   as.numeric()
    beta2 <- apply(t(classInd(cl, Q2)),2,correctBeta)
    
  }
  return(list(list(beta1, beta2)))
}

#' regularization of bipartite graphs
#'
#' @param A adjacency matrix
#' @param tau regularization parameter. Default is tau=3.
#'
#' @return a vector containing a node clustering into K groups
graphRegularization <- function(A, tau, weighted=TRUE){
  n1 <- nrow(A)
  n2 <- ncol(A)
  D_row <- rowSums(abs(A))
  D_row_mean <- mean(D_row)
  alpha_row <- floor(n1/D_row_mean)
  d1hat <- tau*sort(D_row,decreasing = TRUE)[alpha_row]
  
  D_col <- colSums(abs(A))
  D_col_mean <- mean(D_col)
  alpha_col <- floor(n2/D_col_mean)
  d2hat <- tau*sort(D_col,decreasing = TRUE)[alpha_col]
  
  w1 <- sapply(d1hat/D_row, function(a) min(a,1))
  w2 <- sapply(d2hat/D_col, function(a) min(a,1))
  
  A_re <- A*tcrossprod(w1,w2)
  
  return(A_re)
}


#' spectral clustering with absolute values for bipartite graphs
#'
#' performs absolute spectral clustering of an adjacency matrix
#'
#' @param A adjacency matrix
#' @param K1 number of desired clusters over rows
#' @param K2 number of desired clusters over columns
#'
#' @return a list containing a node clustering of rows into K1 groups and a node clustering of columns into K2 groups
spectralClustering <- function(A, K1, K2){
  svd.out <- svd(A,nu=K1,nv=K2)
  clustering_row <- stats::kmeans(svd.out$u[,1:K1],centers = K1,nstart = 100)$cluster
  clustering_col <- stats::kmeans(svd.out$v[,1:K2],centers = K2,nstart = 100)$cluster
  return(list(clustering_row=clustering_row, clustering_col=clustering_col))
}

#' compute initial values of rho
#'
#' for every provided initial point of beta, nbOfPointsPerBeta initial values of rho are computed
#' in the Gamma model also initial values of nu are computed
#'
#' @param beta output of initialBeta()
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss02distr}}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#' @param type type of methods to use when initializing the connectivity matrix
#'
#' @return list of initial points of beta1, beta2 and rho
initialRho <- function(listOfBeta, data, modelFamily, model){
  nbOfBeta <- length(listOfBeta)
  beta <- listOfBeta[[1]]
  beta1 <- beta[[1]]
  beta2 <- beta[[2]]
  Q1 <- nrow(beta1)
  n1 <- ncol(beta1)
  Q2 <- nrow(beta2)
  n2 <- ncol(beta2)

  init <- list(beta=lapply(1:nbOfBeta, function(k) listOfBeta[[k]]),
               rho=lapply(1:nbOfBeta, function(k) array(1, c(Q1, Q2, n1, n2)))  
               # for every beta keep one initialization with all rho equal to 1
  )
  
  if (modelFamily=='Gauss'){
    nu0 <- c(0,1)
    pvalues <- 2*(1-stats::pnorm(abs(data)))
    # pheatmap(pvalues,cluster_rows = F,cluster_cols = F)
    for (k in 1:nbOfBeta){
      k <- 1
      nu <- list() # store the parameters for the alternative distribution
      nu$mean <- matrix(0, Q1, Q2)
      nu$sd <- matrix(1, Q1, Q2)
      pi <- matrix(0, Q1, Q2)
      for (q in 1:Q1){
        for (l in 1:Q2){
          beta_ql <- getBetaql(q, l, init$beta[[k]])
          
          ss <- sum(round(beta_ql))
          if (ss>0)
            pi[q,l] <- max(c(1-2*sum(round(beta_ql)*(pvalues>0.5))/ss,0))
          beta_ql <- beta_ql*(pvalues<=0.5)
          s <- sum(beta_ql)
          if (s>0)
            nu$mean[q,l] <- sum(beta_ql*data)/s
        }
      }
      init$rho[[k]] <- getRho(Q1, Q2, pi, nu0, nu, data, modelFamily)
      
    }
  }
  
  return(init)
}



#' compute rho associated with given values of pi, nu0 and nu
#'
#' @param Q1 number of latent blocks in the noisy bipartite stochastic block model on side one
#' @param Q2 number of latent blocks in the noisy bipartite stochastic block model on side two
#' @param pi connectivity parameter in the noisy bipartite stochastic block model
#' @param nu0 null parameter in the noisy bipartite stochastic block model
#' @param nu alternative parameter in the noisy bipartite stochastic block model
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return a matrix of conditional probabilities of an edge given the node memberships of the interacting nodes
getRho <- function(Q1, Q2, pi, nu0, nu, dataMatrix, modelFamily){
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  rho <-  array(NA, c(Q1,Q2,n1,n2))
  for (q in 1:Q1){
    for (l in 1:Q2){
      rhoNumerator <- modelDensity(dataMatrix, c(nu$mean[q,l],nu$sd[q,l]), modelFamily) * pi[q,l]  # taille N
      rho[q,l, ,] <- rhoNumerator / (rhoNumerator + modelDensity(dataMatrix, nu0, modelFamily)*(1-pi[q,l]))   # taille N
    }
  }
  return(rho)
}




#' convert a clustering into a 0-1-matrix
#'
#' @param cl cluster in vector form
#' @param nbClusters number of clusters
#'
#' @return a 0-1-matrix encoding the clustering
classInd <- function (cl, nbClusters){
  nbClusters <- max(nbClusters, max(cl))
  n <- length(cl)
  x <- matrix(0, n, nbClusters)
  x[(1:n) + n*(cl-1)] <- 1
  return(x)
}


#' Construct initial values with Q groups by splitting groups of a solution obtained with Q-1 groups
#'
#' @param beta_Qm1 beta for a model with Q-1 latent blocks
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss02distr}}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#'
#' @return list of initial points of beta and rho of length nbOfBeta*nbOfPointsPerBeta
initialPointsBySplit <- function(beta_Qm1, nbOfBeta=1, dataMatrix, modelFamily, model){
  beta1 <- betaUp(beta_Qm1[[1]], nbOfBeta)
  beta2 <- betaUp(beta_Qm1[[2]], nbOfBeta)
  nbOfBeta <- min(nbOfBeta, length(beta1))
  listOfBeta <- lapply(1:nbOfBeta, function(k) list(beta1[[k]],beta2[[k]]))
  init <- initialRho(listOfBeta, dataMatrix, modelFamily, model)
  return(init)
}

#' Create new values of beta by splitting groups of provided beta
#'
#' Create nbOfSplits (or all) new values of beta by splitting nbOfSplits (or all) groups of provided beta
#'
#' @param beta a matrix consisting of soft node clustering for row or column
#' @param nbOfSplits number of required splits of blocks
#'
#' @return a list of length nbOfSplits (at most) of initial points for beta
betaUp <- function(beta, nbOfSplits=1){
  n <- ncol(beta)
  Qold <- nrow(beta) # value of Q at previous step (for which a solution is available)
  
  if (Qold==1){
    listOfBeta <- lapply(1:nbOfSplits, function(k) t(gtools::rdirichlet(n, c(0.7,0.7))) )
  }else{ # Qold>=2
    if (nbOfSplits>=Qold){ # split all the groups once
      listOfBeta <- lapply(1:Qold, function(q) addRowToBeta(beta,q))
    }else{ # always split the largest group
      largestGroup <- which.max(apply(beta,1,sum))
      if (nbOfSplits>=2){ # and split the component with largest entropy
        entropy <- beta*log(beta)
        entropy[is.na(entropy)] <- 0
        smallestEntropy <- which.min(apply(entropy,1,sum))
      } else {
        smallestEntropy <- NULL
      }
      groupsToSplit <- if (nbOfSplits>2) sample((1:Qold)[-c(largestGroup,smallestEntropy)],nbOfSplits-2,replace=F) else NULL
      
      listOfBeta <- lapply(c(largestGroup, smallestEntropy, groupsToSplit), function(q) addRowToBeta(beta,q))
    }
  }
  
  return(listOfBeta)
}




#' split group q of provided beta randomly into two
#'
#' @param beta provided beta
#' @param q index of group to split
#'
#' @return new beta
addRowToBeta <- function(beta, q){
  n <- ncol(beta)
  
  newBeta <- beta
  newBeta[q,] <- newBeta[q,]*stats::runif(n)
  newBeta <- rbind(newBeta,1-newBeta[q,])
  newBeta <- apply(newBeta, 2, function(col) col/sum(col))
  return(newBeta)
}

#' Construct initial values with Q groups by meging groups of a solution obtained with Q+1 groups
#'
#' @param beta_Qp1 beta for a model with Q+1 latent blocks
#' @param nbOfBeta number of initializations for the latent block memberships
#' @param nbOfPointsPerBeta number of initializations of the latent binary graph associated with each initial latent block memberships
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param model Implemented models:
#' \describe{
#'   \item{\code{Gauss}}{all Gaussian parameters of the null and the alternative distributions are unknown ; this is the Gaussian model with maximum number of unknown parameters}
#'   \item{\code{Gauss0}}{compared to \code{Gauss}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss01}}{compared to \code{Gauss}, the null distribution is set to N(0,1)}
#'   \item{\code{GaussEqVar}}{compared to \code{Gauss}, all Gaussian variances (of both the null and the alternative) are supposed to be equal, but unknown}
#'   \item{\code{Gauss0EqVar}}{compared to \code{GaussEqVar}, the mean of the null distribution is set to 0}
#'   \item{\code{Gauss0Var1}}{compared to \code{Gauss}, all Gaussian variances are set to 1 and the null distribution is set to N(0,1)}
#'   \item{\code{Gauss02distr}}{the mean of the null distribution is set to 0 and the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{Gauss2distr}}{the alternative distribution is a single Gaussian distribution, i.e. the block memberships of the nodes do not influence on the alternative distribution}
#'   \item{\code{GaussAffil}}{compared to \code{Gauss}, for the alternative distribution, there's a distribution for inter-group and another for intra-group interactions}
#'   \item{\code{Exp}}{the null and the alternatives are all exponential distributions (i.e. Gamma distributions with shape parameter equal to one)  with unknown scale parameters}
#'   \item{\code{ExpGamma}}{the null distribution is an unknown exponential, the alterantive distribution are Gamma distributions with unknown parameters}
#' }
#'
#' @return list of initial points of beta and rho of length nbOfBeta*nbOfPointsPerBeta
initialPointsByMerge <- function(beta_Qp1, nbOfBeta, dataMatrix, modelFamily, model){
  beta1 <- betaDown(beta_Qp1[[1]], nbOfBeta)
  beta2 <- betaDown(beta_Qp1[[2]], nbOfBeta)
  nbOfBeta <- min(nbOfBeta, length(beta1))
  listOfBeta <- lapply(1:nbOfBeta, function(k) list(beta1[[k]],beta2[[k]]))
  init <- initialRho(listOfBeta, dataMatrix, modelFamily, model)
  return(init)
}

#' Create new initial values by merging pairs of groups of provided beta
#'
#' Create nbOfMerges new initial values by merging nbOfMerges (or all possible) pairs of groups of provided beta
#'
#' @param beta soft node clustering
#' @param nbOfMerges number of required merges of blocks
#'
#' @return a list of length nbOfMerges (at most) of initial points for beta
betaDown <- function(beta, nbOfMerges, by = 'row'){
  n <- ncol(beta)
  Qold <- nrow(beta)
  if (Qold==2){
    listOfBeta <- list(matrix(1,1,n))
  } else {
    listOfBeta <- list()
    if (nbOfMerges>=Qold*(Qold-1)/2){ # merge all the possible pairs of groups (q,l) with q<l
      for (q in 1:(Qold-1)){
        for (l in ((q+1):Qold)){
          newBeta <- beta
          newBeta <- rbind(newBeta[-c(q,l),], newBeta[q,]+newBeta[l,])
          listOfBeta <- append(listOfBeta, list(newBeta))
        }
      }
    } else {
      # start by merging the 2 components with largest entropies
      entropy <- beta*log(beta)
      entropy[is.na(entropy)] <- 0
      entropyPerGroup <- apply(entropy,1,sum)
      q <- order(entropyPerGroup)[1]
      l <- order(entropyPerGroup)[2]
      newBeta <- beta
      newBeta <- rbind(newBeta[-c(q,l),], newBeta[q,]+newBeta[l,])
      listOfBeta <- append(listOfBeta, list(newBeta))
      # merge a number nbOfMerges-1 of pairs of groups (q,l) q<l which are randomly chosen
      groupsToMerge <- sample(1:(Qold*(Qold-1)/2), nbOfMerges-1, replace=FALSE)
      for (ind_ql in groupsToMerge){
        newBeta <- beta
        q <- convertGroupPairIdentifier(ind_ql, Qold)[1]
        l <- convertGroupPairIdentifier(ind_ql, Qold)[2]
        newBeta <- rbind(newBeta[-c(q,l),], newBeta[q,]+newBeta[l,])
        listOfBeta <- append(listOfBeta, list(newBeta))
      }
    } # end if/else nbOfMerges>=Q*(Q-1)/2
  } # end if/else Q==2
  return(listOfBeta)
}




