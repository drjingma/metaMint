#' VEM algorithm to adjust the noisy bipartite stochastic block model to an observed dense adjacency matrix
#' a number of passes can be used to deal with potential sensitivity to initialization.
#'
#' \code{fitNobiSBM()} estimates model parameters of the noisy bipartite stochastic block model and provides a clustering of the nodes
#'
#' @details
#' \code{fitNobiSBM()} supports different probability distributions for the edges and can estimate the number of node blocks
#'
#' @param dataMatrix observed dense adjacency matrix
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
#' @param exclusive_Rows_Cols whether to set the number of latent blocks in Rows to be the same as the number of latent blocks in Columns. 
#' @param sbmSize list of parameters determining the size of SBM (the number of latent blocks) to be explored. Currently restrict the number of blocks on both sides to be equal.
#' \describe{
#'   \item{\code{Q1}}{range for the number of latent blocks in Rows}
#'   \item{\code{Q2}}{range for the number of latent blocks in Columns}
#'   \item{\code{explor}}{if \code{Qmax} is not provided, then \code{Qmax} is automatically determined as \code{explor} times the number of blocks where the ICL is maximal}
#' }
#' @param filename results are saved in a file with this name (if provided)
#' @param initParam list of parameters that fix the number of initialization
#' \describe{
#'   \item{\code{nbOfbeta}}{number of initial points for the node clustering (i. e. for the variational parameters \code{beta})}
#'   \item{\code{nbOfPointsPerbeta}}{number of initial points of the latent binary graph}
#'   \item{\code{maxNbOfPasses}}{maximum number of passes through the SBM models, that is, passes from \code{Qmin} to \code{Qmax} or inversely}
#'   \item{\code{minNbOfPasses}}{minimum number of passes through the SBM models}
#' }
#' @param nbCores number of cores used for parallelization
#'
#' @return Returns a list of estimation results for all numbers of latent blocks considered by the algorithm.
#' Every element is a list composed of:
#' \describe{
#'   \item{\code{theta}}{estimated parameters of the noisy bipartite stochastic block model; a list with the following elements:
#'   \describe{
#'     \item{\code{alpha1}}{parameter estimate of alpha1}
#'     \item{\code{alpha2}}{parameter estimate of alpha2}
#'     \item{\code{pi}}{parameter estimate of pi}
#'     \item{\code{nu0}}{parameter estimate of nu0}
#'     \item{\code{nu}}{parameter estimate of nu}
#'    }}
#'   \item{\code{clustering_row}}{node clustering for rows obtained by the noisy bipartite stochastic block model, more precisely, a hard clustering given by the
#'   maximum a posterior estimate of the variational parameters \code{sbmParam$edgeProba[[1]]}}
#'   \item{\code{clustering_col}}{node clustering for columns obtained by the noisy bipartite stochastic block model, more precisely, a hard clustering given by the
#'   maximum a posterior estimate of the variational parameters \code{sbmParam$edgeProba[[2]]}}
#'   \item{\code{sbmParam}}{further results concerning the latent binary stochastic block model. A list with the following elements:
#'   \describe{
#'     \item{\code{Q1,Q2}}{number of latent blocks in the noisy bipartite stochastic block model}
#'     \item{\code{clusterProba}}{a list of length two, consisting of soft clustering given by the conditional probabilities of a node to belong to a given latent block.
#'     In other words, these are the variational parameters \code{beta1}; (Q1 x n1)-matrix and \code{beta2}; (Q2 x n2)-matrix}
#'     \item{\code{edgeProba}}{conditional probabilities \code{rho} of an edges given the node memberships of the interacting nodes; (N_Q x N)-matrix}
#'     \item{\code{ICL}}{value of the ICL criterion at the end of the algorithm}
#'     }}
#'   \item{\code{convergence}}{a list of convergence indicators:
#'   \describe{
#'     \item{\code{J}}{value of the lower bound of the log-likelihood function at the end of the algorithm}
#'     \item{\code{complLogLik}}{value of the complete log-likelihood function at the end of the algorithm}
#'     \item{\code{converged}}{indicates if algorithm has converged}
#'     \item{\code{nbIter}}{number of iterations performed}
#'  }}
#' }
#'
#' @export
#' @examples
#' n1 <- 10
#' n2 <- 15
#' Q1 <- 2
#' Q2 <- 2
#' theta <- list(alpha1 = c(0.5, 0.5), 
#'               alpha2 = c(0.4, 0.6),
#'               nu0 = c(0,.1),
#'               nu = list(mean = matrix(2, Q1, Q2), 
#'                         sd = matrix(0.1, Q1, Q2)),  
#'               pi=matrix(0.01,Q1,Q2) )
#' diag(theta$pi) <- 0.3
#' obs <- rnbisbm(n1, n2, theta, modelFamily='Gauss')
#' res <- fitNobiSBM(obs$dataMatrix, sbmSize = list(Q1 =1:3, Q2=1:3),
#'        initParam=list(nbOfbeta=1, nbOfPointsPerbeta=1), nbCores=1)
fitNobiSBM <- function(dataMatrix, model='Gauss0', 
                       sbmSize = list(Q1 = 1:5, Q2 = 1:5, explor=1.5),
                       exclusive_Rows_Cols = TRUE, 
                       filename=NULL,
                       initParam = list(nbOfbeta=NULL, nbOfPointsPerbeta=NULL,
                                        maxNbOfPasses=NULL, minNbOfPasses=1),
                       nbCores=parallel::detectCores()){
  
  if (!exclusive_Rows_Cols){
    # there exist rows or columns that belong to two clusters
    bestSolutionAtQ <- fitNobiSBM_unequal(dataMatrix, model, sbmSize,filename, nbCores)
  } else {
    # Exclusive rows and columns
    bestSolutionAtQ <- fitNobiSBM_equal(dataMatrix, model, sbmSize, filename, initParam, nbCores)
    
  }
  return(bestSolutionAtQ)
}

fitNobiSBM_unequal <- function(dataMatrix, model='Gauss0',
                               sbmSize = list(Q1 = 1:5, Q2 = 1:5),
                               filename=NULL,
                               nbCores=parallel::detectCores()){
  
  Q1min <- max(c(1, min(sbmSize$Q1)))
  Q1max <- max(sbmSize$Q1)
  
  Q2min <- max(c(1, min(sbmSize$Q2)))
  Q2max <- max(sbmSize$Q2)
  
  possibleQvalues <- expand.grid(Q1min:Q1max,Q2min:Q2max)
  
  if (Sys.info()["sysname"]=="Windows"){
    nbCores <- 1
  } else {
    nbCores <- if (nbCores>1) min(c(round(nbCores), parallel::detectCores())) else 1
  }
  
  doParallelComputing <- (nbCores>1)
  
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  
  if (model %in% c('Gauss','Gauss0','Gauss01','GaussEqVar','Gauss0EqVar','Gauss0Var1','Gauss02distr','Gauss2distr','GaussAffil'))
    modelFamily <- 'Gauss'
  if (model %in% c('Exp','ExpGamma','Gamma'))
    modelFamily <- 'Gamma'
  
  M <- nrow(possibleQvalues)
  bestSolutionAtQ <- lapply(1:M, function(elem) list(convergence=list(J=-Inf))) # list of best solutions for every combinations of Q1 and Q2
  
  init <- lapply(1:M, function(currentQindice){
    Q1 <- possibleQvalues[currentQindice,1]
    Q2 <- possibleQvalues[currentQindice,2]
    init <- initialPoints(Q1, Q2, dataMatrix, modelFamily, model, method='kmeans')
    list(beta = init$beta[[1]], rho = init$rho[[1]])
  })
  ListOfbetaRho <- list(beta = lapply(init, function(a) a$beta),
                        rho = lapply(init, function(a) a$rho))
  
  ## launch VEM for those initial points
  if(doParallelComputing){
    bestSolutionAtQ <- parallel::mclapply(1:M, function(k){
      mainVEM_Q_par(k, ListOfbetaRho, modelFamily, model, dataMatrix)},
      mc.cores=nbCores)
  }else{
    for (s in 1:M){
      # the number of runs is determined by the number of grid points
      cat('Q1 =', possibleQvalues[s,1], 'Q2 =', possibleQvalues[s,2], 'Q1max =', Q1max, 'Q2max =', Q2max, '\n')
      currentInitialPoint <- if (model != 'ExpGamma') list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]]) else
        list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]], nu=ListOfbetaRho$nu[[s]])
      currentSolution <- mainVEM_Q(currentInitialPoint, modelFamily, model, dataMatrix)
      
      if (currentSolution$convergence$J>bestSolutionAtQ[[s]]$convergence$J)
        bestSolutionAtQ[[s]] <- currentSolution
    }
  }
  
  if (!is.null(filename)) save(bestSolutionAtQ, file=filename)
  
  return(bestSolutionAtQ)
}


fitNobiSBM_equal <- function(dataMatrix, model='Gauss0', 
                             sbmSize = list(Q1 = 1:5, Q2 = 1:5, explor=1.5),
                             filename=NULL,
                             initParam = list(nbOfbeta=NULL, nbOfPointsPerbeta=NULL,
                                              maxNbOfPasses=NULL, minNbOfPasses=1),
                             nbCores=parallel::detectCores()){
  
  Qmin <- max(c(1, min(sbmSize$Q1)))
  explor <- if(is.null(sbmSize$explore)) 1.5 else sbmSize$explor
  Qmax <- max(sbmSize$Q1)
  QmaxIsFlexible <- FALSE
  
  minNbOfPasses <- max(c(1,initParam$minNbOfPasses))
  
  maxNbOfPasses <- if (is.null(initParam$maxNbOfPasses)) max(c(10,initParam$minNbOfPasses)) else initParam$maxNbOfPasses
  
  if (Sys.info()["sysname"]=="Windows"){
    nbCores <- 1
  } else {
    nbCores <- if (nbCores>1) min(c(round(nbCores), parallel::detectCores())) else 1
  }
  
  doParallelComputing <- (nbCores>1)
  
  n1 <- nrow(dataMatrix)
  n2 <- ncol(dataMatrix)
  
  if (model %in% c('Gauss','Gauss0','Gauss01','GaussEqVar','Gauss0EqVar','Gauss0Var1','Gauss02distr','Gauss2distr','GaussAffil'))
    modelFamily <- 'Gauss'
  if (model %in% c('Exp','ExpGamma','Gamma'))
    modelFamily <- 'Gamma'
  
  possibleQvalues <- Qmin:Qmax
  bestSolutionAtQ <- lapply(1:(Qmax-Qmin+1), function(elem) list(convergence=list(J=-Inf))) # list of best solutions for every Q
  
  currentQindice <- 1
  notConverged <- TRUE
  up <- TRUE  # redundant since up=TRUE iff nbOfPasses is odd
  nbOfPasses <- 1
  while (notConverged){
    Q1 <- Q2 <- Q <- possibleQvalues[currentQindice]
    N_Q <- Q1*Q2
    cat("-- pass =", nbOfPasses, 'Q1 = Q2 =', Q, 'Qmax =', Qmax, '\n')
    
    ## create list of initial points
    if (nbOfPasses==1){
      ListOfbetaRho <- initialPoints(Q1, Q2, dataMatrix, modelFamily, model, method='kmeans')   # use kmeans, random initializations for gamma parameters, for gamma use more kmeans perturbations??
      if (Q>Qmin){ # add split up solution if Q>Qmin
        additionalbetaRho <- initialPointsBySplit(bestSolutionAtQ[[currentQindice-1]]$sbmParam$clusterProba, nbOfBeta = round(Q/2), dataMatrix, modelFamily, model)
        ListOfbetaRho$beta <- append(ListOfbetaRho$beta, additionalbetaRho$beta)
        ListOfbetaRho$rho <- append(ListOfbetaRho$rho, additionalbetaRho$rho)
        if(model=='ExpGamma'){
          ListOfbetaRho$nu <- append(ListOfbetaRho$nu, additionalbetaRho$nu)
        }
      }
    }else{ # nbOfPasses>1
      if (up){ # add split up solution if Q>Qmin
        ListOfbetaRho <- initialPointsBySplit(bestSolutionAtQ[[currentQindice-1]]$sbmParam$clusterProba, nbOfBeta= max(c(3,round(Q/2))), dataMatrix, modelFamily, model)
      }else{
        # merge solutions with betaDown
        ListOfbetaRho <- initialPointsByMerge(bestSolutionAtQ[[currentQindice+1]]$sbmParam$clusterProba, nbOfBeta= max(c(5,round(Q/2))), dataMatrix, modelFamily, model)
      }
    }
    
    ## launch VEM for those initial points
    M <- length(ListOfbetaRho$beta)
    # M
    if(doParallelComputing){
      ListOfSolutions <- parallel::mclapply(1:M, function(k){
        mainVEM_Q_par(k, ListOfbetaRho, modelFamily, model, dataMatrix)},
        mc.cores=nbCores)
      ListOfJ <- lapply(ListOfSolutions, function(solutionThisRun) solutionThisRun$convergence$J)
      currentSolution <- ListOfSolutions[[which.max(ListOfJ)]]
      if (currentSolution$convergence$J>bestSolutionAtQ[[currentQindice]]$convergence$J)
        bestSolutionAtQ[[currentQindice]] <- currentSolution
    }else{
      for (s in 1:M){
        # the number of runs is determined by the number of initial points
        cat("-- pass = ", nbOfPasses, "Q1 = Q2 = ",Q,'run = ',s,"\n")
        currentInitialPoint <- if (model != 'ExpGamma') list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]]) else
          list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]], nu=ListOfbetaRho$nu[[s]])
        print(dim(currentInitialPoint$rho))
        currentSolution <- mainVEM_Q(currentInitialPoint, modelFamily, model, dataMatrix)
        print(currentSolution$convergence$J)
        if (currentSolution$convergence$J>bestSolutionAtQ[[currentQindice]]$convergence$J)
          bestSolutionAtQ[[currentQindice]] <- currentSolution
      }
    }
    
    ## evaluate results of the current pass and decide whether to continue or not
    if (nbOfPasses==1){
      if (Q<Qmax){
        currentQindice <- currentQindice + 1
      }else{ #Q==Qmax
        if (!QmaxIsFlexible){ # fixed Qmax
          if ((Qmax==Qmin)){
            notConverged <- FALSE
          }else{ # fixed Qmax, Qmax>Qmin -> prepare 2nd pass
            bestQ <- getBestQ(bestSolutionAtQ)
            nbOfPasses <- 2
            currentQindice <- currentQindice - 1
            up <- FALSE
          }
        }else{ #nbOfPasses=1, flexible Qmax
          bestQ <- getBestQ(bestSolutionAtQ)  # list $ICL, $Q
          if(bestQ$Q<=round(Qmax/explor)){ # if max is attained for a 'small' Q -> go to 2nd pass
            nbOfPasses <- 2
            currentQindice <- currentQindice - 1
            up <- FALSE
          }else{ # increase Qmax and continue first pass
            oldQmax <- Qmax
            Qmax <- round(bestQ$Q*explor)
            possibleQvalues <- Qmin:Qmax
            bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(convergence=list(J=-Inf))))
            # bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(J=-Inf)))
            currentQindice <- currentQindice + 1
          }
        }
      }
    }else{ # nbOfPasses>1
      if((Q>Qmin)&(Q<Qmax)){
        currentQindice <- currentQindice + up - !up
      }else{ # at the end of a pass
        bestQInThisPass <- getBestQ(bestSolutionAtQ)  # list $ICL, $Q
        if(QmaxIsFlexible){
          if(bestQInThisPass$ICL>bestQ$ICL){ # in this pass the solution was improved
            if(bestQInThisPass$Q<=round(Qmax/explor) ){ # the improvement is achieved for a "small" Q -> go to next pass
              if ((Q==Qmin)){ # only stop when Q==Qmin
                bestQ <- bestQInThisPass
                notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
              }
              nbOfPasses <- nbOfPasses + 1
              currentQindice <- currentQindice - up + !up
              up <- !up
            }else{ # improvement is achieved for a "large" Q -> increase Qmax and continue pass
              oldQmax <- Qmax
              Qmax <- round(bestQ$Q*explor)
              possibleQvalues <- Qmin:Qmax
              bestSolutionAtQ <- append(bestSolutionAtQ, lapply((oldQmax+1):Qmax, function(elem) list(convergence=list(J=-Inf))))
              currentQindice <- currentQindice + 1
              if (!up){ # go to next pass
                bestQ <- bestQInThisPass
                notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
                nbOfPasses <- nbOfPasses + 1
                up <- TRUE
              }
            }
          }else{ # no improvement -> exit
            notConverged <- (nbOfPasses<minNbOfPasses) #FALSE
          }
        }else{ # fixed Qmax
          if (bestQInThisPass$ICL>bestQ$ICL){ # in this pass the solution was improved -> prepare next pass
            if (Q==Qmin){ # only stop when Q==Qmin
              bestQ <- bestQInThisPass
              notConverged <- (nbOfPasses<maxNbOfPasses) # nb of max passes attained ?
            }
            nbOfPasses <- nbOfPasses + 1
            currentQindice <- currentQindice - up + !up
            up <- !up
          }else{ # no improvement
            if (nbOfPasses>=minNbOfPasses) # min nb of passes attained -> exit
              notConverged <- FALSE
            else{ # continue
              nbOfPasses <- nbOfPasses + 1
              currentQindice <- currentQindice - up + !up
              up <- !up
            }
          }
        }
      }
    }
  }
  if (!is.null(filename)) save(bestSolutionAtQ, file=filename)
  
  return(bestSolutionAtQ)
}



#' main function of VEM algorithm with fixed number of SBM blocks
#'
#' @param init list of initial points for the algorithm
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
#' @param dataMatrix observed dense adjacency matrix
#'
#' @return list of estimated model parameters and a node clustering; like the output of fitNobiSBM()
mainVEM_Q <- function(init, modelFamily, model, dataMatrix){
  
  nb.iter <- 500
  epsilon <- 1e-6
  
  Q1 <- nrow(init$beta[[1]])
  Q2 <- nrow(init$beta[[2]])
  N_Q <- Q1*Q2
  n1 <- ncol(init$beta[[1]])
  n2 <- ncol(init$beta[[2]])
  
  VE <- list(beta1 = init$beta[[1]], beta2 = init$beta[[2]], rho=init$rho)
  if (modelFamily=='Gauss'){
    theta.init <- list(pi=matrix(NA, Q1, Q2),
                       nu0=c(0,1),
                       nu=list(mean=matrix(0,Q1,Q2),
                               sd = matrix(1,Q1,Q2)))
  }else{
    if(model=='Exp'){
      theta.init <- list(pi=rep(NA, N_Q),
                         nu0=c(1,1),
                         nu=matrix(c(1,1), N_Q, 2, byrow=TRUE))
    }
    if(model=='ExpGamma'){
      theta.init <- list(pi=rep(NA, N_Q),
                         nu0=c(1,1),
                         nu=init$nu)
    }
  }
  theta.init$alpha1 <- if (Q1>1) rowMeans(VE$beta1) else 1
  theta.init$alpha2 <- if (Q2>1) rowMeans(VE$beta2) else 1
  mstep <- theta.init
  J.old <- -Inf
  Jeval <- list(J=-Inf, complLogLik=-Inf)
  
  convergence <- list(converged=FALSE, nb.iter.achieved=FALSE)
  it <- 0
  while (sum(convergence==TRUE)==0){
    # cat('current iteration ...', it, '...\n')
    it <- it+1
    mstep <- Mstep(VE, mstep, model, dataMatrix, modelFamily)
    VE <- VEstep(VE, mstep, dataMatrix, modelFamily)
    Jeval <- JEvalMstep(VE, mstep, dataMatrix, modelFamily)
    convergence$nb.iter.achieved <- (it > nb.iter+1)
    convergence$converged <-  (abs((Jeval$J-J.old)/Jeval$J)< epsilon)
    J.old <- Jeval$J
  }
  
  solutionThisRun <- list(theta=list(pi=mstep$pi, 
                                     nu0=mstep$nu0, 
                                     nu=mstep$nu, 
                                     alpha1 = mstep$alpha1, 
                                     alpha2 = mstep$alpha2),
                          clustering_row=NULL, clustering_col=NULL,
                          sbmParam=list(Q1=Q1, Q2=Q2, 
                                        clusterProba=list(VE$beta1,VE$beta2), 
                                        edgeProba=VE$rho, ICL=NULL),
                          convergence = list(J=Jeval$J, complLogLik=Jeval$complLogLik,
                                             converged=convergence$converged, nbIter=it))
  solutionThisRun$sbmParam$ICL <- ICL_Q(solutionThisRun, model)
  solutionThisRun$clustering_row <- if (Q1>1) apply(VE$beta1,2,which.max) else rep(1, n1)
  solutionThisRun$clustering_col <- if (Q2>1) apply(VE$beta2,2,which.max) else rep(1, n2)
  
  return(solutionThisRun)
}



#' main function of VEM algorithm for fixed number of latent blocks in parallel computing
#'
#' runs the VEM algorithm the provided initial point
#'
#' @param s indice of initial point in ListOfbetaRho to be used for this run
#' @param ListOfbetaRho a list of initial points
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
#' @param dataMatrix observed dense adjacency matrix
#'
#' @return list of estimated model parameters and a node clustering; like the output of fitNobiSBM()
mainVEM_Q_par <- function(s, ListOfbetaRho, modelFamily, model, dataMatrix){
  currentInitialPoint <- if (model != 'ExpGamma') list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]]) else
    list(beta=ListOfbetaRho$beta[[s]], rho=ListOfbetaRho$rho[[s]], nu=ListOfbetaRho$nu[[s]])
  currentSolution <- mainVEM_Q(currentInitialPoint, modelFamily, model, dataMatrix)
  return(currentSolution)
}


#' VE-step (updated for biSBM)
#'
#' performs one VE-step, that is, update of beta and rho
#'
#' @param VE list with variational parameters beta and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#' @param fix.iter maximal number of iterations for fixed point equation
#'
#' @return updated list \code{VE} with variational parameters beta and rho
VEstep <- function(VE, mstep, dataMatrix, modelFamily, fix.iter=5){
  Q1 <- nrow(VE$beta1)
  Q2 <- nrow(VE$beta2)
  epsilon <- 1e-6
  # compute rho using pi
  ind <- 0
  for (q in 1:Q1){
    for (l in 1:Q2){
      rho_numerateur <- modelDensity(dataMatrix, c(mstep$nu$mean[q,l],
                                                   mstep$nu$sd[q,l]), modelFamily) * mstep$pi[q,l]  # taille N
      VE$rho[q,l,,] <-  rho_numerateur / (rho_numerateur + modelDensity(dataMatrix, mstep$nu0, modelFamily)*(1-mstep$pi[q,l]))
    }
  }
  
  log.pi <- log(mstep$pi)
  log.pi[mstep$pi==0] <- 0
  log.1mpi <- log(1-mstep$pi)
  log.1mpi[mstep$pi==1] <- 0
  
  # solve the fixed point equation
  converged <- FALSE
  it <- 0
  while ((!converged)&(it<fix.iter)){
    it <- it+1
    beta.old <- list(VE$beta1,VE$beta2) # a list of length 2
    fit <- betaUpdate(beta.old, log.pi, log.1mpi, dataMatrix, VE, mstep, modelFamily)
    VE$beta1 <- fit$beta1
    VE$beta2 <- fit$beta2
    converged <- ((max(abs(c(VE$beta1,VE$beta2)-unlist(beta.old)))<epsilon))
  }
  
  return(VE)
}

#' Compute one iteration to solve the fixed point equation in the VE-step
#'
#' @param beta.list current value of beta in a list of two
#' @param log.pi value of log(w)
#' @param log.1mpi value of log(1-w)
#' @param dataMatrix observed dense adjacency matrix
#' @param VE list with variational parameters beta and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return updated value of beta
betaUpdate <- function (beta.list, log.pi, log.1mpi, dataMatrix, VE, mstep, modelFamily){
  beta1 <- beta.list[[1]]
  beta2 <- beta.list[[2]]
  
  Q1 <- nrow(beta1)
  n1 <- ncol(beta1)
  Q2 <- nrow(beta2)
  n2 <- ncol(beta2)
  
  beta1.new <- beta1
  beta2.new <- beta2
  
  rho.mat <- VE$rho
  log.rho.mat <- log(rho.mat)
  log.rho.mat[rho.mat==0] <- 0
  log.1mrho.mat <- log(1-rho.mat)
  log.1mrho.mat[rho.mat==1] <- 0
  
  # first compute the D matrix which is Q1 by Q2 by n1 by n2
  D.mat <- array(0, c(Q1,Q2,n1,n2))
  for (q in 1:Q1){
    for (l in 1:Q2){
      D.mat[q, l, , ] = rho.mat[q,l,,] * (log.pi[q,l] + log(modelDensity(dataMatrix, c(mstep$nu$mean[q,l],mstep$nu$sd[q,l]), modelFamily)) - log.rho.mat[q,l,,]) +
        (1-rho.mat[q,l,,])*( log.1mpi[q,l] +  log(modelDensity(dataMatrix, mstep$nu0, modelFamily)) - log.1mrho.mat[q,l, ,] )
    }
  }
  
  if (Q1>1){
    logbeta1 <- log(beta1)
    
    for (i in 1:n1){
      for (q in 1:Q1){
        logbeta1[q,i] <- log(mstep$alpha1[q]) + sum(beta2 * D.mat[q, , i, ])
      }
      
      ## Normalizing in the log space to avoid numerical problems
      ## and going back to exponential with the same normalization
      logbeta1[,i] <- logbeta1[,i]-max(logbeta1[,i])
      beta1.new[,i] <- exp(logbeta1[,i])
      beta1.new[,i] <- correctBeta(beta1.new[,i])
    }
  }
  
  if (Q2>1){
    logbeta2 <- log(beta2)
    
    for (j in 1:n2){
      for (l in 1:Q2){
        logbeta2[l,j] <- log(mstep$alpha2[l]) + sum(beta1 * D.mat[, l, , j])
      }
      
      ## Normalizing in the log space to avoid numerical problems
      ## and going back to exponential with the same normalization
      logbeta2[,j] <- logbeta2[,j]-max(logbeta2[,j])
      beta2.new[,j] <- exp(logbeta2[,j])
      beta2.new[,j] <- correctBeta(beta2.new[,j])
    }
  }
  return(list(beta1 = beta1.new, beta2 = beta2.new))
}


#' M-step
#'
#' performs one M-step, that is, update of pi, alpha1, alpha2, nu, nu0
#'
#' @param VE list with variational parameters beta and rho
#' @param mstep list with current model parameters and additional auxiliary terms
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
#' @param dataMatrix observed dense adjacency matrix 
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return updated list \code{mstep} with current model parameters and additional auxiliary terms
Mstep <- function(VE, mstep, model, dataMatrix, modelFamily){
  Q1 <- nrow(VE$beta1)
  n1 <- ncol(VE$beta1)
  Q2 <- nrow(VE$beta2)
  n2 <- ncol(VE$beta2)
  N <- n1*n2
  
  mstep$alpha1 <- if (Q1>1) rowMeans(VE$beta1) else 1
  mstep$alpha2 <- if (Q2>1) rowMeans(VE$beta2) else 1
  mstep$sum_beta <- matrix(0,Q1,Q2)
  mstep$sum_kappa <- matrix(0,Q1,Q2)
  
  denom_1 <- denom_2 <- denom_3 <- numer_1 <- numer_2 <- numer_3 <- numer_4 <- numer_5 <- numer_6 <- 0
  default_var <- sqrt(log(min(n1,n2))) # default variance for each block
  for (q in 1:Q1){
    for (l in 1:Q2){
      beta_ql <- getBetaql(q, l, list(VE$beta1,VE$beta2))
      kappa_ql <- beta_ql*VE$rho[q,l, ,] 
      
      s <- sum(kappa_ql) 
      mstep$sum_beta[q,l] <- sum(beta_ql)
      mstep$sum_kappa[q,l] <- s
      
      calcul_pi <- s/mstep$sum_beta[q,l]
      mstep$pi[q,l] <- if (calcul_pi <=.Machine$double.eps) .Machine$double.eps else calcul_pi
      mstep$pi[q,l] <- if (calcul_pi >= 1-.Machine$double.eps) 1-.Machine$double.eps else calcul_pi
      
      # mstep$pi[q,l]
      if (modelFamily=='Gauss'){
        if (model=='Gauss'){
          if (s > 0){
            mstep$nu$mean[q,l] <- sum(kappa_ql*dataMatrix)/s
            sigma <- sqrt(sum(kappa_ql*(dataMatrix-mstep$nu$mean[q,l])^2)/s)
            mstep$nu$sd[q,l] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu$mean[q,l] <- 0
            mstep$nu$sd[q,l] <- default_var
          }
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          
          # the following are quantities needed to calculate nu0
          numer_1 <- numer_1 + sum(beta_ql_1mrho*dataMatrix) # mean
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2) # variance
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
        }
        
        if(model=='Gauss0'){
          if (s > 0){
            mstep$nu$mean[q,l] <- sum(kappa_ql*dataMatrix)/s
            sigma <- sqrt(sum(kappa_ql*(dataMatrix-mstep$nu$mean[q,l])^2)/s)
            mstep$nu$sd[q,l] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu$mean[q,l] <- 0
            mstep$nu$sd[q,l] <- default_var
          }
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2)
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
        }
        
        if(model=='Gauss01'){
          if (s > 0){
            mstep$nu$mean[q,l] <- sum(kappa_ql*dataMatrix)/s
            sigma <- sqrt(sum(kappa_ql*(dataMatrix-mstep$nu$mean[q,l])^2)/s)
            mstep$nu$sd[q,l] <- if (sigma <.Machine$double.eps) default_var else sigma
          }else{ # s==0
            mstep$nu$mean[q,l] <- 0
            mstep$nu$sd[q,l] <- default_var
          }
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
        }
        
        if (model=='GaussEqVar'){
          mstep$nu$mean[q,l] <- if (s>0) sum(kappa_ql*dataMatrix)/s else 0
          numer_4 <- numer_4 + sum(kappa_ql*(dataMatrix-mstep$nu$mean[q,l])^2)
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
          numer_1 <- numer_1 + sum(beta_ql_1mrho*dataMatrix)
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2)
        }
        
        if (model=='Gauss0EqVar'){
          mstep$nu$mean[q,l] <- if (s>0) sum(kappa_ql*dataMatrix)/s else 0
          numer_2 <- numer_2 + sum(kappa_ql*(dataMatrix-mstep$nu$mean[q,l])^2)
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2)
        }
        
        if (model=='Gauss0Var1'){
          mstep$nu$mean[q,l] <- if (s>0) sum(kappa_ql*dataMatrix)/s else 0
        }
        
        if (model=='Gauss02distr'){
          numer_3 <- numer_3 + sum(kappa_ql*dataMatrix)  # for mean estimate mu_1
          numer_4 <- numer_4 + sum(kappa_ql*dataMatrix^2) # for variance estimate sigma_1
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2) # for variance estimate sigma_0
        }
        
        if (model=='Gauss2distr'){
          numer_3 <- numer_3 + sum(kappa_ql*dataMatrix)  # for mean estimate mu_1
          numer_4 <- numer_4 + sum(kappa_ql*dataMatrix^2) # for variance estimate sigma_1
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
          numer_1 <- numer_1 + sum(beta_ql_1mrho*dataMatrix) # for mean estimate mu_0
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2) # for variance estimate sigma_0
        }
        
        if (model=='GaussAffil'){
          if (q==l){
            numer_3 <- numer_3 + sum(kappa_ql*dataMatrix)  # for mean estimate mu_in
            numer_4 <- numer_4 + sum(kappa_ql*dataMatrix^2) # for variance estimate sigma_in
            denom_2 <- denom_2 + sum(kappa_ql)
          }else{
            numer_5 <- numer_5 + sum(kappa_ql*dataMatrix)  # for mean estimate mu_out
            numer_6 <- numer_6 + sum(kappa_ql*dataMatrix^2) # for variance estimate sigma_out
            denom_3 <- denom_3 + sum(kappa_ql)
          }
          beta_ql_1mrho <- beta_ql*(1-VE$rho[q,l,,])
          denom_1 <- denom_1 + sum(beta_ql_1mrho)
          numer_1 <- numer_1 + sum(beta_ql_1mrho*dataMatrix) # for mean estimate mu_0
          numer_2 <- numer_2 + sum(beta_ql_1mrho*dataMatrix^2) # for variance estimate sigma_0
        }
      }
    } # end for l
  }# end for q
  if (modelFamily=='Gauss'){
    if (model=='GaussEqVar'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
      }
      sigma <- if (denom_1>0) sqrt((denom_1*(numer_2/denom_1-mstep$nu0[1]^2) + numer_4)/N) else default_var
      mstep$nu0[2] <- sigma
      mstep$nu$sd <- matrix(sigma,Q1,Q2)
    }
    
    if (model=='Gauss0EqVar'){
      mstep$nu0[2] <- sqrt(numer_2/N)
      mstep$nu$sd <- matrix(mstep$nu0[2],Q1,Q2) 
    }
    
    if (model=='Gauss02distr'){
      mstep$nu0[2] <- if (denom_1>0) sqrt(numer_2/denom_1) else default_var
      
      if (denom_1<N){
        nu.mean <- numer_3/(N-denom_1)
        mstep$nu$mean <- matrix(nu.mean,Q1,Q2)
        nu.sd <- sqrt(numer_4/(N-denom_1) - nu.mean^2)
        mstep$nu$sd <- matrix(nu.sd,Q1,Q2)
      }else{
        mstep$nu$mean <- matrix(0,Q1,Q2)
        mstep$nu$sd <- matrix(default_var,Q1,Q2)
      }
    }
    
    if (model=='Gauss2distr'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
        mstep$nu0[2] <- sqrt(numer_2/denom_1 - mstep$nu0[1]^2)
      }else{
        mstep$nu0 <- c(0,default_var)
      }
      if (denom_1<N){
        nu.mean <- numer_3/(N-denom_1)
        mstep$nu$mean <- matrix(nu.mean,Q1,Q2)
        nu.sd <- sqrt(numer_4/(N-denom_1) - nu.mean^2)
        mstep$nu$sd <- matrix(nu.sd,Q1,Q2)
      }else{
        mstep$nu$mean <- matrix(0,Q1,Q2)
        mstep$nu$sd <- matrix(default_var,Q1,Q2)
      }
    }
    if (model=='Gauss0'){
      mstep$nu0[2] <- if (denom_1>0) sqrt(numer_2/denom_1) else default_var
    }
    
    if (model=='Gauss'){
      if (denom_1>0){
        mstep$nu0[1] <- numer_1/denom_1
        mstep$nu0[2] <- sqrt(numer_2/denom_1 - mstep$nu0[1]^2)
      }else{
        mstep$nu0 <- c(0,default_var)
      }
    }
    
    if ((mstep$nu0[2]<1e-7)|(is.na(mstep$nu0[2])) )
    {print(paste("nu0=", mstep$nu0[2])) ;
      mstep$nu0[2] <- default_var
    }
    varZero <- ((mstep$nu$sd < 1e-7) | (is.na(mstep$nu$sd)))
    if (sum(varZero) > 0)
      mstep$nu$sd[varZero] <- default_var
  }
  
  return(mstep)
}



#' Perform one iteration of the Newton-Raphson to compute the MLE of the parameters of the Gamma distribution
#'
#' @param param current parameters of the Gamma distribution
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#'
#' @return updated parameters of the Gamma distribution
# update_newton_gamma <- function(param, L, M){
#   gradient <- c(log(param[2])-digamma(param[1])+L, param[1]/param[2]-M)
#   w <- trigamma(param[1])
#   detH <- 1/(1-param[1]*w)
#   a.new <- min(c(param[1] - detH*sum(param*gradient), 150))
#   b.new <- param[2] - detH*sum(c(param[2],param[2]^2*w)*gradient)
#   return(c(a.new,b.new))
# }

#' evaluate the objective in the Gamma model
#'
#' @param param parameters of the Gamma distribution
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#'
#' @return value of the lower bound of the log-likelihood function
# J.gamma <- function(param, L, M){
#   return(param[1]*log(param[2]) - log(gamma(param[1])) + (param[1]-1)*L - param[2]*M)
# }


#' compute the MLE in the Gamma model using the Newton-Raphson method
#'
#' @param L weighted mean of log(data)
#' @param M weighted mean of the data
#' @param param.old parameters of the Gamma distribution
#' @param epsilon threshold for the stopping criterion
#' @param nb.iter.max maximum number of iterations
#'
#' @return updated parameters of the Gamma distribution
# emv_gamma <- function(L, M, param.old, epsilon=1e-3, nb.iter.max=10){
#   logvrais.old <- J.gamma(param.old, L, M)
#   notConverged <- TRUE
#   nb.iter <- 0
#   while ((notConverged) & (nb.iter <= nb.iter.max)){
#     param.new <- update_newton_gamma(param.old, L, M)
#     if (sum(param.new>0)==2){
#       # check convergence criterion
#       logvrais.new <- J.gamma(param.new, L, M)
#       notConverged <- (abs((logvrais.new-logvrais.old)/logvrais.new)>epsilon)
#       param.old <- param.new
#       logvrais.old <- logvrais.new
#     }
#     else{
#       notConverged <- FALSE
#     }
#     nb.iter <- nb.iter + 1
#   }
#   return(param.old)
# }


#' evaluation of the objective in the Gauss model
#' this evaluates both the ELBO and the complete data log-likelihood
#'
#' @param VE list with variational parameters beta and rho
#' @param mstep list with current model parameters and additional auxiliary terms
#' @param dataMatrix observed dense adjacency matrix
#' @param modelFamily probability distribution for the edges. Possible values:
#'       \code{Gauss}, \code{Gamma}
#'
#' @return value of the ELBO and the complete log likelihood function
JEvalMstep  <- function(VE, mstep, dataMatrix, modelFamily){
  Q1 <- nrow(VE$beta1)
  n1 <- ncol(VE$beta1)
  Q2 <- nrow(VE$beta2)
  n2 <- ncol(VE$beta2)
  
  log.beta1 <- log(VE$beta1)
  log.beta1[VE$beta1 ==0] <- 0
  log.beta2 <- log(VE$beta2)
  log.beta2[VE$beta2 ==0] <- 0
  
  log.alpha1 <- log(mstep$alpha1)
  log.alpha1[mstep$alpha1 ==0] <- 0
  log.alpha2 <- log(mstep$alpha2)
  log.alpha2[mstep$alpha2 ==0] <- 0
  
  log.rho <- log(VE$rho)
  log.rho[VE$rho==0] <- 0
  log.1mrho <- log(1-VE$rho)
  log.1mrho[VE$rho==1] <- 0
  
  log.pi <- log(mstep$pi)
  log.pi[mstep$pi==0] <- 0
  log.1mpi <- log(1-mstep$pi)
  log.1mpi[mstep$pi==1] <- 0
  
  term1 <- sum(log.1mpi*mstep$sum_beta + (log.pi-log.1mpi)*mstep$sum_kappa)
  term2 <- term3 <- matrix(NA, Q1, Q2)
  
  for (q in 1:Q1){
    for (l in 1:Q2){
      beta_ql <- getBetaql(q,l, list(VE$beta1, VE$beta2))
      
      logf1 <- log(modelDensity(dataMatrix, c(mstep$nu$mean[q,l],
                                              mstep$nu$sd[q,l]), modelFamily))
      logf0 <- log(modelDensity(dataMatrix, mstep$nu0, modelFamily))
      
      term2[q,l] <- sum(beta_ql*VE$rho[q,l,,] * logf1 + beta_ql*(1-VE$rho[q,l,,])*logf0)
      
      term3[q,l] <- sum(beta_ql*(VE$rho[q,l,,] * log.rho[q,l,,] + (1-VE$rho[q,l,,])*log.1mrho[q,l, ,]))
    }
  }
  term2_sum <- sum(term2)
  term3_sum <- sum(term3)
  
  J <- sum( log.alpha1 %*% VE$beta1 )- sum( VE$beta1 * log.beta1 ) +
    sum( log.alpha2 %*% VE$beta2 )- sum( VE$beta2 * log.beta2 ) + term1 + term2_sum - term3_sum
  
  complLogLik <- sum( log.alpha1 %*% VE$beta1 ) + sum( log.alpha2 %*% VE$beta2 ) + term1 + term2_sum
  
  return(list(J=J, complLogLik=complLogLik))
}



#' computation of the Integrated Classification Likelihood criterion
#'
#' computation of the Integrated Classification Likelihood criterion for a result provided by mainVEM_Q()
#'
#' @param solutionThisRun result provided by mainVEM_Q()
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
#' @return value of the ICL criterion
ICL_Q <- function(solutionThisRun, model){
  Q1 <- solutionThisRun$sbmParam$Q1
  Q2 <- solutionThisRun$sbmParam$Q2
  n1 <- ncol(solutionThisRun$sbmParam$clusterProba[[1]])
  n2 <- ncol(solutionThisRun$sbmParam$clusterProba[[2]])
  N_Q <- Q1*Q2
  N <- n1*n2 
  dimPi <- N_Q
  if (model %in% c('Gauss0', 'ExpGamma')){
    dimH0 <- 1
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss'){
    dimH0 <- 2
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss01'){
    dimH0 <- 0
    dimH1 <- 2*N_Q
  }
  if (model=='Gauss0Var1'){
    dimH0 <- 0
    dimH1 <- N_Q
  }
  if (model=='GaussEqVar'){
    dimH0 <- 2 # = 1 gaussian mean under H0 + 1 common variance for all gaussians in the model
    dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (model=='Gauss0EqVar'){
    dimH0 <- 1 # 1 common variance for all gaussians in the model
    dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (model=='GaussAffil'){
    dimH0 <- 2
    dimH1 <- 4
  }
  if (model=='Gauss02distr'){
    dimH0 <- 1
    dimH1 <- 2
  }
  if (model=='Gauss2distr'){
    dimH0 <- 2
    dimH1 <- 2
  }
  if (model=='Exp'){
    dimH0 <- 1
    dimH1 <- N_Q
  }
  dimParam <- dimPi + dimH0 + dimH1
  
  penalty <- (Q1-1)*log(n1) + (Q2-1)*log(n2) + dimParam*log(N)
  ICL <- 2*solutionThisRun$convergence$complLogLik - penalty
  return(ICL)
}



#' optimal number of SBM blocks
#'
#' returns the number of SBM blocks that maximizes the ICL
#'
#' @param bestSolutionAtQ output of \code{fitNobiSBM()}, i.e. a list of estimation results for varying number of latent blocks
#' @param exclusive whether the number of row clusters is equal to the number of column clusters. 
#' @return a list the maximal value of the ICL criterion among the provided solutions along with the best number of latent blocks
#' @export
#'
#' @examples
#' # res_gauss is the output of a call of fitNobiSBM()
#' getBestQ(res_gauss)

getBestQ <- function(bestSolutionAtQ, exclusive = TRUE){
  
  if (exclusive){
    L <- length(bestSolutionAtQ)
    Qmin <- bestSolutionAtQ[[1]]$sbmParam$Q1
    Qmax <- bestSolutionAtQ[[L]]$sbmParam$Q1
    possibleQvalues <- Qmin:Qmax
    listICL <- sapply(bestSolutionAtQ, function(elem) elem$sbmParam$ICL)
    Q1=possibleQvalues[which.max(listICL)]
    out <- list(Q1 = Q1, Q2 =Q1, ICL = max(listICL))
  } 
  
  if (!exclusive){
    dfICL <- data.frame(ICL = sapply(bestSolutionAtQ, function(el) el$sbmParam$ICL),
                        Q1 = sapply(bestSolutionAtQ, function(el) el$sbmParam$Q1),
                        Q2 = sapply(bestSolutionAtQ, function(el) el$sbmParam$Q2))
    id <- which.max(dfICL$ICL)
    out <- list(Q1 = dfICL$Q1[id], Q2 = dfICL$Q2[id], ICL = dfICL$ICL[id])
  }
  
  return(out)
}

#' plot ICL curve
#'
#' @param res output of fitNobiSBM()
#' @param exclusive whether the number of row clusters is equal to the number of column clusters. 
#'
#' @return figure of ICL curve
#' @export
#' @examples
#' # res_gauss is the output of a call of fitNobiSBM()
#' plotICL(res_gauss)
plotICL <- function(res, exclusive = TRUE){
  hatQ <- getBestQ(res, exclusive)

  if (exclusive){
    dfICL <- data.frame(ICL = sapply(res, function(el) el$sbmParam$ICL),
                        Q1 = sapply(res, function(el) el$sbmParam$Q1))
    
    g <- ggplot2::ggplot(dfICL, ggplot2::aes(x=Q1, y=ICL)) +
      ggplot2::geom_line() + 
      ggplot2::scale_color_continuous(type = "viridis") +
      ggplot2::ggtitle(paste('best number of blocks: Q = ', hatQ$Q1))    
  }
  
  if (!exclusive){
    dfICL <- data.frame(ICL = sapply(res, function(el) el$sbmParam$ICL),
                        Q1 = factor(sapply(res, function(el) el$sbmParam$Q1)),
                        Q2 = factor(sapply(res, function(el) el$sbmParam$Q2)))
    g <- ggplot2::ggplot(dfICL,
                         ggplot2::aes(x=Q1, y=Q2)) +
      ggplot2::geom_point(aes(colour=ICL)) + 
      ggplot2::scale_color_continuous(type = "viridis") +
      ggplot2::ggtitle(paste('best number of blocks: Q1 =', hatQ$Q1, ', Q2 =', hatQ$Q2))  
  }
  
  return(g)
}

