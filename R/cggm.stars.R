
# Stability approach to regularization selection
#' @param est An object from fitting one of the CGGM methods.
#' @param stars.thresh The variability threshold in stars. The default value is \code{0.1}. An alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}.
#' @param stars.subsample.ratio The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size. Only applicable when \code{criterion = "stars"}.

cggm.stars <- function(est, stars.thresh=0.1,stars.subsample.ratio=NULL, rep.num=20, verbose=TRUE){
  n = nrow(est$data)
  d = ncol(est$data)
  nlambda = length(est$lambda)

  if(is.null(stars.subsample.ratio))
  {
    if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
    if(n<=144) stars.subsample.ratio = 0.8
  }

  est$merge = list()
  for(i in 1:nlambda) est$merge[[i]] = Matrix(0,d,d)

  for(i in 1:rep.num)
  {
    if(verbose)
    {
      mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/rep.num), "%"), collapse="")
      cat(mes, "\r")
      flush.console()
    }
    ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)

    if(est$method == "direct")
      tmp = cggm.direct(est$data[ind.sample,], lambda = est$lambda, method='glasso')$path
    if(est$method == "spring")
      tmp = cggm.spring(est$data[ind.sample,], lambda = est$lambda, method='glasso')$path
    if(est$method == "glasso")
      tmp = cggm.glasso(est$data[ind.sample,], lambda = est$lambda)$path
    if(est$method == "gcoda")
      tmp = cggm.gcoda(est$data[ind.sample,], lambda = est$lambda, counts = T)$path

    for(i in 1:nlambda)
      est$merge[[i]] = est$merge[[i]] + tmp[[i]]

    rm(ind.sample,tmp)
    gc()
  }

  if(verbose){
    mes = "Conducting Subsampling....done.                 "
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }

  est$variability = rep(0,nlambda)
  for(i in 1:nlambda){
    est$merge[[i]] = est$merge[[i]]/rep.num
    est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
  }

  est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
  est$refit = est$path[[est$opt.index]]
  est$opt.lambda = est$lambda[est$opt.index]
  est$opt.sparsity = est$sparsity[est$opt.index]
  est$opt.icov = est$icov[[est$opt.index]]

  class(est) = "select"
  return(est)
}
