#' @param mu A 2 by 1 vector of mean
#' @param sd A 2 by 1 vector of standard deviation
#' @param rho A scalar of correlation
cgm.bivariate.loglik.single <- function(y, y_L, mu, sd, rho){
  S <- diag(sd^2)
  S[1,2] <- S[2,1] <- rho*prod(sd)
  if (sum(y-y_L>1e-08)==2){
    dmvnorm(y,mean=mu,sigma=S,log = TRUE)
  } else if (sum(y-y_L>1e-08)==0){
    log(pmvnorm(upper=y_L,mean=mu,sigma=S))
  } else if (y[2] > y_L[2]){
    mu_cond <- mu[1] + (y[2] - mu[2]) * rho * sd[1] / sd[2]
    sig_cond <- sqrt(1-rho^2) * sd[1]
    log(pnorm(y_L[1],mu_cond,sig_cond) * dnorm(y[2], mean = mu[2], sd = sd[2]))
  } else if (y[1] > y_L[1]){
    mu_cond <- mu[2] + (y[1] - mu[1]) * rho * sd[2] / sd[1]
    sig_cond <- sqrt(1-rho^2) * sd[2]
    log(pnorm(y_L[2],mu_cond,sig_cond) * dnorm(y[1], mean = mu[1], sd = sd[1]))
  }
}
cgm.mixed.loglik.single <- function(y, y_L, mu, sd, rho){
  S <- diag(sd^2)
  S[1,2] <- S[2,1] <- rho*prod(sd)
  if (y[1]>y_L[1]){
    dmvnorm(y,mean=mu,sigma=S,log = TRUE)
  } else {
    mu_cond <- mu[1] + (y[2] - mu[2]) * rho * sd[1] / sd[2]
    sig_cond <- sqrt(1-rho^2) * sd[1]
    log(pnorm(y_L[1],mu_cond,sig_cond) * dnorm(y[2], mean = mu[2], sd = sd[2]))
  }
}

cgm.marginal <- function(y){
  # set to be the smallest value if the smallest value appears more than 10% of the total.
  n <- nrow(y)
  p <- ncol(y)
  y_L <- rep(-Inf,p)
  for (j in 1:p){
    y_L[j] <- ifelse(length(which(y[,j]==min(y[,j])))>0, min(y[,j]), -Inf)
  }

  # estimate marginal mean and sd
  mu <- rep(0,p)
  sig <- rep(0,p)
  for (j in 1:p){
    if (is.infinite(y_L[j])){
      mu[j] <- mean(y[,j]); sig[j] <- sd(y[,j])
    } else {
      ddd <- tibble::tibble(y=y[,j])
      est.obj <- censReg::censReg( y ~ 1, data = ddd, left=y_L[j])
      mu[j] <- as.numeric(est.obj$estimate[1])
      sig[j] <- exp(as.numeric(est.obj$estimate[2]))
    }
  }
  return(list(mu=mu, sd=sig, y_L = y_L))
}

# Function to estimate the correlation for bivariate cgm model
#' @param y An n by 2 data matrix
#' @param y_L A 2 by 1 vector indicating the left limit for each of the censored dependent variable. Default to be NULL.
#' @param rho The correlation between the two variables
cgm.bivariate <- function(y, y_L, mu, sig){
  rho.init <- cor(y[,1],y[,2])

  if (rho.init>0){
    rho.interval <- c(rho.init - 0.5,rho.init + (1-rho.init)/2)
  } else{
    rho.interval <- c(rho.init - (rho.init+1)/2, rho.init+0.5)
  }
  if (sum(is.infinite(y_L))==0){
    myf <- function(b) sum(apply(y,1,function(a) cgm.bivariate.loglik.single(a,y_L,mu,sig,b)))
  } else if (is.infinite(y_L[2])){
    myf <- function(b) sum(apply(y,1,function(a) cgm.mixed.loglik.single(a,y_L,mu,sig,b)))
  } else if (is.infinite(y_L[1])){
    # reverse the order
    y <- y[,2:1]
    y_L <- y_L[2:1]
    mu <- mu[2:1]
    sig <- sig[2:1]
    myf <- function(b) sum(apply(y,1,function(a) cgm.mixed.loglik.single(a,y_L,mu,sig,b)))
  }

  op <- tryCatch(optimize(myf, rho.interval, maximum = TRUE)[1], error=function(e) 100)

  ifelse(op==100, 0, unlist(op))
}

# Covariance estimation based on the cgm model
# type1 = "continuous", type2 = "censored"
cgm.covariance.mixed <- function(X1=NULL, X2, use.nearPD=TRUE){
  if (is.null(X1)) {
    p <- 0; Q <- ncol(X2);
    y = X2;
  } else {
    p <- ncol(X1); Q <- ncol(X1) + ncol(X2)
    y <- cbind(X1,X2);
  }
  obj <- cgm.marginal(y)
  corr_mat <- diag(1,Q)
  if (!is.null(X1)){
    corr_mat[1:p,1:p] <- cor(X1)
  }
  for (j in (p+1):Q){
    # cat('index...',j,'...\n')
    if (j==1){next;}
    for (k in 1:(j-1)){
      corr_mat[j,k] <- corr_mat[k,j] <- cgm.bivariate(y[,c(j,k)],y_L = obj$y_L[c(j,k)],mu=obj$mu[c(j,k)],sig = obj$sd[c(j,k)])
    }
  }

  if ( use.nearPD == TRUE & min(eigen(corr_mat)$values) < 0 ) {
    message(" minimum eigenvalue of correlation estimator is ", min(eigen(corr_mat)$values), "\n nearPD is used")
    corr_mat <- as.matrix(Matrix::nearPD(corr_mat, corr = TRUE, posd.tol=5*1e-04)$mat)
  }

  return(list(mu=obj$mu, sigma=obj$sd, corr=corr_mat))
}


