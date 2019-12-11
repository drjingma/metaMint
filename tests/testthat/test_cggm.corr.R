test_that("a simulation example", {

  library(MASS)
  library(glasso)
  p <- 20
  S <- diag(1,p)
  for (j in 1:(p-1)){
    S[j,j+1] <- S[j+1,j] <- 0.5
  }
  X <- mvrnorm(n=100, mu=rep(0,p), Sigma=S)
  X_cens <- apply(X,2,function(a) ifelse(a>1,a,0))
  S_hat <- cggm.corr(X_cens)
  Sinv_hat <- cggm.pcorr(X_cens,c(0.2,0.1,0.05),'glasso')
  sel_models <- cggm.stars(Sinv_hat, rep.num=5)

})

test_that("a simulation example on mclr", {

  set.seed(1)
  x <- c(rep(0,10),rpois(10,100))
  y <- c(rep(0,5),rpois(15,10))
  z <- rbind(x,y)
  z_mclr <- mclr(z)
  # compare with clr
  z_clr <- apply(z,1,compositions::clr)

})

