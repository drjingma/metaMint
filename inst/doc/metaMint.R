## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",eval=FALSE
)

## ----setup---------------------------------------------------------------
#  library(metaMint)
#  library(dplyr)
#  library(GGally)

## ------------------------------------------------------------------------
#  data(BV, package = "metaMint")

## ----clr,fig.cap='Illustration of how clr and mclr transformation affects the marginal distribution of one OTU.'----
#  OTUs_clr <- apply(OTUs,2,compositions::clr)
#  OTUs_mclr <- t(clr.epsilon(t(OTUs)))
#  
#  ggpairs(data.frame(clr=OTUs_clr[22,],
#                     mclr=OTUs_mclr[22,]),
#          diag = list(continuous = "barDiag",binwidth = 0.5))

## ----cor,eval=FALSE------------------------------------------------------
#  est_Cor <- cggm.corr(t(OTUs_mclr))

## ----invcor,eval=FALSE---------------------------------------------------
#  est_InvCor <- cggm.pcorr(t(OTUs_mclr))

## ----MINT,eval=FALSE-----------------------------------------------------
#  dat_combined <- t(rbind(OTUs_mclr,metabs))
#  est_InvCor_combined <- cggm.pcorr(dat_combined)

## ----stars,eval=FALSE----------------------------------------------------
#  est_InvCor <- cggm.pcorr(t(OTUs_mclr))
#  stars_obj <- cggm.stars(est_InvCor,stars.thresh = 0.1, rep.num = 20)

