## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# load package: SLIP
library(SLIP)
set.seed(1234)

## -----------------------------------------------------------------------------
  N = 90
  p = 200
  data = SLIP.scp.generator(N, p, dist = "t", param = 5)

## -----------------------------------------------------------------------------
  alpha = 0.2
  
  # SLIP-thresh-C
  sig.thrsh.c = SLIP.thresh.c(data$dat, alpha)$sig
  
  # SLIP-thresh-D
  sig.thrsh.d = SLIP.thresh.d(data$dat, alpha)$sig
  
  # SLIP-lasso
  sig.lasso = SLIP.lasso(data$dat, alpha)$sig
  
  # SLIP-indep
  sig.indep = SLIP.indep(data$dat, alpha)$sig
  
  # BH-simul
  ECDF = bootstrap.cusum(N)
  sig.simul = BH.simul(data$dat, alpha, ECDF)$sig

  # BH-asymp
  sig.asymp = BH.asymp(data$dat, alpha)$sig

## -----------------------------------------------------------------------------
  # false discovery proportion (FDP) and true discovery proportion (TDP)
  sigList = list(sig.thrsh.c, sig.thrsh.d, sig.lasso, sig.indep, sig.simul, sig.asymp)
  FDP = sapply(sigList, function(sig){ length(setdiff(sig, data$index))/max(1, length(sig)) })
  TDP = sapply(sigList, function(sig){ length(intersect(sig, data$index))/length(data$index) })
  Proceudre = c("SLIP.thresh.c", "SLIP.thresh.d", "SLIP.lasso", "SLIP.indep", "BH.simul", "BH.asymp")
  res = data.frame(Proceudre, FDP = round(FDP, 4), TDP = round(TDP, 4))
  knitr::kable(t(res))

## -----------------------------------------------------------------------------
  library(SLIP)
  data = fmri.data
  (dimA = c(data$dimx, data$dimy, data$dimz))
  
  # load the ROI data used in the paper
  (dim(data$dat))
  dat = apply(data$dat[-c(1:5, 356:360), ], 2, function(X){ colMeans(matrix(X, nrow = 10)) })
  (dim(dat))

## -----------------------------------------------------------------------------
  alpha = 0.2
  sigInfo = SLIP.thresh.d(dat, alpha, estMthd = "POET", outputW = TRUE, outputCP = TRUE)
  
  # The threshold L 
  (sigInfo$L)
  
  # The discovery set
  (names(sigInfo$sig))
  
  # The estimated FDP
  (sigInfo$estFDP)
  
  # The estimated change-point location (ratio)
  (sigInfo$cps)
  
  # The W-statistics
  (sigInfo$W)

