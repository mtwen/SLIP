#' Generate n-by-p matrix with one changes in each abnormal datastream
#'
#' @param n the number of observations
#' @param p the number of datastreams
#' @param cov covariance type; \code{c("CS", "Factor", "AR1", "Block", "Sparse",
#'     "Identity", "Given")}
#' @param rho optional when using \code{cov} from \code{c("CS", "AR1", "Block")},
#'     rho in [0, 1]
#' @param dist noise distribution \code{c("Gaussian", "t", "exp", "chisq")}
#' @param ratio abnormal streams/total streams = ratio, ratio in [0, 1]
#' @param delta the magnitudes of changes (delta >= 0)
#'     lie in [delta-0.1, delta+0.1] with equally probable sign from {+, -}
#' @param Sigma the covariance matrix, optional only when \code{cov=("Given")}
#' @param varrho the parameter avoiding the boundary problem, [0, 0.5)
#' @param param the parameter when using \code{dist} from
#'     \code{c("t", "chisq", "exp")}
#'
#' @return data with changes in some datastreasm
#'     \item{dat}{n-by-p data matrix}
#'     \item{index}{those datastreams containing changes}
#'     \item{mu}{changes}
#'     \item{loc}{locations of changes}
#' @export
#' @importFrom stats rnorm rt rchisq rexp rbinom runif
#'
SLIP.scp.generator <- function(n, p, cov = "Identity", dist="Gaussian",
                               rho = NULL, ratio=0.15, delta=1, varrho = 0.05,
                               Sigma = NULL, param = NULL){

  if (varrho < 0 | varrho >= 0.5) { stop(list("Please specify varrho in [0, 0.5)."))}
  if (ratio < 0 | ratio > 1) { stop(list("Please specify ratio in [0, 1].")) }
  if (delta < 0) {stop(list("Please specify delta >= 0."))}

  if (cov == "CS") {
    if (is.numeric(rho)){
      if (rho<0 | rho>1) {stop(list("Please specify rho in [0, 1]"))}
      Sig0 = (1-rho)*diag(p) + rho * rep(1, p) %*% t(rep(1, p))
    } else { stop(list("Please input rho in [0, 1] when using CS structure.")) }
  } else if (cov == "Factor") { Sig0 = corrStructFactor(p)
  } else if (cov == "AR1") {
    if (is.numeric(rho)){
      if (rho<0 | rho>1) {stop(list("Please specify rho in [0, 1]"))}
      Sig0 = rho^(abs(row(diag(p)) - col(diag(p))))
    } else { stop(list("Please input rho in [0, 1] when using AR(1) structure.")) }
  } else if (cov == "Block") {
    if (is.numeric(rho)){
      if (rho<0 | rho>1) {stop(list("Please specify rho in [0, 1]"))}
      Sig0 = corrStructBlock(p, rho)
    } else { stop(list("Please input rho in [0, 1] when using Block structure.")) }
  } else if (cov == "Sparse") { Sig0 = corrStructSparse(p)
  } else if (cov == "Identity") { Sig0 = diag(p)
  } else if (cov == "Given"){
    if (is.null(Sigma)){ stop(list("Please specify the Sigma.")) }
    if (nrow(Sigma) != ncol(Sigma) | nrow(Sigma) != p) {
      stop(list("Check the dimension of Sigma."))
    }
    Sig0 = Sigma
  } else { stop(list("Please choose one optional covariance type.")) }

  LMtx = chol(Sig0)
  if (dist=='Gaussian'){ dat = matrix(rnorm(n*p), nrow = n) %*% LMtx }
  else if (dist=='t'){
    if (is.null(param)){
      stop(list("Please give the df of t-distribution."))
    } else {
      k = param; dat = matrix((rt(n*p, df = k))/sqrt(k/(k-2)), nrow = n) %*% LMtx
    }
  }
  else if(dist=='chisq'){
    if (is.null(param)){
      stop(list("Please give the df of chisquare distribution."))
    } else {
      k = param; dat = matrix((rchisq(n*p, df = k) - k)/sqrt(2*k), nrow = n) %*% LMtx
    }
  }
  else if(dist=='exp'){
    if (is.null(param)){
      stop(list("Please give the rate of exponential distribution."))
    } else {
      rate = param; dat = matrix((rexp(n*p, rate = rate) - 1/rate)/(1/rate), nrow = n) %*% LMtx
    }
  }
  else { stop(list("Please choose one optional distribution of noise.")) }

  p1 = floor(p*ratio)
  index = sort(sample(1:p, p1, replace = F))
  mu = runif(p1, delta-0.1, delta+0.1) * (2*rbinom(p1, 1, 0.5) - 1)
  loc = sample((floor(n*varrho)+1):(n-1 - floor(n*varrho)), p1, replace = T)
  tdat = sapply(1:p1, function(x){
    c(rep(0, loc[x]), rep(mu[x], (n-loc[x])))
  })
  dat[, index] = dat[, index] + tdat

  return(list(dat = dat, index = index, mu = mu, loc = loc))

}
