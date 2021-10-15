#' SLIP with thresholding (d version)
#'
#' @description Use SLIP with mean screening to detect abnormal data streams
#' each of which occurs at least one change.
#'
#' @param dat n x p matrix (p features, n observations)
#' @param alpha FDR nominal level
#' @param r splitting ratio, (r-1) pieces versus 1 piece
#' @param covEst Estimate covariance or not (logical); T for Est
#' @param estMthd optional estimation methods \code{c("Cholesky", "POET")}
#' @param upperPi retained proportion after thresholding \code{0 < d < 1}(default 0.5)
#' @param trueCov the true covariance matrix; only optional when \code{covEst=F}
#' @param outputW a logical parameter \code{FALSE}(default); if \code{TRUE}, the
#' W-statistics and the threshold will be returned.
#' @param outputCP logical parameter \code{FALSE}(default); if \code{TRUE}, the
#' change-point location in (0, 1) corresponding to signals will be returned.
#'
#' @return A list contains:
#'     \item{sig}{indices of signals}
#'     \item{FDP}{estiamted FDP}
#'     \item{W}{W-statistic, optional only when \code{W = TRUE}}
#'     \item{L}{threshold, optional only when \code{W = TRUE}}
#'     \item{cps}{change-points, optional only when \code{outputCP = TRUE}}
#' @export
#' @importFrom POET POETKhat POET
#' @importFrom CovTools PreEst.2017Lee
#' @examples
#'   N = 120; p = 200
#'   data = SLIP.scp.generator(N, p)
#'   SLIP.thresh.d(data$dat, 0.1)
#'
SLIP.thresh.d <- function(dat, alpha, upperPi = 0.5, r=3, covEst = T,
                          estMthd = "Cholesky", trueCov = NULL, outputW = FALSE,
                          outputCP = FALSE){

  if (is.null(dat)) {stop(list("Please input the data matrx."))}
  if(is.null(alpha)){stop(list("Please input the nominal FDR level."))}

  p = ncol(dat)
  N = nrow(dat)

  if (!covEst) {
    if (nrow(trueCov)==ncol(trueCov) & nrow(trueCov)==p){
      Sig = trueCov
    } else {
      stop(list("Something wrong with Covariance matrix! Please Check it."))
    }
  }

  SPLITE = (1:floor(N/r))*r
  m = length(SPLITE)
  n = N - m
  dat_1 = dat[-SPLITE, ]
  dat_2 = dat[SPLITE, ]

  sum.dat1 = cbind(dat_1[1, ], sapply(2:n, function(x){ colSums(dat_1[1:x, ])})) ## p x n
  cusum = abs(t(sum.dat1[, -n] - matrix(rep((1:(n-1))/n, p), byrow = T, ncol = n-1) * sum.dat1[, n])) * sqrt(n/((1:(n-1))*(n - 1:(n-1)))) ## n x p
  tau_1 = (2*r-2) + apply(cusum[(2*r-1):(n-2*r), ], 2, which.max)
  tau_2 = ceiling(tau_1/(r-1))

  xi_1 = sqrt(tau_1*(n-tau_1)/n) * sapply(1:p, function(x){
    mean(dat_1[(tau_1[x]+1):n, x]) - mean(dat_1[1:tau_1[x], x])
  })
  xi_2 = sqrt(tau_2*(m-tau_2)/m) * sapply(1:p, function(x){
    mean(dat_2[(tau_2[x]+1):m, x]) - mean(dat_2[1:(tau_2[x]-1), x])
  })

  ## Estimate covariance matrix
  if (covEst){
    dat_c = sapply(1:p, function(x){
      dat_1[, x] - c(rep(mean(dat_1[1:tau_1[x], x]), tau_1[x]), rep(mean(dat_1[(tau_1[x]+1):n, x]), (n - tau_1[x])))
    })
    if (estMthd == "Cholesky"){
      Omega = PreEst.2017Lee(dat_c, upperK = 10L)$C
      Sig = chol2inv(chol(Omega))
    } else if (estMthd == "POET"){
      K = POETKhat(t(dat_c))$K1HL
      Sig = POET(t(dat_c), K)$SigmaY
    } else {
      stop(list("Please specify an estimation method from c('POET', 'Cholesky')."))
    }
  }

  D_sig = diag(Sig)
  Gamma = GammaCalcu(Sig, tau_2, m)
  rm(Sig); gc()

  yt = as.vector(Gamma %*% xi_2)

  w1 = xi_1/sqrt(D_sig)
  sw1 = sort(abs(w1), decreasing = T)
  idx = which(abs(w1) >= sw1[floor(p*upperPi)])

  if (length(idx) == 0){
    W = rep(0, p)
  } else {
    Q = chol2inv(chol(t(Gamma[, idx])%*%Gamma[, idx]))
    bt = Q %*% (t(Gamma[, idx]) %*% yt)
    w2 = rep(0, p)
    w2[idx] = bt/sqrt(diag(Q))
    W = w1 * w2
    W[is.na(W)] = 0
  }

  s = unique(sort(abs(W), decreasing = F))
  eFDP = sapply(s, function(x){
    sum(W < -x) / max(1, sum(W > x))
  })
  if (sum(eFDP <= alpha) > 0){
    L = min(s[which(eFDP <= alpha)])
    sig = which(W > L)
  } else {
    L = Inf
    sig = vector(length = 0L)
  }

  if (length(sig)) {cps = tau_1[sig]/n} else {cps = NULL}
  res = list(sig = sig, estFDP = ifelse(length(sig), sum(W < -L)/length(sig), 0))
  if (outputW) {
    if (outputCP){
      return(c(res, list(W = W, L = L, cps = cps)))
    } else {  return(c(res, list(W = W, L = L))) }
  } else {
    if (outputCP){
      return(c(res, list(cps = cps)))
    } else { return(res) }
  }
}
