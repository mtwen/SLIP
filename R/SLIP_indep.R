#' SLIP without incorporating dependence
#'
#' @description Use SLIP without incorporating dependence to detect abnormal
#' data streams each of which occurs at least one change.
#'
#' @param dat n x p matrix (p features, n observations)
#' @param alpha FDR nominal level
#' @param varEst logical 0 or 1; estimate variance or not (TRUE for esitmation)
#' @param D_sig p x 1 vector, true variances of streams,
#'   optional only when \code{varEst = FALSE}
#' @param r splitting ratio, (r-1) pieces versus 1 piece
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
#' @importFrom stats var
#' @importFrom CovTools PreEst.2017Lee
#' @importFrom POET POETKhat POET
#' @examples
#'   N = 120; p = 200
#'   data = SLIP.scp.generator(N, p)
#'   SLIP.indep(data$dat, 0.1)
#'
SLIP.indep <- function(dat, alpha, r=3, varEst = T, D_sig = NULL, outputW = FALSE,
                       outputCP = FALSE){

  if (is.null(dat)) {stop(list("Please input the data matrx."))}
  if(is.null(alpha)){stop(list("Please input the nominal FDR level."))}

  p = ncol(dat)
  N = nrow(dat)

  SPLITE = (1:floor(N/r))*r
  m = length(SPLITE)
  n = N - m
  dat_1 = dat[-SPLITE, ]
  dat_2 = dat[SPLITE, ]

  sum.dat1 = cbind(dat_1[1, ], sapply(2:n, function(x){ colSums(dat_1[1:x, ])})) ## p x n
  cusum = abs(t(sum.dat1[, -n] - matrix(rep((1:(n-1))/n, p), byrow = T, ncol = n-1) * sum.dat1[, n])) * sqrt(n/((1:(n-1))*(n - 1:(n-1)))) ## n x p
  tau_1 = (2*r-2) + apply(cusum[(2*r-1):(n-2*r), ], 2, which.max)
  tau_2 = ceiling(tau_1/(r-1))

  D1 = sqrt(tau_1*(n-tau_1)/n)
  D2 = sqrt(tau_2*(m-tau_2)/m)
  gam_1 = sapply(1:p, function(x){
    mean(dat_1[(tau_1[x]+1):n, x]) - mean(dat_1[1:tau_1[x], x])
  })
  gam_2 = sapply(1:p, function(x){
    mean(dat_2[(tau_2[x]+1):m, x]) - mean(dat_2[1:(tau_2[x]-1), x])
  })
  xi_1 = D1*gam_1
  xi_2 = D2*gam_2

  ## Estimate covariance matrix
  if (varEst){
    dat_c = sapply(1:p, function(x){
      dat_1[, x] - c(rep(mean(dat_1[1:tau_1[x], x]), tau_1[x]), rep(mean(dat_1[(tau_1[x]+1):n, x]), (n - tau_1[x])))
    })
    D_sig = apply(dat_c, 2, var)
  } else if (is.null(D_sig)){
    stop(list("Please give the D_sig or set varEst=TRUE."))
  }

  W = (xi_1)*(xi_2)/D_sig

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
