#' BH with p-values calculated by Bootstrap
#'
#' @param dat n x p data matrix (p features, n observations)
#' @param alpha FDR nominal level
#' @param ECDF empirical null distribution of CUSUM
#' @param varEst logical 0 or 1; estimate variance or not (TRUE for esitmation)
#' @param D_sig p x 1 vector, true variances of streams,
#'   optional only when \code{varEst = FALSE}
#' @param outputCP logical parameter \code{FALSE}(default); if \code{TRUE}, the
#' change-point location in (0, 1) corresponding to signals will be returned.
#'
#' @return A list contains:
#'     \item{sig}{indices of signals}
#'     \item{FDP}{estiamted FDP}
#'     \item{cps}{change-points, optional only when \code{outputCP = TRUE}}
#' @export
#'
#' @examples
#'   N = 120; p = 200; B = 1000
#'   data = SLIP.scp.generator(N, p)
#'   ECDF = bootstrap.cusum(N, B)
#'   BH.simul(data$dat, 0.1, ECDF)
#'
BH.simul <- function(dat, alpha, ECDF, varEst = T, D_sig=NULL, outputCP = FALSE) {

  if (is.null(dat)) {stop(list("Please input the data matrx."))}
  if(is.null(alpha)){stop(list("Please input the nominal FDR level."))}
  if(is.null(ECDF)){
    ECDF = bootstrap.cusum(nrow(dat))
  }

  p = ncol(dat)
  N = nrow(dat)
  B = length(ECDF)

  sum.dat = cbind(dat[1, ], sapply(2:N, function(x){ colSums(dat[1:x, ])})) ## p x n
  cusum = abs(t(sum.dat[, -N] - matrix(rep((1:(N-1))/N, p), byrow = T, ncol = N-1) * sum.dat[, N])) * sqrt(N/((1:(N-1))*(N - 1:(N-1)))) ## n x p
  tau = apply(cusum, 2, which.max)

  if (varEst){
    dat_c = sapply(1:p, function(x){
      dat[, x] - c(rep(mean(dat[1:tau[x], x]), tau[x]), rep(mean(dat[(tau[x]+1):N, x]), (N-tau[x])))
    })
    D_sig = apply(dat_c, 2, var)
  } else if (is.null(D_sig)){
    stop(list("Please give the D_sig or set varEst=TRUE."))
  }

  ds = sapply(1:p, function(x){
    sqrt(tau[x]*(N-tau[x])/N)*abs(mean(dat[1:tau[x], x]) - mean(dat[(tau[x]+1):N, x]))
  })
  xi = ds/sqrt(D_sig)
  pvalues = sapply(xi, function(x){
    sum(ECDF >= x)/B
  })

  temp = (1:p)/p
  spv = sort(pvalues, decreasing = F)
  if (sum(spv <= temp*alpha)) {
    thrsh = spv[max(which(spv <= temp*alpha))]
    sig = which(pvalues <= thrsh)
  } else {
    sig = vector(length = 0L)
  }

  if (length(sig)){ cps = tau[sig]/N} else{ cps = NULL }
  res = list(sig = sig, estFDP = ifelse(length(sig), (p*thrsh)/length(sig), 0))
  if (outputCP){
    return(c(res, list(cps = cps)))
  } else { return(res)}
}
