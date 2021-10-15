#' Bootstrap for CUSUM statistics
#'
#' @param N the length of the stream
#' @param B the number of simulations
#'
#' @return B x 1 vector; empirical distribution of the CUSUM
#' @export
#' @importFrom stats rnorm
#'
bootstrap.cusum <- function(N, B=20000){
  dat = matrix(rnorm(N*B), nrow = N)
  sum.dat = cbind(dat[1, ], sapply(2:N, function(x){ colSums(dat[1:x, ])})) ## p x n
  cusum = abs(t(sum.dat[, -N] - matrix(rep((1:(N-1))/N, B), byrow = T, ncol = N-1) * sum.dat[, N])) * sqrt(N/((1:(N-1))*(N - 1:(N-1)))) ## n x p
  ds = apply(cusum, 2, max)
  return(sort(ds, decreasing = F))
}
