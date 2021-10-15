corrStructFactor <- function(p){
  b1 = rnorm(p)
  b2 = rnorm(p)
  b3 = rnorm(p)
  B = cbind(b1, b2, b3)
  Sigma = B %*% t(B) + 0.5 * diag(p)
  Sig0 = diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
  return(Sig0)
}
