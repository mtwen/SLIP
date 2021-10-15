corrStructSparse <- function(p){
  U = (2*rbinom(p, 1, 0.5)-1) * runif(p, 1, 2)
  Gamma = matrix(0, p, p)
  for (i in 1:p) Gamma[i, sample(1:p, 1)]=U[i]
  S = Gamma%*%t(Gamma)+diag(p)
  Sig0 = diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
  return(Sig0)
}
