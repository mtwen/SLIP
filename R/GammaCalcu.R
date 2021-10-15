GammaCalcu <- function(Sig, tau_2, m){
  p = length(tau_2)

  temp.row = tau_2[as.vector(row(diag(p)))]
  temp.col = tau_2[as.vector(col(diag(p)))]
  INDI = which(temp.row > temp.col)
  ma = rep(0, p^2)
  ma[INDI] = temp.row[INDI]
  ma[-INDI] = temp.col[-INDI]
  mi = rep(0, p^2)
  mi[INDI] = temp.col[INDI]
  mi[-INDI] = temp.row[-INDI]

  Sigma = matrix(sqrt(mi*(m - ma)/(ma*(m - mi)))*as.vector(Sig), nrow = p)
  Gamma = chol(chol2inv(chol(Sigma)))

  return(Gamma)
}
