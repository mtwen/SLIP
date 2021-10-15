corrStructBlock <- function(p, rho){
  UU = matrix(c(1, rho, rho, 1), nrow = 2)
  S = simex::diag.block(UU, floor(p/2))
  return(S)
}
