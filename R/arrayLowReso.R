#' Reduce resolution of image
#'
#' @param A an array
#' @param pixel a vector determining how much pixels of initial image one pixel
#' of the output image contains; such as \code{c(2, 2, 2)} for a 3D array \code{A}
#' @param func the operation executed to get the output image; such as \code{mean}
#' or \code{max}
#'
#' @return an array with dim = \code{ceiling(dim(A)/pixel)}
#' @export
#'
#' @examples
#' A = array(1:27, c(3, 3, 3))
#' pixel = c(2, 2, 2)
#' B = arrayLowReso(A, pixel, mean)
#'
arrayLowReso <- function(A, pixel, func){

  dima = dim(A)
  if (length(dima) != length(pixel)) { stop(list("Dimensions mismatch! Please check again.")) }

  ldim = ceiling(dima/pixel)
  ndim = length(dima)

  ## generate coordinates
  coord = matrix(0, prod(dima), ndim)
  temp  = 1:prod(dima)
  for (k in ndim:1){
    if (k == 1){
      coord[, 1] = temp
    } else {
      coord[, k] = ceiling(temp/prod(dima[1:(k-1)]))
      temp = temp - (coord[, k]-1)*prod(dima[1:(k-1)])
    }
  }

  lcoord = ceiling(t(t(coord)/pixel))
  temp = rep(0, prod(dima))
  for (k in 1:ndim) {
    if (k == 1){
      temp = temp + lcoord[, k]
    } else {
      temp = temp + (lcoord[, k]-1)*prod(ldim[1:(k-1)])
    }
  }
  mask = array(temp, dima)

  return(array(tapply(A, mask, func), ldim))
}
