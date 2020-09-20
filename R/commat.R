#' @title Commutation matrix
#'
#' @description This function returns the (\code{(r c) x (r c)}) commutation matrix \code{P} such that
#' \deqn{P vec(Y) = vec(Y'),}
#' where \code{Y} is a (\code{r x c}) matrix.
#'
#' @param r Number of rows of \code{Y}.
#' @param c Number of columns of \code{Y}.
#'
#' @return A sparse (\code{(r c) x (r c)}) matrix.
#'
#' @references Magnus, J.R., Neudecker, H. (2019), Matrix Differential Calculus with Applications
#' in Statistics and Econometrics, third edition, New York, Wiley, pp. 54-55.
#'
#' @keywords utilities
#' @examples
#' Y <- matrix(rnorm(30), 5, 6)
#' P <- commat(5, 6)
#' P %*% as.vector(Y) == as.vector(t(Y)) # check
#' @import Matrix
#' @export
commat <- function(r, c) {
  I <- Matrix(1:(r * c), c, r, byrow = T) # initialize a matrix of indices of size (c x r)
  I <- as.vector(I) # vectorize the required indices
  P <- Diagonal(r * c) # Initialize an identity matrix
  P <- P[I, ] # Re-arrange the rows of the identity matrix
  return(P)
}
