#' @title Commutation matrix
#'
#' @description
#' \loadmathjax
#' This function returns the (\mjseqn{r c \times r c})
#' commutation matrix \mjseqn{\mathbf{P}} such that
#' \mjsdeqn{\mathbf{P} \mbox{vec}(\mathbf{Y}) = \mbox{vec}(\mathbf{Y}'),}
#' where \mjseqn{\mathbf{Y}} is a (\mjseqn{r \times c}) matrix.
#'
#' @param r Number of rows of \mjseqn{\mathbf{Y}}.
#' @param c Number of columns of \mjseqn{\mathbf{Y}}.
#'
#' @return A sparse (\mjseqn{r c \times r c}) matrix, \mjseqn{\mathbf{P}}.
#'
#' @references Magnus, J.R., Neudecker, H. (2019), Matrix Differential Calculus
#' with Applications in Statistics and Econometrics, third edition, New York,
#' Wiley, pp. 54-55.
#'
#' @keywords utilities
#' @family utilities
#'
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
