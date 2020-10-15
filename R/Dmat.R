# @title Maps a vectorized matrix of the reconciled forecasts into a differently organized vector
#
# @description This function returns the [\code{hn(ks+m) x hn(ks+m)}]
# permutation matrix transforming vec(Y') into vec(Y^*) (see notaH.pdf)
#
# @param h forecast horizon for the lowest frequency (most temporally agregated)
# time series;
# @param kset set of p factors of \code{m}, from \code{m} to 1;
# @param n number of variables (\code{n} = \code{n_a} + \code{n_b}).
#
# @return a sparse [\code{hn(ks+m) x hn(ks+m)}] matrix
#
Dmat <- function(h, kset, n) {
  I <- .sparseDiagonal(h * sum(kset) * n)
  i <- rep(rep(rep(1:h, length(kset)), rep(rev(kset), each = h)), n)
  D <- I[order(i), ]
  return(D)
}
