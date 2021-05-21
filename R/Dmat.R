# @title Maps a vectorized matrix of the reconciled forecasts into a differently organized vector
#
# @description
# \loadmathjax
# This function returns the [\mjseqn{hn(k^\ast+m) \times hn(k^\ast+m)}]
# permutation matrix transforming \mjseqn{\mbox{vec}(\mathbf{Y}')} into
# \mjseqn{\mbox{vec}(\mathbf{Y}_h')}:
# \mjsdeqn{\mathbf{D}_h \mbox{vec}(\mathbf{Y}') = \mbox{vec}(\mathbf{Y}_h')}
# where \mjseqn{\mathbf{Y}_h'} is a
# \mjseqn{h \times n(k^\ast+m)} with column ordered as [lowest_freq' ...  highest_freq']
# repeat for \mjseqn{n} series.
#
# @param h forecast horizon for the lowest frequency (most temporally agregated)
# time series;
# @param m Highest available sampling frequency per seasonal cycle (max. order
# of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
# of \mjseqn{m}.
# @param n number of the cross-sectional variables (\code{n} = \code{n_a} + \code{n_b}).
#
# @return A sparse [\code{hn(ks+m) x hn(ks+m)}] matrix.
#
# @examples
# data(FoReco_data)
# # An example in the temporal frameworks
# # top ts residuals ([h*lowest_freq' ...  h*highest_freq']')
# topres <- FoReco_data$res[1, ]
# D1 <- FoReco:::Dmat(m = c(1,2,3,4,6,12), n = 1, h = 14)
# Yh_t <- matrix(D1%*%topres, nrow = 14, byrow = TRUE)
# # check
# all(as.vector(D1%*%topres) == as.vector(t(Yh_t)))
#
# # An example in the cross-temporal frameworks
# Dn <- FoReco:::Dmat(m = c(1,2,3,4,6,12), n = 8, h = 14)
# Yh_ct <- matrix(Dn%*%as.vector(t(FoReco_data$res)), nrow = 14, byrow = TRUE)
# # check
# all(Dn%*%as.vector(t(FoReco_data$res)) == as.vector(t(Yh_ct)))
#
Dmat <- function(h, m, n){
  if(length(m)==1){
    kset <- 1
  }else{
    kset <- sort(m, decreasing = TRUE)
    m <- max(m)
  }

  I <- .sparseDiagonal(h * sum(m/kset) * n)
  i <- rep(rep(rep(1:h, length(kset)), rep(m/kset, each = h)), n)
  D <- I[order(i), ]
  return(D)
}
