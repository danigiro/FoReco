#' @title Bottom-up cross-temporal forecast reconciliation
#'
#' @description
#' \loadmathjax
#' Cross temporal reconciled forecasts for all series at any temporal
#' aggregation level are computed by appropriate summation of the high-frequency
#' bottom base forecasts \mjseqn{\widehat{\mathbf{b}}_i, i = 1,...,n_b}, according to a
#' bottom-up procedure like what is currently done in both the cross-sectional
#' and temporal frameworks.
#'
#' @param Bmat (\mjseqn{n_b \times h m}) matrix of high-frequency bottom time
#' series base forecasts (\mjseqn{\widehat{\mathbf{B}}^{[1]}}).
#' \mjseqn{h} is the forecast horizon for the lowest frequency (most temporally aggregated)
#' time series.
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of the \mjseqn{p} factors
#' of \mjseqn{m}.
#'
#' @details
#' Denoting by \mjseqn{\ddot{\mathbf{Y}}} the (\mjseqn{n \times h (k^\ast + m)}) matrix containing
#' the bottom-up cross temporal reconciled forecasts, it is:
#' \mjsdeqn{\ddot{\mathbf{Y}} = \left[\begin{array}{cc}
#' \mathbf{C}\widehat{\mathbf{B}}^{[1]}\mathbf{K}_1' & \mathbf{C}\widehat{\mathbf{B}}^{[1]} \cr
#' \widehat{\mathbf{B}}^{[1]} \mathbf{K}_1' & \widehat{\mathbf{B}}^{[1]}
#' \end{array}\right],}
#' where \mjseqn{\mathbf{C}} is the cross-sectional (contemporaneous) aggregation matrix,
#' \mjseqn{\mathbf{K}_1} is the temporal aggregation matrix with \mjseqn{h=1}, and
#' \mjseqn{\widehat{\mathbf{B}}^{[1]}} is the matrix containing the high-frequency bottom
#' time series base forecasts. This expression is equivalent to
#' \mjseqn{\mbox{vec}(\ddot{\mathbf{Y}}') = \widetilde{\mathbf{S}}
#' \mbox{vec}(\widehat{\mathbf{Y}}')} for \mjseqn{h = 1}, where
#' \mjseqn{\widetilde{\mathbf{S}}} is the cross-temporal summing matrix for
#' \mjseqn{\mbox{vec}(\widehat{\mathbf{Y}}')}, and \mjseqn{\widehat{\mathbf{Y}}}
#' is the (\mjseqn{n \times h (k^\ast + m)}) matrix containing all the base forecasts
#' at any temporal aggregation order.
#'
#'
#' @return The function returns a (\mjseqn{n \times h (k^\ast + m)}) matrix of
#' bottom-up cross-temporally reconciled forecasts, \mjseqn{\ddot{\mathbf{Y}}}.
#'
#' @references
#' Di Fonzo, T., and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39(1), 39-57.
#'
#' @examples
#' data(FoReco_data)
#' # monthly base forecasts
#' hfbts <- t(FoReco2matrix(FoReco_data$base, m = 12)$k1[, -c(1:3), drop = FALSE])
#' obj <- ctbu(Bmat = hfbts, m = 12, C = FoReco_data$C)
#' rownames(obj) <- rownames(FoReco_data$base)
#'
#' @keywords bottom-up
#' @family reconciliation procedures
#' @import Matrix
#' @export
ctbu <- function(Bmat, m, C) {
  if (!(is.matrix(C) | is(C, "Matrix"))) stop("C must be a matrix", call. = FALSE)
  if (!(is.matrix(Bmat) | is(Bmat, "Matrix"))) stop("Bmat must be a matrix", call. = FALSE)

  tools <- thf_tools(m)
  kset <- tools$kset
  m <- tools$m

  if (NCOL(Bmat) %% m != 0) {
    stop("The number of Bmat's columns doesn't match with m*h", call. = FALSE)
  } else {
    h <- NCOL(Bmat) / m
  }

  # Temporal tools for the block matrix Y
  tools <- thf_tools(m = kset, h = h)
  kset <- tools$kset
  p <- tools$p
  kt <- tools$kt
  ks <- tools$ks
  K <- tools$K

  na <- NROW(C)
  nb <- NCOL(C)
  n <- na + nb

  if (NROW(Bmat) != nb) {
    stop("The number of Bmat's rows doesn't match with number of C's columns", call. = FALSE)
  }

  # Calculation of each single block of the Y matrix
  out <- Matrix(NA, nrow = n, ncol = h * (ks + m))
  out[1:na, (h * ks + 1):(h * (ks + m))] <- C %*% Bmat
  out[1:na, 1:(h * ks), drop = FALSE] <- out[1:na, (h * ks + 1):(h * (ks + m))] %*% t(K)
  out[(na + 1):n, 1:(h * ks)] <- Bmat %*% t(K)
  out[(na + 1):n, (h * ks + 1):(h * (ks + m))] <- Bmat

  out <- as.matrix(out)

  rownames(out) <- paste("serie", 1:n, sep = "")
  colnames(out) <- paste("k", rep(kset, h * (m/kset)), "h",
    do.call("c", as.list(sapply(
      (m/kset) * h,
      function(x) seq(1:x)
    ))),
    sep = ""
  )
  return(out)
}
