#' @title Bottom-up Cross-temporal forecast reconciliation
#'
#' @description Cross temporal reconciled forecasts for all series at any temporal
#' aggregation level can be easily computed by appropriate summation of the high-frequency
#' bottom base forecasts \eqn{\hat{\textbf{b}}_i, i = 1,...,n_b, \;}{} according to a
#' bottom-up procedure like what is currently done in both the cross-sectional
#' and temporal frameworks.
#'
#' @param Bmat (\code{nb x (h m)}) matrix of high-frequency bottom time series base forecasts.
#' \code{h} is the forecast horizon for the lowest frequency (most temporally aggregated)
#' time series.
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param m Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation).
#'
#' @return The function returns a (\code{n x (h (k* + m))}) matrix of
#' bottom-up cross-temporally reconciled forecasts.
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' @examples
#' data(FoReco_data)
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base),
#'                                     split = "_"))[1, ] == "k1")
#' hfbts <- FoReco_data$base[-c(1:3), id]
#' obj <- ctbu(Bmat = hfbts, m = 12, C = FoReco_data$C)
#' rownames(obj) <- rownames(FoReco_data$base)
#'
#' @keywords reconciliation bottom-up
#' @import Matrix
#' @export
ctbu <- function(Bmat, m, C) {
  if (!(is.matrix(C))) stop("C must be a matrix", call. = FALSE)
  if (!(is.matrix(Bmat))) stop("Bmat must be a matrix", call. = FALSE)

  if (NCOL(Bmat) %% m != 0) {
    stop("The number of Bmat's columns doesn't match with m*h", call. = FALSE)
  } else {
    h <- NCOL(Bmat) / m
  }

  # Temporal tools for the block matrix Y
  tools <- thf_tools(m, h = h)
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
  colnames(out) <- paste("k", rep(kset, h * rev(kset)), "h",
    do.call("c", as.list(sapply(
      rev(kset) * h,
      function(x) seq(1:x)
    ))),
    sep = ""
  )
  return(out)
}
