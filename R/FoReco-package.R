#' FoReco: point forecast reconciliation
#'
#' An R package offering classical (bottom-up and top-down), and modern (optimal and heuristic combination)
#' forecast reconciliation procedures for cross-sectional, temporal, and cross-temporal
#' linearly constrained time series.
#'
#' @author Tommaso Di Fonzo and Daniele Girolimetto, Department of Statistical Sciences, University of Padua (Italy).
#'
#' @details
#' The \code{FoReco} package is designed for point forecast reconciliation, a
#' post-forecasting process aimed to improve the accuracy of the base
#' forecasts for a system of linearly constrained (e.g. hierarchical/grouped) time series.
#' The main functions are:
#'
#' \describe{
#'   \item{\code{\link[FoReco]{htsrec}():}}{cross-sectional (contemporaneous) forecast reconciliation.}
#'   \item{\code{\link[FoReco]{thfrec}():}}{forecast reconciliation for a single time series through temporal hierarchies.}
#'   \item{\code{\link[FoReco]{lccrec}():}}{level conditional forecast reconciliation for genuine hierarchical/grouped time series.}
#'   \item{\code{\link[FoReco]{tdrec}():}}{top-down (cross-sectional, temporal, cross-temporal) forecast reconciliation for genuine hierarchical/grouped time series.}
#'   \item{\code{\link[FoReco]{ctbu}():}}{bottom-up cross-temporal forecast reconciliation.}
#'   \item{\code{\link[FoReco]{tcsrec}():}}{heuristic first-temporal-then-cross-sectional cross-temporal forecast reconciliation.}
#'   \item{\code{\link[FoReco]{cstrec}():}}{heuristic first-cross-sectional-then-temporal cross-temporal forecast reconciliation.}
#'   \item{\code{\link[FoReco]{iterec}():}}{heuristic iterative cross-temporal forecast reconciliation.}
#'   \item{\code{\link[FoReco]{octrec}():}}{optimal combination cross-temporal forecast reconciliation.}
#' }
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' Di Fonzo, T., Girolimetto, D. (2021), \emph{Forecast combination based forecast reconciliation: insights and extensions} (mimeo).
#'
#' @import mathjaxr
#' @docType package
#' @name FoReco-package
#' @keywords package
NULL

.onAttach <- function(...) {
  packageStartupMessage(cli::rule(
    right = paste0("FoReco ", utils::packageVersion("FoReco"))
  ))
}


