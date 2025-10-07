#' Informations on the reconciliation process
#'
#' @description
#' This function extracts reconciliation information from the output of any
#' reconciled function implemented by \pkg{FoReco}.
#'
#' @param x An output from any reconciliation function implemented by
#' \pkg{FoReco}.
#' @param verbose If \code{TRUE} (\emph{defaults}), reconciliation information
#' are printed.
#'
#' @returns A list containing the following reconciliation process information:
#'   \item{rfun}{the reconciliation function.}
#'   \item{cs_n}{the cross-sectional number of variables.}
#'   \item{te_set}{the set of temporal aggregation orders.}
#'   \item{forecast_horizon}{the forecast horizon (in temporal and
#'   cross-temporal frameworks, for the most temporally aggregated series).}
#'   \item{framework}{the reconciliation framework (cross-sectional, temporal
#'   or cross-temporal).}
#'   \item{info}{non-negative reconciled forecast convergence information.}
#'   \item{lcc}{list of level conditional reconciled forecasts (+ BU) for
#'   [cslcc], [telcc] and [ctlcc].}
#'   \item{nn}{if \code{TRUE}, all the forecasts are not negative.}
#'   \item{comb}{the covariance approximation.}
#'
#' @family Utilities
#' @export
#'
recoinfo <- function(x, verbose = TRUE) {
  if (is.null(attr(x, "FoReco"))) {
    cli_warn(c("!" = "No information available."), call = NULL)
    invisible(NULL)
  } else {
    out <- as.list(attr(x, "FoReco"))

    if (is.numeric(x)) {
      out$nn <- all(x >= 0)
      intro <- ""
    } else {
      out$nn <- out$nn
      intro <- "Probabilistic "
    }

    if (verbose) {
      if (out$rfun %in% c("cslcc", "telcc", "ctlcc")) {
        title <- "Level Conditional Coherent "
      } else if (out$rfun %in% c("csrec", "terec", "ctrec")) {
        title <- "Optimal "
      } else if (out$rfun %in% c("iterec", "tcsrec", "cstrec")) {
        title <- "Heuristic "
      } else if (out$rfun %in% c("ctbu", "csbu", "tebu")) {
        title <- "Bottom-up "
      } else if (out$rfun %in% c("cttd", "cstd", "tetd")) {
        title <- "Top-down "
      } else if (out$rfun %in% c("ctmo", "csmo", "temo")) {
        title <- "Middle-out "
      } else {
        title <- " "
      }
      cli_alert_success(
        "{.emph {intro}}{.emph {title}}{.strong {out$framework}} Forecast Reconciliation"
      )
      if (!is.null(out$rfun)) {
        cli_alert_info("{.pkg FoReco} function: {.strong {.code {out$rfun}}}")
      }

      if (!is.null(out$comb)) {
        if (length(out$comb) > 1) {
          tmp <- paste(names(out$comb), out$comb, sep = "-")
        } else {
          tmp <- out$comb
        }
        cli_alert_info("Covariance approximation: {.strong {.code {tmp}}}")
      }
      if (!is.null(out$nn)) {
        cli_alert_info("Non-negative forecasts: {.strong {.code {out$nn}}}")
      }
    }

    invisible(out)
  }
}

#' Low-level construction for reconcilied forecasts attribute foreco_info class
#'
#' `new_foreca_info()` is the contructor for the `foreca_info` class, which acompany
#' the output from the reconciliation functions in the attribute `FoReco`.
#' This is exported for extension purposes and for expert use only.
#' @aliases print.foreco_info
#' @param x A list of information related to the reconcilied forecasts.
#' @examples
#' new_foreco_info(list(
#'   framework = "Cross-sectional",
#'   forecast_horizon = 1,
#'   comb = "shr",
#'   cs_n = 3,
#'   rfun = "csrec"
#' ))
#' @export
new_foreco_info <- function(x = list()) {
  structure(x, class = "foreco_info")
}

#' @export
#' @method print foreco_info
print.foreco_info <- function(x, ...) {
  invisible(x)
}
