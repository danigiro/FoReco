#' @title Heuristic first-cross-sectional-then-temporal cross-temporal forecast reconciliation
#'
#' @description
#' \loadmathjax
#' Cross-temporal forecast reconciliation according to the heuristic procedure by
#' Kourentzes and Athanasopoulos (2019), where the order of application of the
#' two reconciliation steps (temporal-first-then-cross-sectional, as in the function
#' \code{\link[FoReco]{tcsrec}()}), is inverted. The function
#' \code{\link[FoReco]{cstrec}()} performs cross-sectional reconciliation
#' (\code{\link[FoReco]{htsrec}()}) first, then temporal reconciliation
#' (\code{\link[FoReco]{thfrec}()}), and finally applies the average of the
#' projection matrices obtained in the second step to the one dimensional
#' reconciled values obtained in the first step.
#'
#' @param basef  (\mjseqn{n \times h(k^\ast+m)}) matrix of base forecasts to be
#' reconciled, \mjseqn{\widehat{\mathbf{Y}}}; \mjseqn{n} is the total number of variables,
#' \mjseqn{m} is the highest time frequency, \mjseqn{k^\ast} is the sum of (a
#' subset of) (\mjseqn{p-1}) factors of \mjseqn{m}, excluding \mjseqn{m}, and
#' \mjseqn{h} is the forecast horizon for the lowest frequency time series.
#' Each row identifies a time series, and the forecasts are ordered as
#' [lowest_freq' ...  highest_freq']'.
#' @param hts_comb,thf_comb Type of covariance matrix (respectively
#' (\mjseqn{n \times n}) and (\mjseqn{(k^\ast + m) \times (k^\ast + m)})) to
#' be used in the cross-sectional and temporal reconciliation, see more in
#' \code{comb} param of \code{\link[FoReco]{htsrec}()} and
#' \code{\link[FoReco]{thfrec}()}.
#' @param res (\mjseqn{n \times N(k^\ast + m)}) matrix containing the residuals
#' at all the temporal frequencies ordered as [lowest_freq' ...  highest_freq']'
#' (columns) for each variable (row), needed to estimate the covariance matrix
#' when \code{hts_comb =} \code{\{"wls",} \code{"shr",} \code{"sam"\}} and/or
#' \code{hts_comb =} \code{\{"wlsv",} \code{"wlsh",} \code{"acov",}
#' \code{"strar1",} \code{"sar1",} \code{"har1",} \code{"shr",} \code{"sam"\}}.
#' The rows must be in the same order as \code{basef}.
#' @param ... any other options useful for \code{\link[FoReco]{htsrec}()} and
#' \code{\link[FoReco]{thfrec}()}, e.g. \code{m}, \code{C} (or \code{Ut} and
#' \code{nb}), \code{nn} (for non-negative reconciliation only at the first step),
#' \code{mse}, \code{corpcor}, \code{type}, \code{sol}, \code{settings},
#' \code{W}, \code{Omega},...
#'
#' @details
#' \strong{Warning},
#' the two-step heuristic reconciliation allows non negativity constraints only in the first step.
#' This means that it is not guaranteed the non-negativity of the final reconciled values.
#'
#' @return
#' The function returns a list with two elements:
#' \item{\code{recf}}{(\mjseqn{n \times h(k^\ast + m)}) reconciled forecasts matrix, \mjseqn{\widetilde{\mathbf{Y}}}.}
#' \item{\code{M}}{Matrix which transforms the uni-dimensional reconciled forecasts of step 1 (projection approach) .}
#'
#' @references
#' Di Fonzo, T., and Girolimetto, D. (2021), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, in press.
#'
#' Kourentzes, N., Athanasopoulos, G. (2019), Cross-temporal coherent forecasts
#' for Australian tourism, \emph{Annals of Tourism Research}, 75, 393-409.
#'
#' \enc{Schäfer}{Schafer}, J.L., Opgen-Rhein, R., Zuber, V., Ahdesmaki, M.,
#' Duarte Silva, A.P., Strimmer, K. (2017), \emph{Package `corpcor'}, R
#' package version 1.6.9 (April 1, 2017), \href{https://CRAN.R-project.org/package=corpcor}{https://CRAN.R-project.org/package= corpcor}.
#'
#' \enc{Schäfer}{Schafer}, J.L., Strimmer, K. (2005), A Shrinkage Approach to Large-Scale Covariance
#' Matrix Estimation and Implications for Functional Genomics, \emph{Statistical
#' Applications in Genetics and Molecular Biology}, 4, 1.
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd, S. (2020). OSQP:
#' An Operator Splitting Solver for Quadratic Programs, \emph{Mathematical Programming Computation},
#' 12, 4, 637-672.
#'
#' Stellato, B., Banjac, G., Goulart, P., Boyd, S., Anderson, E. (2019), OSQP:
#' Quadratic Programming Solver using the 'OSQP' Library, R package version 0.6.0.3
#' (October 10, 2019), \href{https://CRAN.R-project.org/package=osqp}{https://CRAN.R-project.org/package=osqp}.
#'
#' @keywords heuristic
#' @family reconciliation procedures
#'
#' @examples
#' data(FoReco_data)
#' obj <- cstrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'               hts_comb = "shr", thf_comb = "acov", res = FoReco_data$res)
#'
#' @usage cstrec(basef, hts_comb, thf_comb, res, ...)
#'
#' @export
cstrec <- function(basef, hts_comb, thf_comb, res, ...) {

  arg_input <- list(...)

  if (missing(basef)) {
    stop("the argument basef is not specified", call. = FALSE)
  }
  if (all(names(arg_input) != "m")) {
    stop("the argument m is not specified", call. = FALSE)
  } else {
    m <- arg_input$m
  }
  if (missing(thf_comb)) {
    stop("the argument thf_comb is not specified", call. = FALSE)
  }
  if (missing(hts_comb)) {
    stop("the argument hts_comb is not specified", call. = FALSE)
  }

  tools <- thf_tools(m)
  kset <- tools$kset
  m <-tools$m
  kt <- tools$kt

  # Step 1
  arg_hts <- names(as.list(args(htsrec)))
  arg_hts <- arg_hts[!(arg_hts %in% c("basef", "keep", "res", "", "comb", "bounds"))]

  # Create list with lenght p, with time by time temporally reconciled forecasts matrices
  h <- NCOL(basef) / kt
  Y <- lapply(kset, function(x) basef[, rep(kset, (m/kset) * h) == x, drop = FALSE])

  if (missing(res)) {
    Y1 <- lapply(Y, function(x){
      obj <- do.call("htsrec", c(list(basef = t(x), comb = hts_comb),
                                 arg_input[which(names(arg_input) %in% arg_hts)]))
      return(t(obj[["recf"]]))
    })
  } else {
    # Create list with lenght p, with time by time temporally reconciled residuals matrices
    r <- NCOL(res) / kt
    E <- lapply(kset, function(x) res[, rep(kset, (m/kset) * r) == x, drop = FALSE])

    ## list of time by time cross sectional M matrix
    Y1 <- mapply(function(Y, E){
      obj <- do.call("htsrec", c(list(basef = t(Y), comb = hts_comb, res = t(E)),
                                 arg_input[which(names(arg_input) %in% arg_hts)]))
      return(t(obj[["recf"]]))
    }, Y = Y, E = E)
  }

  Y1 <- do.call("cbind", Y1)

  # Step 2
  arg_thf <- names(as.list(args(thfrec)))
  arg_thf <- arg_thf[!(arg_thf %in% c("basef", "keep", "res", "", "comb", "m", "nn", "bounds", "nn_type"))]

  if (missing(res)) {
    M <- lapply(split(Y1, row(Y1)), function(x) {
      obj <- do.call("thfrec", c(list(basef = x, m = kset, comb = thf_comb),
                                 arg_input[which(names(arg_input) %in% arg_thf)]))
      return(obj[["M"]])
    })
  } else {
    M <- mapply(function(Y, X) {
      obj <- do.call("thfrec", c(list(basef = Y, m = kset, comb = thf_comb, res = X),
                                 arg_input[which(names(arg_input) %in% arg_thf)]))
      return(obj[["M"]])
    },
    Y = split(Y1, row(Y1)), X = split(res, row(res)), SIMPLIFY = FALSE)
  }

  ## Step 3: Cross-Temporal reconciled forecasts with heuristic
  #meanM <- apply(simplify2array(lapply(M, as.matrix)), c(1, 2), sum) / length(M)
  meanM <- Reduce("+", M)/length(M)

  if (h == 1) {
    Y3 <- Y1 %*% t(meanM)
  } else {
    h_pos <- rep(rep(1:h, length(kset)), rep((m/kset), each = h))
    D <- Dmat(h = h, m = kset, n = 1)

    Y3 <- list()
    for (i in 1:h) {
      Y3[[i]] <- Y1[, which(h_pos == i)] %*% t(meanM)
    }
    Y3 <- do.call("cbind", Y3)
    Y3 <- Y3 %*% D
  }

  rec_sol <- list()
  rec_sol$recf <- Y3
  rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:NROW(rec_sol$recf), sep = "") else rownames(basef)
  colnames(rec_sol$recf) <- paste("k", rep(kset, h * (m/kset)), "h",
                                  do.call("c", as.list(sapply(
                                    (m/kset) * h,
                                    function(x) seq(1:x)
                                  ))),
                                  sep = ""
  )
  rec_sol$M <- meanM
  rec_sol$nn_check <- sum(rec_sol$recf < 0)
  return(rec_sol)
}
