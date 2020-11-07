#' @title Heuristic first-cross-sectional-then-temporal cross-temporal forecast reconciliation
#'
#' @description
#' The order of application of the two reconciliation steps proposed by Kourentzes and Athanasopoulos (2019),
#' implemented in the function \code{\link[FoReco]{tcsrec}}, may be inverted. The function
#' \code{\link[FoReco]{cstrec}} performs cross-sectional reconciliation (\code{\link[FoReco]{htsrec}})
#' first, then temporal reconciliation (\code{\link[FoReco]{thfrec}}), and finally applies the
#' average of the projection matrices obtained in the second step to the one dimensional reconciled
#' values obtained in the first step.
#'
#' @param basef  (\code{n x h(k* + m)}) matrix of base forecasts to be reconciled;
#' \code{n} is the total number of variables, \code{m} is the highest frequency,
#' \code{k*} is the sum of (\code{p-1}) factors of \code{m}, excluding \code{m},
#' and \code{h} is the forecast horizon. Each row identifies, a time series, and the forecasts
#' are ordered as [lowest_freq' ...  highest_freq']'.
#' @param hts_comb,thf_comb Type of covariance matrix (respectively (\code{n x n}) and
#' (\code{(k* + m) x (k* + m)})) to be used in the cross-sectional and temporal reconciliation,
#' see more in \code{comb} param of \code{\link[FoReco]{htsrec}} and \code{\link[FoReco]{thfrec}}.
#' @param res (\code{n x N(k* + m)}) matrix containing the residuals at all the
#' temporal frequencies, ordered [lowest_freq' ...  highest_freq']' (columns) for
#' each variable (row), needed to estimate the covariance matrix when \code{hts_comb =}
#' \code{\{"wls",} \code{"shr",} \code{"sam"\}} and/or \code{hts_comb =} \code{\{"wlsv",}
#' \code{"wlsh",} \code{"acov",} \code{"strar1",} \code{"sar1",} \code{"har1",}
#' \code{"shr",} \code{"sam"\}}. The rows must be in the same order as \code{basef}.
#' @param ... any other options useful for \code{\link[FoReco]{htsrec}} and
#' \code{\link[FoReco]{thfrec}}, e.g. \code{m}, \code{C} (or \code{Ut} and \code{nb}),
#' \code{nn} (for non negativity reconciliation only at first step), ...
#'
#' @return
#' The function returns a list with two elements:
#' \item{\code{recf}}{(\code{n x h(k* + m)}) reconciled forecasts matrix.}
#' \item{\code{M}}{Projection matrix (projection approach).}
#'
#' @references
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
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
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A., Boyd, S. (2018). OSQP:
#' An Operator Splitting Solver for Quadratic Programs, \href{https://arxiv.org/abs/1711.08013}{arXiv:1711.08013}.
#'
#' Stellato, B., Banjac, G., Goulart, P., Boyd, S., Anderson, E. (2019), OSQP:
#' Quadratic Programming Solver using the 'OSQP' Library, R package version 0.6.0.3
#' (October 10, 2019), \href{https://CRAN.R-project.org/package=osqp}{https://CRAN.R-project.org/package=osqp}.
#'
#' @keywords reconciliation heuristic
#' @examples
#' data(FoReco_data)
#' obj <- cstrec(FoReco_data$base, m = 12, C = FoReco_data$C, thf_comb = "acov",
#'               hts_comb = "shr", res = FoReco_data$res)
#'
#' @usage cstrec(basef, thf_comb, hts_comb, res, ...)
#'
#' @export
cstrec <- function(basef, thf_comb, hts_comb, res, ...) {

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
  m <- max(kset)
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
  arg_thf <- arg_thf[!(arg_thf %in% c("basef", "keep", "res", "", "comb", "m", "nn", "bounds"))]

  if (missing(res)) {
    M <- lapply(split(Y1, row(Y1)), function(x) {
      obj <- do.call("thfrec", c(list(basef = x, m = kset, comb = thf_comb),
                                 arg_input[which(names(arg_input) %in% arg_thf)]))
      out <- obj[["M"]]
      if(is.null(out)){
        return(obj[["G"]])
      }else{
        return(out)
      }
    })
  } else {
    M <- mapply(function(Y, X) {
      obj <- do.call("thfrec", c(list(basef = Y, m = kset, comb = thf_comb, res = X),
                                 arg_input[which(names(arg_input) %in% arg_thf)]))
      out <- obj[["M"]]
      if(is.null(out)){
        return(obj[["G"]])
      }else{
        return(out)
      }
    },
    Y = split(Y1, row(Y1)), X = split(res, row(res)), SIMPLIFY = FALSE)
  }

  ## Step 3: Cross-Temporal reconciled forecasts with heuristic
  meanM <- apply(simplify2array(lapply(M, as.matrix)), c(1, 2), sum) / length(M)

  if (h == 1) {
    Y3 <- Y1 %*% t(meanM)
  } else {
    h_pos <- rep(rep(1:h, length(kset)), rep((m/kset), each = h))
    D <- Dmat(h = h, kset = kset, n = 1)

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
