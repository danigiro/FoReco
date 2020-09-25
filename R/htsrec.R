#' @title Cross-sectional (contemporaneous) forecast reconciliation
#'
#' @description
#' Cross-sectional (contemporaneous) forecast reconciliation of hierarchical and
#' grouped time series. The reconciled forecasts are calculated either through a
#' projection approach (Byron, 1978, see also van Erven and Cugliari, 2015, and
#' Wickramasuriya et al., 2019), or the equivalent structural approach by Hyndman
#' et al. (2011). Moreover, the classic bottom-up approach is available.
#'
#' @usage htsrec(basef, comb, C, res, Ut, nb, mse = TRUE, corpcor = FALSE, W,
#'        type = "M", sol = "direct", nn = FALSE, keep = "list",
#'        settings = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'        eps_rel = 1e-5, polish_refine_iter = 100, polish=TRUE))
#'
#' @param basef (\code{h x n}) matrix of base forecasts to be reconciled;
#' \code{h} is the forecast horizon and \code{n} is the total number of time series.
#' @param comb Type of the reconciliation. Except for Bottom-up, each
#' option corresponds to a specified (\code{n x n}) covariance matrix:
#' \itemize{
#'     \item \bold{bu} (Bottom-up);
#'     \item \bold{ols} (Identity);
#'     \item \bold{struc} (Structural variances);
#'     \item \bold{wls} (Series variances) - uses res;
#'     \item \bold{shr} (Shrunk covariance matrix - MinT-shr) - uses res;
#'     \item \bold{sam} (Sample covariance matrix - MinT-sam) - uses res;
#'     \item \bold{w} use your personal matrix W in param \code{W}.
#'   }
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param res (\code{N x n}) in-sample residuals matrix needed when \code{comb =}
#' \code{\{"wls",} \code{"shr",} \code{"sam"\}}. The columns must be in
#' the same order as \code{basef}.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \eqn{(\textbf{U}'\textbf{Y} = \mathbf{0})}{} spanning the null space valid for the reconciled
#' forecasts. It can be used instead of parameter \code{C}, but in this case \code{nb} (n = na + nb) is needed. If
#' the hierarchy admits a structural representation, \code{Ut} has dimension (\code{na x n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb} is not used.
#' @param W This option permits to directly enter the covariance matrix:
#'   \enumerate{
#'     \item \code{W} must be a p.d. (\code{n x n}) matrix;
#'     \item if \code{comb} is different from "\code{w}", \code{W} is not used.
#'   }
#' @param mse Logical value: \code{TRUE} (\emph{default}) calculates the
#' covariance matrix of the in-sample residuals (when necessary) according to the original
#' \pkg{hts} and \pkg{thief} formulation: no mean correction, T as denominator.
#' @param corpcor Logical value: \code{TRUE} if \pkg{corpcor} (\enc{Schäfer}{Schafer} et
#' al., 2017) must be used to shrink the sample covariance matrix according to
#' \enc{Schäfer}{Schafer} and Strimmer (2005), otherwise the function uses the same
#' implementation as package \pkg{hts}.
#' @param type Approach used to compute the reconciled forecasts: \code{"M"} for
#' the projection approach with matrix M (\emph{default}), or \code{"S"} for the
#' structural approach with summing matrix S.
#' @param keep Return a list object of the reconciled forecasts at all levels.
#' @param sol Solution technique for the reconciliation problem: either \code{"direct"} (\emph{default}) for the direct
#' solution or \code{"osqp"} for the numerical solution (solving a linearly constrained quadratic
#' program using \code{\link[osqp]{solve_osqp}}).
#' @param nn Logical value: \code{TRUE} if non-negative reconciled forecasts are wished.
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}). The default options
#' are: \code{verbose = FALSE}, \code{eps_abs = 1e-5}, \code{eps_rel = 1e-5},
#' \code{polish_refine_iter = 100} and \code{polish = TRUE}. For details, see the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2019).
#'
#' @details
#' In case of non-negativity constraints, there are two ways:
#' \enumerate{
#'   \item \code{sol = "direct"} and \code{nn = TRUE}: the base forecasts
#'   will be reconciled at first without non-negativity constraints, then, if negative reconciled
#'   values are present, the \code{"osqp"} solver is used.
#'   \item \code{sol = "osqp"} and \code{nn = TRUE}: the base forecasts will be
#'   reconciled through the \code{"osqp"} solver.
#' }
#'
#' @return
#' If the parameter \code{keep} is equal to \code{"recf"}, then the function
#' returns only the (\code{h x n}) reconciled forecasts matrix, otherwise (\code{keep="all"})
#' it returns a list that mainly depends on what type of representation (\code{type})
#' and methodology (\code{sol}) have been used:
#' \item{\code{recf}}{(\code{h x n}) reconciled forecasts matrix.}
#' \item{\code{W}}{Covariance matrix used for forecast reconciliation.}
#' \item{\code{nn_check}}{Number of negative values (if zero, there are no values below zero).}
#' \item{\code{rec_check}}{Logical value: has the hierarchy been respected?}
#' \item{\code{M} (\code{type="M"} and \code{type="direct"})}{Projection matrix (projection approach)}
#' \item{\code{G} (\code{type="S"} and \code{type="direct"})}{Projection matrix (structural approach).}
#' \item{\code{S} (\code{type="S"} and \code{type="direct"})}{Cross-sectional summing matrix, \eqn{\textbf{S}}{S}.}
#' \item{\code{info} (\code{type="osqp"})}{matrix with information in columns
#' for each forecast horizon \code{h} (rows): run time (\code{run_time}) number of iteration,
#' norm of primal residual (\code{pri_res}), status of osqp's solution (\code{status}) and
#' polish's status (\code{status_polish}).}
#'
#' Only if \code{comb = "bu"}, the function returns \code{recf}, \code{S} and \code{M}.
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social accounts matrices,
#' \emph{Journal of the Royal Statistical Society A}, 141, 3, 359-367.
#'
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G., Shang, H.L.(2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
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
#' van Erven, T., Cugliari, J. (2015), Game-theoretically Optimal Reconciliation
#' of Contemporaneous Hierarchical Time Series Forecasts, in Antoniadis, A.,
#' Poggi, J.M., Brossat, X. (eds.), \emph{Modeling and Stochastic Learning for
#' Forecasting in High Dimensions}, Berlin, Springer, 297-317.
#'
#' Wickramasuriya, S.L., Athanasopoulos, G., Hyndman, R.J. (2019), Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace minimization,
#' \emph{Journal of the American Statistical Association}, 114, 526, 804-819.
#'
#' @keywords reconciliation
#' @examples
#' data(FoReco_data)
#' # monthly base forecasts
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$base),
#'                                     split = "_"))[1, ] == "k1")
#' mbase <- t(FoReco_data$base[, id])
#' # monthly residuals
#' id <- which(simplify2array(strsplit(colnames(FoReco_data$res),
#'                                     split = "_"))[1, ] == "k1")
#' mres <- t(FoReco_data$res[, id])
#' obj <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)
#'
#' @export
#'
#' @import Matrix osqp methods
#'
htsrec <- function(basef, comb, C, res, Ut, nb, mse = TRUE, corpcor = FALSE, W,
                   type = "M", sol = "direct", nn = FALSE, keep = "list",
                   settings = osqpSettings(
                     verbose = FALSE, eps_abs = 1e-5, eps_rel = 1e-5,
                     polish_refine_iter = 100, polish = TRUE
                   )) {
  if (missing(comb)) {
    stop("The argument comb is not specified.")
  }
  comb <- match.arg(comb, c("bu", "ols", "struc", "wls", "shr", "sam", "w"))

  type <- match.arg(type, c("M", "S"))
  keep <- match.arg(keep, c("list", "recf"))

  # base forecasts condition
  if (missing(basef)) {
    stop("The argument basef is not specified.")
  }
  if (NCOL(basef) == 1) {
    basef <- t(basef)
  }

  n <- NCOL(basef)

  # Using Ut AND nb or C
  if (missing(C)) {
    if (missing(Ut) | missing(nb)) {
      stop("Please, give C (or Ut AND nb).", call. = FALSE)
    }

    if (comb == "bu" | comb == "struc") {
      stop("bu and struc use C matrix.", call. = FALSE)
    }

    if (type == "S") {
      stop("Type = S use C matrix.", call. = FALSE)
    }

    if (n < 3) {
      stop("The hierarchy must have, at least, three series.", call. = FALSE)
    }

    if (NCOL(Ut) != n) {
      stop("Incorrect dimension of Ut or basef (they don't have same columns).", call. = FALSE)
    }

    if (n <= nb) {
      stop("n <= nb, total number of TS is less (or equal) than the number of bottom TS.", call. = FALSE)
    }

    na <- n - nb
  } else {
    if (!(is.matrix(C) | is(C, "Matrix"))) stop("C must be a matrix.", call. = FALSE)

    C <- Matrix(C, sparse = TRUE)

    nb <- NCOL(C)
    na <- NROW(C)

    if (comb != "bu" & nb != n) {
      if (n < 3) {
        stop("The hierarchy must have, at least, three series.", call. = FALSE)
      }

      if (n <= nb) {
        stop("n <= nb, total number of TS is less (or equal) than the number of bottom TS.", call. = FALSE)
      }

      if (na + nb != n) {
        stop("na + nb != n, matrix C or basef has incorrect dimension.", call. = FALSE)
      }
    }

    S <- rbind(C, .sparseDiagonal(nb))
    Ut <- cbind(.sparseDiagonal(na), -C)
  }

  if (any(comb == c("wls", "shr", "sam"))) {
    if (missing(res)) {
      stop("Don't forget residuals!", call. = FALSE)
    } else if (NCOL(res) != n) {
      stop(paste("The number of columns of res must be", n), call. = FALSE)
    } else if (comb == "sam" & NROW(res) < n) {
      stop("The number of rows of res is less than columns: \n it could lead to singularity problems if comb == 'sam'.", call. = FALSE)
    }
  }

  if (mse) {
    cov_mod <- function(x, ...) crossprod(stats::na.omit(x), ...) / NROW(stats::na.omit(x))
  } else {
    cov_mod <- function(x, ...) stats::var(x, na.rm = TRUE, ...)
  }

  if (corpcor) {
    shr_mod <- function(x, ...) corpcor::cov.shrink(x, verbose = FALSE, ...)
  } else {
    shr_mod <- function(x, ...) shrink_estim(x, minT = mse)[[1]]
  }

  switch(comb,
    bu = {
      if (n == nb) {
        outf <- basef %*% t(S)
      } else {
        outf <- basef[, (na + 1):n] %*% t(S)
      }

      rownames(outf) <- paste("h", 1:NROW(outf), sep = "")
      colnames(outf) <- if (is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)

      if (keep == "list") {
        return(list(recf = outf, S = S, M = S %*% cbind(matrix(0, nb, na), diag(nb))))
      } else {
        return(outf)
      }
    },
    ols =
      W <- .sparseDiagonal(n),
    struc =
      W <- .sparseDiagonal(x = rowSums(S)),
    wls = {
      diagW <- diag(cov_mod(res))
      W <- .sparseDiagonal(x = diagW)
    },
    shr = {
      W <- shr_mod(res)
    },
    sam = {
      W <- cov_mod(res)
    },
    w = {
      if (missing(W)) {
        stop("Please, put in option W your covariance matrix", call. = FALSE)
      }
      W <- W
    }
  )

  b_pos <- c(rep(0, na), rep(1, nb))

  if (type == "S") {
    rec_sol <- recoS(
      basef = basef, W = W, S = S, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos
    )
  } else {
    rec_sol <- recoM(
      basef = basef, W = W, H = Ut, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos
    )
  }

  if (keep == "list") {
    rec_sol$nn_check <- sum(rec_sol$recf < 0)
    rec_sol$rec_check <- all(rec_sol$recf %*% t(Ut) < 1e-6)

    rec_sol$recf <- as.matrix(rec_sol$recf)
    rownames(rec_sol$recf) <- paste("h", 1:NROW(rec_sol$recf), sep = "")
    colnames(rec_sol$recf) <- if (is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)

    rec_sol$W <- as.matrix(rec_sol$W)
    dimnames(rec_sol$W) <- NULL

    names_all_list <- c("recf", "W", "Omega", "nn_check", "rec_check", "varf", "M", "G", "S", "info")
    names_list <- names(rec_sol)
    rec_sol <- rec_sol[names_list[order(match(names_list, names_all_list))]]
    return(rec_sol)
  } else {
    if (length(rec_sol) == 1) {
      rec_sol$recf <- as.matrix(rec_sol$recf)
      rownames(rec_sol$recf) <- paste("h", 1:NROW(rec_sol$recf), sep = "")
      colnames(rec_sol$recf) <- if (is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
      return(rec_sol$recf)
    } else {
      rec_sol$recf <- as.matrix(rec_sol$recf)
      rownames(rec_sol$recf) <- paste("h", 1:NROW(rec_sol$recf), sep = "")
      colnames(rec_sol$recf) <- if (is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)
      return(rec_sol)
    }
  }
}
