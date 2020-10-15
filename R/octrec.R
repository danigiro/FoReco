#' @title Optimal combination cross-temporal forecast reconciliation
#'
#' @description
#' Optimal (in least squares sense) combination cross-temporal forecast reconciliation.
#' The reconciled forecasts are calculated either through a projection approach
#' (Byron, 1978), or the equivalent structural approach by Hyndman et al. (2011).
#'
#' @usage octrec(basef, m, C, comb, res, Ut, nb, W, Sstruc,
#'        mse = TRUE, corpcor = FALSE, type = "M", sol = "direct",
#'        nn = FALSE, keep = "list",
#'        settings = osqpSettings(verbose = FALSE, eps_abs = 1e-5,
#'        eps_rel = 1e-5, polish_refine_iter = 100, polish=TRUE))
#'
#' @param basef  (\code{n x h(k* + m)}) matrix of base forecasts to be reconciled;
#' \code{n} is the total number of variables, \code{m} is the highest frequency,
#' \code{k*} is the sum of (\code{p-1}) factors of \code{m}, excluding \code{m},
#' and \code{h} is the forecast horizon. Each row identifies, a time series, and the forecasts
#' are ordered as [lowest_freq' ...  highest_freq']'.
#' @param m Highest available sampling frequency per seasonal cycle (max. order of temporal aggregation).
#' @param comb Type of the reconciliation, it corrispond to a different covariance
#' matrix (\code{n(k* + m) x n(k* + m)}), \code{k*} is the sum of (\code{p-1})
#' factors of \code{m} (exclude \code{m} as factors) and \code{n} is the number
#' of variables:
#' \itemize{
#'   \item \bold{ols} (Identity);
#'   \item \bold{struc} (Cross-temporal summing matrix, use the \code{Sstruc} param to reduce computation time);
#'   \item \bold{wlsh} (Hierarchy variances matrix);
#'   \item \bold{wlsv} (Series variances matrix);
#'   \item \bold{bdshr} (Shrunk cross-covariance matrix, cross-sectional framework);
#'   \item \bold{bdsam}  (Sample cross-covariance matrix, cross-sectional framework);
#'   \item \bold{acov} (Series auto-covariance matrix);
#'   \item \bold{Sshr} (Series shrunk cross-covariance matrix);
#'   \item \bold{Ssam} (Series cross-covariance matrix);
#'   \item \bold{shr}  (Shrunk cross-covariance matrix);
#'   \item \bold{sam} (Sample cross-covariance matrix);
#'   \item \bold{w} use your personal matrix W in param \code{W}.
#' }
#' @param C (\code{na x nb}) cross-sectional (contemporaneous) matrix mapping the bottom
#' level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \eqn{(\textbf{U}'\textbf{Y} = \mathbf{0})}{} spanning the null space valid for the reconciled
#' forecasts. It can be used instead of parameter \code{C}, but in this case \code{nb} (n = na + nb) is needed. If
#' the hierarchy admits a structural representation, \code{Ut} has dimension (\code{na x n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb} is not used.
#' @param res (\code{n x N(k* + m)}) matrix containing the residuals at all the
#' temporal frequencies ordered [lowest_freq' ...  highest_freq']' (columns) for
#' each variable (row), needed to estimate the covariance matrix when \code{comb =}
#'  \code{\{"sam",} \code{"wlsv",} \code{"wlsh",} \code{"acov",} \code{"Ssam",}
#'  \code{"Sshr",} \code{"Sshr1",} \code{"shr"\}}.
#' @param W This option permits to directly enter the covariance matrix:
#' \enumerate{
#'   \item \code{W} must be a p.d. (\code{n(k* + m) x n(k* + m)}) matrix;
#'   \item if \code{comb} is different from "\code{w}", \code{W} is not used.
#' }
#' @param Sstruc Cross-temporal summing matrix (structural representation)\eqn{,\; \check{\textbf{S}}}{};
#' can be obtained through the function \link[FoReco]{ctf_tools}.
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
#' returns only the (\code{n x h(k* + m)}) reconciled forecasts matrix, otherwise (\code{keep="all"})
#' it returns a list that mainly depends on what type of representation (\code{type})
#' and methodology (\code{sol}) have been used:
#' \item{\code{recf}}{(\code{n x h(k* + m)}) reconciled forecasts matrix.}
#' \item{\code{Omega}}{Covariance matrix used for reconciled forecasts (\code{vec}(\code{Y}\enc{'}{'}) representation).}
#' \item{\code{W}}{Covariance matrix used for reconciled forecasts (\code{vec(Y)} representation).}
#' \item{\code{nn_check}}{Number of negative values (if zero, there are no values below zero).}
#' \item{\code{rec_check}}{Logical value: has the hierarchy been respected?}
#' \item{\code{M} (\code{type="M"} and \code{type="direct"})}{Projection matrix (projection approach).}
#' \item{\code{G} (\code{type="S"} and \code{type="direct"})}{Projection matrix (structural approach).}
#' \item{\code{S} (\code{type="S"} and \code{type="direct"})}{Cross-temporal summing matrix\eqn{, \; \textbf{Q}\check{\textbf{S}}}{} (\code{vec}(\code{Y}\enc{'}{'}) representation).}
#' \item{\code{info} (\code{type="osqp"})}{matrix with some useful indicators (columns)
#' for each forecast horizon \code{h} (rows): run time (\code{run_time}) number of iteration,
#' norm of primal residual (\code{pri_res}), status of osqp's solution (\code{status}) and
#' polish's status (\code{status_polish}).}
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social accounts matrices,
#' \emph{Journal of the Royal Statistical Society A}, 141, 3, 359-367.
#'
#' Di Fonzo, T., Girolimetto, D. (2020), Cross-Temporal Forecast Reconciliation:
#' Optimal Combination Method and Heuristic Alternatives, Department of Statistical
#' Sciences, University of Padua, \href{https://arxiv.org/abs/2006.08570}{arXiv:2006.08570}.
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
#' @keywords reconciliation
#' @examples
#' data(FoReco_data)
#' obj <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'               comb = "bdshr", res = FoReco_data$res)
#'
#' @export
#'
#' @import Matrix osqp methods
#'
octrec <- function(basef, m, C, comb, res, Ut, nb, W, Sstruc, mse = TRUE,
                   corpcor = FALSE, type = "M", sol = "direct", nn = FALSE, keep = "list",
                   settings = osqpSettings(
                     verbose = FALSE, eps_abs = 1e-5, eps_rel = 1e-5,
                     polish_refine_iter = 100, polish = TRUE
                   )) {
  if (missing(m)) {
    stop("The argument m is not specified")
  }

  if (missing(comb)) {
    stop("The argument comb is not specified")
  }
  comb <- match.arg(comb, c(
    "ols", "struc", "sam", "wlsv", "wlsh", "shr",
    "acov", "Ssam", "Sshr", "bdshr", "bdsam", "w"
  ))

  type <- match.arg(type, c("M", "S"))
  keep <- match.arg(keep, c("list", "recf"))

  if (missing(basef)) {
    stop("The argument basef is not specified")
  }

  tools <- thf_tools(m, sparse = TRUE)
  kset <- tools$kset
  p <- tools$p
  kt <- tools$kt
  ks <- tools$ks

  # matrix
  Zt <- tools$Zt

  n <- NROW(basef)
  if (missing(C)) { # Using Ut AND nb or C
    if (missing(Ut) | missing(nb)) {
      stop("Please, give C (or Ut AND nb)", call. = FALSE)
    }

    if (comb == "struc") {
      stop("struc use C matrix", call. = FALSE)
    }

    if (n < 3) {
      stop("The hierarchy must have, at least, three series", call. = FALSE)
    }

    if (NCOL(Ut) != n) {
      stop("Incorrect dimension of Ut or basef (they don't have same columns)", call. = FALSE)
    }

    if (n <= nb) {
      stop("n <= nb, total number of TS is less (or equal) than the number of bottom TS", call. = FALSE)
    }

    na <- n - nb
    Ut <- Matrix(Ut, sparse = TRUE)
  } else {
    if (!(is.matrix(C) | is(C, "Matrix"))) stop("C must be a matrix", call. = FALSE)

    C <- Matrix(C, sparse = TRUE)

    nb <- NCOL(C)
    na <- NROW(C)

    if (n < 3) {
      stop("The hierarchy must have, at least, three series", call. = FALSE)
    }

    if (n <= nb) {
      stop("n <= nb, total number of TS is less (or equal) than the number of bottom TS", call. = FALSE)
    }

    if (na + nb != n) {
      stop("na + nb != n, matrix C or basef has incorrect dimension", call. = FALSE)
    }

    Ut <- cbind(Diagonal(na), -C)
  }

  # H
  P <- commat(n, kt)
  Us <- cbind(Matrix(0, NROW(Ut) * m, n * ks), kronecker(Diagonal(m), Ut)) %*% t(P)
  Ht <- rbind(Us, kronecker(Diagonal(n), Zt))

  # Base forecast
  if (NCOL(basef) %% kt != 0) {
    stop("basef has a number of row not in line with the frequency of the series", call. = FALSE)
  }

  h <- NCOL(basef) / kt
  Dh <- Dmat(h = h, kset = kset, n = n)
  Ybase <- matrix(Dh %*% as.vector(t(basef)), nrow = h, byrow = TRUE)

  # In-sample errors
  if (any(comb == c("sam", "wlsv", "wlsh", "shr", "acov", "Ssam", "Sshr", "bdshr", "bdsam"))) {
    # residual condition
    if (missing(res)) {
      stop("Don't forget residuals!", call. = FALSE)
    }
    if (NCOL(res) == 1) {
      stop("res must be a matrix", call. = FALSE)
    }
    if (NCOL(res) %% kt != 0) {
      stop("res has a number of columns not in line with frequency of the series", call. = FALSE)
    }

    if (NCOL(res) %% kt != 0) {
      stop("res has a number of columns not in line with frequency of the series", call. = FALSE)
    }
    if (NROW(res) != n) {
      stop(paste("The number of rows of res must be", n), call. = FALSE)
    }

    N <- NCOL(res) / kt

    DN <- Dmat(h = N, kset = kset, n = n)
    E <- matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)

    if (comb == "sam" & N < n * kt) {
      stop("N < n(k* + m): it could lead to singularity problems if comb == sam", call. = FALSE)
    }
    if (comb == "acov" & N < m) {
      stop("N < m: it could lead to singularity problems if comb == acov", call. = FALSE)
    }
    if (comb == "Ssam" & N < kt) {
      stop("N < (k* + m): it could lead to singularity problems if comb == Ssam", call. = FALSE)
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

  if (type == "S" | comb == "struc") {
    Qtilde <- commat(nb, kt) %*% bdiag(commat(ks, nb), commat(m, nb))
    Q <- bdiag(Diagonal(na * kt), Qtilde)
    if (missing(Sstruc)) {
      Htstruc <- Matrix(pracma::rref(as.matrix(Ht %*% Q)), sparse = TRUE)
      Cstruc <- -Htstruc[, (nrow(Htstruc) + 1):ncol(Htstruc)]
      Sstruc <- rbind(Cstruc, Diagonal(m * nb))
    }
    S <- Q %*% Sstruc
  }

  switch(comb,
    ols = {
      Omega <- .sparseDiagonal(n * kt)
    },
    struc = {
      Omega <- .sparseDiagonal(x = rowSums(S))
    },
    sam = {
      Omega <- cov_mod(E)
    },
    wlsh = {
      Omega <- .sparseDiagonal(x = diag(cov_mod(E)))
    },
    wlsv = {
      var_freq <- apply(res, 1, function(z) sapply(kset, function(x) cov_mod(z[rep(kset, rev(kset) * N) == x])))
      Omega <- .sparseDiagonal(x = rep(as.vector(var_freq), rep(rev(kset), n)))
    },
    shr = {
      Omega <- shr_mod(E)
    },
    acov = {
      mat1 <- bdiag(rep(lapply(rev(kset), function(x) matrix(1, nrow = x, ncol = x)), n))
      cov <- cov_mod(E)
      Omega <- cov * mat1
    },
    Ssam = {
      mat1 <- bdiag(replicate(n, matrix(1, nrow = kt, ncol = kt), simplify = FALSE))
      cov <- cov_mod(E)
      Omega <- cov * mat1
    },
    Sshr = {
      shrink <- lapply(1:n, function(x) shr_mod(E[, rep(1:n, each = kt) == x]))
      Omega <- bdiag(shrink)
    },
    w = {
      if (missing(W)) {
        stop("Please, put in option W your covariance matrix", call. = FALSE)
      }
      Omega <- P %*% W %*% t(P)
    },
    bdshr = {
      blockW <- lapply(kset, function(x) shr_mod(t(res[, rep(kset, N * rev(kset)) == x])))
      blockW <- rep(blockW, rev(kset))
      P <- commat(n, kt)
      Omega <- P %*% bdiag(blockW) %*% t(P)
    },
    bdsam = {
      blockW <- lapply(kset, function(x) cov_mod(t(res[, rep(kset, N * rev(kset)) == x])))
      blockW <- rep(blockW, rev(kset))
      P <- commat(n, kt)
      Omega <- P %*% bdiag(blockW) %*% t(P)
    }
  )

  b_pos <- c(rep(0, na * kt), rep(rep(kset, rev(kset)), nb) == 1)

  if (type == "S") {
    rec_sol <- recoS(
      basef = Ybase, W = Omega, S = S, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos
    )
  } else {
    rec_sol <- recoM(
      basef = Ybase, W = Omega, H = Ht, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos
    )
  }

  if (keep == "list") {
    rec_sol$nn_check <- sum(rec_sol$recf < 0)
    rec_sol$rec_check <- all(rec_sol$recf %*% t(Ht) < 1e-6)

    rec_sol$recf <- matrix(t(Dh) %*% as.vector(t(rec_sol$recf)), nrow = n, byrow = TRUE)
    names(rec_sol)[names(rec_sol) == "W"] <- "Omega"
    rec_sol$W <- as.matrix(t(P) %*% rec_sol$Omega %*% P)
    dimnames(rec_sol$W) <- NULL
    rec_sol$Omega <- as.matrix(rec_sol$Omega)
    dimnames(rec_sol$Omega) <- NULL

    names_all_list <- c("recf", "W", "Omega", "nn_check", "rec_check", "varf", "M", "G", "S", "info")
    names_list <- names(rec_sol)
    rec_sol <- rec_sol[names_list[order(match(names_list, names_all_list))]]

    rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
    colnames(rec_sol$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
      do.call("c", as.list(sapply(
        rev(kset) * h,
        function(x) seq(1:x)
      ))),
      sep = ""
    )
    return(rec_sol)
  } else {
    if (length(rec_sol) == 1) {
      rec_sol$recf <- matrix(t(Dh) %*% as.vector(t(rec_sol$recf)), nrow = n, byrow = TRUE)
      rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
      colnames(rec_sol$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
        do.call("c", as.list(sapply(
          rev(kset) * h,
          function(x) seq(1:x)
        ))),
        sep = ""
      )
      return(rec_sol$recf)
    } else {
      rec_sol$recf <- matrix(t(Dh) %*% as.vector(t(rec_sol$recf)), nrow = n, byrow = TRUE)
      rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
      colnames(rec_sol$recf) <- paste("k", rep(kset, h * rev(kset)), "h",
        do.call("c", as.list(sapply(
          rev(kset) * h,
          function(x) seq(1:x)
        ))),
        sep = ""
      )
      return(rec_sol)
    }
  }
}
