#' @title Forecast reconciliation through temporal hierarchies (temporal reconciliation)
#'
#' @description
#' Forecast reconciliation of one time series through temporal hierarchies
#' (Athanasopoulos et al., 2017). The reconciled forecasts are calculated
#' either through a projection approach (Byron, 1978), or the equivalent
#' structural approach by Hyndman et al. (2011). Moreover, the classic
#' bottom-up approach is available.
#'
#' @usage thfrec(basef, m, comb, res, mse = TRUE, corpcor = FALSE,
#'        type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
#'         nn_type = "osqp", settings = list(), bounds = NULL, Omega = NULL)
#'
#' @param basef (\mjseqn{h(k^\ast + m) \times 1}) vector of base forecasts to be
#' reconciled, containing the forecasts at all the needed temporal frequencies
#' ordered as [lowest_freq' ...  highest_freq']'.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param comb Type of the reconciliation. Except for bottom up, all other
#' options correspond to a different (\mjseqn{(k^\ast + m) \times (k^\ast + m)})
#' covariance matrix, \mjseqn{k^\ast} is the sum of (\mjseqn{p-1}) factors of
#' \mjseqn{m} (excluding \mjseqn{m}):
#' \itemize{
#'   \item \bold{bu} (Bottom-up);
#'   \item \bold{ols} (Identity);
#'   \item \bold{struc} (Structural variances);
#'   \item \bold{wlsv} (Series variances);
#'   \item \bold{wlsh} (Hierarchy variances);
#'   \item \bold{acov} (Auto-covariance matrix);
#'   \item \bold{strar1} (Structural Markov);
#'   \item \bold{sar1} (Series Markov);
#'   \item \bold{har1} (Hierarchy Markov);
#'   \item \bold{shr} (Shrunk cross-covariance matrix);
#'   \item \bold{sam} (Sample cross-covariance matrix);
#'   \item \bold{omega} use your personal matrix Omega in param \code{Omega}.
#' }
#' @param res vector containing the in-sample residuals at all the temporal
#' frequencies ordered as \code{basef}, i.e. [lowest_freq' ...  highest_freq']',
#' needed to estimate the covariance matrix when \code{comb =} \code{\{"wlsv",}
#' \code{"wlsh",} \code{"acov",} \code{"strar1",} \code{"sar1",} \code{"har1",}
#' \code{"shr",} \code{"sam"\}}.
#' @param Omega This option permits to directly enter the covariance matrix:
#' \enumerate{
#'   \item \code{Omega} must be a p.d. (\mjseqn{(k^\ast + m) \times (k^\ast + m)})
#'   matrix or a list of \mjseqn{h} matrix (one for each forecast horizon);
#'   \item if \code{comb} is different from "\code{omega}", \code{Omega} is
#'   not used.
#' }
#' @param mse Logical value: \code{TRUE} (\emph{default}) calculates the
#' covariance matrix of the in-sample residuals (when necessary) according to
#' the original \pkg{hts} and \pkg{thief} formulation: no mean correction,
#' T as denominator.
#' @param corpcor Logical value: \code{TRUE} if \pkg{corpcor} (\enc{Schäfer}{Schafer} et
#' al., 2017) must be used to shrink the sample covariance matrix according to
#' \enc{Schäfer}{Schafer} and Strimmer (2005), otherwise the function uses the
#' same implementation as package \pkg{hts}.
#' @param type Approach used to compute the reconciled forecasts: \code{"M"} for
#' the projection approach with matrix M (\emph{default}), or \code{"S"} for the
#' structural approach with temporal summing matrix R.
#' @param keep Return a list object of the reconciled forecasts at all levels
#' (if keep = "list") or only the reconciled forecasts matrix (if keep = "recf").
#' @param sol Solution technique for the reconciliation problem: either
#' \code{"direct"} (\emph{default}) for the closed-form matrix solution, or
#' \code{"osqp"} for the numerical solution (solving a linearly constrained
#' quadratic program using \code{\link[osqp]{solve_osqp}}).
#' @param nn Logical value: \code{TRUE} if non-negative reconciled forecasts
#' are wished.
#' @param nn_type "osqp" (default), "KAnn" (only \code{type == "M"}) or "sntz".
#' @param settings Settings for \pkg{osqp} (object \code{\link[osqp]{osqpSettings}}).
#' The default options are: \code{verbose = FALSE}, \code{eps_abs = 1e-5},
#' \code{eps_rel = 1e-5}, \code{polish_refine_iter = 100} and \code{polish = TRUE}.
#' For details, see the \href{https://osqp.org/}{\pkg{osqp} documentation}
#' (Stellato et al., 2019).
#' @param bounds (\mjseqn{(k^\ast + m) \times 2}) matrix with temporal bounds: the
#' first column is the lower bound, and the second column is the upper bound.
#' @param v vector index of the fixed base forecast (\mjseqn{\mbox{min}(v) > 0}
#' and \mjseqn{\mbox{max}(v) < (k^\ast + m)}).
#'
#' @details
#' \loadmathjax
#' Let \mjseqn{m} be the highest available
#' sampling frequency per seasonal cycle, and denote
#' \mjseqn{{\cal K} = \left\lbrace k_m, k_{p-1}, \ldots, k_{2}, k_1\right\rbrace}
#' the \mjseqn{p} factors of \mjseqn{m}, in descending order, where
#' \mjseqn{k_p=m}, and \mjseqn{k_1=1}. Define \mjseqn{\mathbf{K}}
#' the \mjseqn{\left( k^\ast \times m\right)} temporal aggregation matrix
#' converting the high-frequency observations into lower-frequency
#' (temporally aggregated) ones:
#' \mjsdeqn{ \mathbf{K} = \left[\begin{array}{c}
#' \mathbf{1}_m' \cr
#' \mathbf{I}_{\frac{m}{k_{p-1}}} \otimes \mathbf{1}_{k_{p-1}}' \cr
#' \vdots \cr
#' \mathbf{I}_{\frac{m}{k_{2}}} \otimes \mathbf{1}_{k_{2}}' \cr
#' \end{array}\right].}
#' Denote
#' \mjseqn{\mathbf{R} = \left[\begin{array}{c}
#' \mathbf{K} \cr \mathbf{I}_m
#' \end{array}\right]} the \mjseqn{\left[(k^\ast+m) \times m \right]}
#' \emph{temporal summing} matrix, and
#' \mjseqn{\mathbf{Z}' = \left[ \mathbf{I}_{k^\ast} \; -\mathbf{K} \right]}
#' the zero constraints kernel matrix.
#'
#' Suppose we have the \mjseqn{\left[(k^\ast+m) \times 1\right]} vector
#' \mjseqn{\widehat{\mathbf{y}}} of unbiased base forecasts for the
#' \mjseqn{p} temporal aggregates of a single time series \mjseqn{Y}
#' within a complete time cycle, i.e. at the forecast horizon \mjseqn{h=1}
#' for the lowest (most aggregated) time frequency. If the base forecasts
#' have been independently obtained, generally they do not fulfill the
#' temporal aggregation constraints, i.e. \mjseqn{\mathbf{Z}'
#' \widehat{\mathbf{y}} \ne \mathbf{0}_{(k^\ast \times 1)}}.
#' By adapting the general point forecast reconciliation according to
#' the projection approach (\code{type = "M"}),
#' the vector of temporally reconciled forecasts
#' is given by:
#' \mjsdeqn{\widetilde{\mathbf{y}} = \widehat{\mathbf{y}} -
#' \mathbf{\Omega}\mathbf{Z}\left(\mathbf{Z}'\mathbf{\Omega}
#' \mathbf{Z}\right)^{-1}\mathbf{Z}'\widehat{\mathbf{y}},}
#' where \mjseqn{\mathbf{\Omega}} is a \mjseqn{\left[(k^\ast+m)
#' \times (k^\ast+m)\right]} p.d. matrix, assumed known. The alternative
#' equivalent solution (\code{type = "S"}) following the
#' structural reconciliation approach by Athanasopoulos et al. (2017) is given by:
#' \mjsdeqn{\widetilde{\mathbf{y}} = \mathbf{R}\left(\mathbf{R}'
#' \mathbf{\Omega}^{-1}\mathbf{R}\right)^{-1}\mathbf{R}'
#' \mathbf{\Omega}^{-1}\widehat{\mathbf{y}}.}
#'
#' \strong{Bounds on the reconciled forecasts}
#'
#' When the reconciliation makes use of the optimization package osqp,
#' the user may impose bounds on the reconciled forecasts.
#' The parameter \code{bounds} permits to consider lower (\mjseqn{\mathbf{a}}) and
#' upper (\mjseqn{\mathbf{b}}) bounds like \mjseqn{\mathbf{a} \leq
#' \widetilde{\mathbf{y}} \leq \mathbf{b}} such that:
#' \mjsdeqn{ \begin{array}{c}
#' a_1 \leq \widetilde{y}_1 \leq b_1 \cr
#' \dots \cr
#' a_{(k^\ast + m)} \leq \widetilde{y}_{(k^\ast + m)} \leq b_{(k^\ast + m)} \cr
#' \end{array} \Rightarrow
#' \mbox{bounds} = [\mathbf{a} \; \mathbf{b}] =
#' \left[\begin{array}{cc}
#' a_1 & b_1 \cr
#' \vdots & \vdots \cr
#' a_{(k^\ast + m)} & b_{(k^\ast + m)} \cr
#' \end{array}\right],}
#' where \mjseqn{a_i \in [- \infty, + \infty]} and \mjseqn{b_i \in [- \infty, + \infty]}.
#' If \mjseqn{y_i} is unbounded, the \mjseqn{i}-th row of \code{bounds} would be equal
#' to \code{c(-Inf, +Inf)}.
#' Notice that if the \code{bounds} parameter is used, \code{sol = "osqp"} must be used.
#' This is not true in the case of non-negativity constraints:
#' \itemize{
#'   \item \code{sol = "direct"}: first the base forecasts
#'   are reconciled without non-negativity constraints, then, if negative reconciled
#'   values are present, the \code{"osqp"} solver is used;
#'   \item \code{sol = "osqp"}: the base forecasts are
#'   reconciled using the \code{"osqp"} solver.
#' }
#'
#' In this case it is not necessary to build a matrix containing
#' the bounds, and it is sufficient to set \code{nn = "TRUE"}.
#'
#' Non-negative reconciled forecasts may be obtained by setting \code{nn_type} alternatively as:
#' \itemize{
#'   \item \code{nn_type = "KAnn"} (Kourentzes and Athanasopoulos, 2021)
#'   \item \code{nn_type = "sntz"} ("set-negative-to-zero")
#'   \item \code{nn_type = "osqp"} (Stellato et al., 2020)
#' }
#'
#' @return
#' If the parameter \code{keep} is equal to \code{"recf"}, then the function
#' returns only the (\mjseqn{h(k^\ast + m) \times 1}) reconciled forecasts vector, otherwise (\code{keep="all"})
#' it returns a list that mainly depends on what type of representation (\code{type})
#' and solution technique (\code{sol}) have been used:
#' \item{\code{recf}}{(\mjseqn{h(k^\ast + m) \times 1}) reconciled forecasts vector, \mjseqn{\widetilde{\mathbf{y}}}.}
#' \item{\code{Omega}}{Covariance matrix used for forecast reconciliation, \mjseqn{\mathbf{\Omega}}.}
#' \item{\code{nn_check}}{Number of negative values (if zero, there are no values below zero).}
#' \item{\code{rec_check}}{Logical value: has the hierarchy been respected?}
#' \item{\code{varf} (\code{type="direct"})}{(\mjseqn{(k^\ast + m) \times 1}) reconciled forecasts variance vector for \mjseqn{h=1}, \mjseqn{\mbox{diag}(\mathbf{MW}}).}
#' \item{\code{M} (\code{type="direct"})}{Projection matrix, \mjseqn{\mathbf{M}} (projection approach).}
#' \item{\code{G} (\code{type="S"} and \code{type="direct"})}{Projection matrix, \mjseqn{\mathbf{G}} (structural approach, \mjseqn{\mathbf{M}=\mathbf{R}\mathbf{G}}).}
#' \item{\code{S} (\code{type="S"} and \code{type="direct"})}{Temporal summing matrix, \mjseqn{\mathbf{R}}.}
#' \item{\code{info} (\code{type="osqp"})}{matrix with information in columns
#' for each forecast horizon \mjseqn{h} (rows): run time (\code{run_time}),
#' number of iteration (\code{iter}), norm of primal residual (\code{pri_res}),
#' status of osqp's solution (\code{status}) and polish's status
#' (\code{status_polish}). It will also be returned with \code{nn = TRUE} if
#' a solver (see \code{nn_type}) will be use.}
#'
#' Only if \code{comb = "bu"}, the function returns \code{recf}, \code{R} and \code{M}.
#'
#' @references
#' Athanasopoulos, G., Hyndman, R.J., Kourentzes, N., Petropoulos, F. (2017),
#' Forecasting with Temporal Hierarchies, \emph{European Journal of Operational
#' Research}, 262, 1, 60-74.
#'
#' Byron, R.P. (1978), The estimation of large social accounts matrices,
#' \emph{Journal of the Royal Statistical Society A}, 141, 3, 359-367.
#'
#' Di Fonzo, T., and Girolimetto, D. (2021), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, in press.
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G., Shang, H.L. (2011), Optimal combination
#' forecasts for hierarchical time series, \emph{Computational Statistics & Data
#' Analysis}, 55, 9, 2579-2589.
#'
#' Nystrup, P.,  \enc{Lindström}{Lindstrom}, E., Pinson, P., Madsen, H. (2020),
#' Temporal hierarchies with autocorrelation for load forecasting,
#' \emph{European Journal of Operational Research}, 280, 1, 876-888.
#'
#' \enc{Schäfer}{Schafer}, J.L., Opgen-Rhein, R., Zuber, V., Ahdesmaki, M.,
#' Duarte Silva, A.P., Strimmer, K. (2017), Package `corpcor', R
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
#' @keywords bottom-up
#' @family reconciliation procedures
#' @examples
#' data(FoReco_data)
#' # top ts base forecasts ([lowest_freq' ...  highest_freq']')
#' topbase <- FoReco_data$base[1, ]
#'  # top ts residuals ([lowest_freq' ...  highest_freq']')
#' topres <- FoReco_data$res[1, ]
#' obj <- thfrec(topbase, m = 12, comb = "acov", res = topres)
#'
#' @export
#'
#' @import Matrix osqp
#'
thfrec <- function(basef, m, comb, res, mse = TRUE, corpcor = FALSE,
                   type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
                   nn_type = "osqp", settings = list(), bounds = NULL, Omega = NULL) {

  if(missing(comb)){
    stop("The argument comb is not specified.", call. = FALSE)
  }else if(comb != "omega"){
    Omega <-  NULL
  }

  UseMethod("thfrec", Omega)
}

#' @export
thfrec.default <- function(basef, m, comb, res, mse = TRUE, corpcor = FALSE,
                           type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
                           nn_type = "osqp", settings = list(), bounds = NULL, Omega) {
  # m condition
  if (missing(m)) {
    stop("The argument m is not specified", call. = FALSE)
  }
  tools <- thf_tools(m)
  kset <- tools$kset
  m <- tools$m
  p <- tools$p
  kt <- tools$kt
  ks <- tools$ks

  # matrix
  K <- tools$K
  R <- tools$R
  Zt <- tools$Zt

  comb <- match.arg(comb, c(
    "bu", "ols", "struc", "wlsv", "wlsh", "acov",
    "strar1", "sar1", "har1", "shr", "sam", "omega"
  ))
  type <- match.arg(type, c("M", "S"))
  keep <- match.arg(keep, c("list", "recf"))

  nn_type <- match.arg(nn_type, c("osqp", "KAnn", "fbpp", "sntz"))

  if(nn){
    if(nn_type == "fbpp" | nn_type == "KAnn"){
      type = "M"
    }
  }

  # base forecasts condition
  if (missing(basef)) {
    stop("The argument basef is not specified", call. = FALSE)
  }

  if (NCOL(basef) != 1) {
    stop("basef must be a vector", call. = FALSE)
  }

  # Base Forecasts matrix
  if (comb == "bu" & length(basef) %% m == 0) {
    h <- length(basef) / m
    Dh <- Dmat(h = h, m = kset, n = 1)
    BASEF <- matrix(basef, h, m, byrow = T)
  } else if (length(basef) %% kt != 0) {
    stop("basef vector has a number of elements not in line with the frequency of the series", call. = FALSE)
  } else {
    h <- length(basef) / kt
    Dh <- Dmat(h = h, m = kset, n = 1)
    BASEF <- matrix(Dh %*% basef, nrow = h, byrow = T)
  }

  # Residuals Matrix
  if (any(comb == c("wlsv", "wlsh", "acov", "strar1", "sar1", "har1", "sGlasso", "hGlasso", "shr", "sam"))) {
    # residual condition
    if (missing(res)) {
      stop("Don't forget residuals!", call. = FALSE)
    }
    if (NCOL(res) != 1) {
      stop("res must be a vector", call. = FALSE)
    }
    if (length(res) %% kt != 0) {
      stop("res vector has a number of row not in line with frequency of the series", call. = FALSE)
    }

    N <- length(res) / kt
    DN <- Dmat(h = N, m = kset, n = 1)
    RES <- matrix(DN %*% res, nrow = N, byrow = T)

    # singularity problems
    if (comb == "sam" & N < kt) {
      stop("N < (k* + m): it could lead to singularity problems if comb == sam", call. = FALSE)
    }

    if (comb == "acov" & N < m) {
      stop("N < m: it could lead to singularity problems if comb == acov", call. = FALSE)
    }
  }

  if (mse) {
    cov_mod <- function(x, ...) crossprod(stats::na.omit(x), ...) / NROW(stats::na.omit(x))
  } else {
    cov_mod <- function(x, ...) stats::var(x, na.rm = TRUE, ...)
  }

  if (corpcor) {
    shr_mod <- function(x, ...) Matrix(unclass(corpcor::cov.shrink(x, verbose = FALSE, lambda.var=0, ...)))
  } else {
    shr_mod <- function(x, ...) shrink_estim(x, minT = mse)[[1]]
  }

  # Reconciliation

  switch(comb,
         bu = {
           if (NCOL(BASEF) != m) {
             BASEF <- BASEF[, (ks + 1):kt]
           }

           if(nn){
             BASEF <- BASEF * (BASEF > 0)
           }
           OUTF <- BASEF %*% t(R)

           outf <- as.vector(t(Dh) %*% as.vector(t(OUTF)))

           outf <- stats::setNames(outf, paste("k", rep(kset, h * (m/kset)), "h",
                                               do.call("c", as.list(sapply(
                                                 (m/kset) * h,
                                                 function(x) seq(1:x)
                                               ))),
                                               sep = ""
           ))
           if (keep == "list") {
             return(list(
               recf = outf, R = R,
               M = R %*% cbind(matrix(0, m, ks), diag(m))
             ))
           } else {
             return(outf)
           }
         },
         ols =
           Omega <- .sparseDiagonal(kt),
         struc =
           Omega <- .sparseDiagonal(x = rowSums(R)),
         wlsv = {
           var_freq <- sapply(kset, function(x) cov_mod(res[rep(kset, (m/kset) * N) == x]))
           Omega <- .sparseDiagonal(x = rep(var_freq, (m/kset)))
         },
         wlsh = {
           diagO <- diag(cov_mod(RES))
           Omega <- .sparseDiagonal(x = diagO)
         },
         acov = {
           Omega <- Matrix::bdiag(lapply(kset, function(x) cov_mod(RES[, rep(kset, (m/kset)) == x])))
         },
         strar1 = {
           rho <- lapply(kset, function(x) stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]),
                                                      1, plot = F)$acf[2, 1, 1])
           expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

           Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
           Ostr2 <- .sparseDiagonal(x = apply(R, 1, sum))^0.5
           Omega <- Ostr2 %*% Gam %*% Ostr2
         },
         sar1 = {
           rho <- lapply(kset, function(x) stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]),
                                                      1, plot = F)$acf[2, 1, 1])
           expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

           Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
           var_freq <- sapply(kset, function(x) cov_mod(res[rep(kset, (m/kset) * N) == x]))
           Os2 <- .sparseDiagonal(x = rep(var_freq, (m/kset)))^0.5
           Omega <- Os2 %*% Gam %*% Os2
         },
         har1 = {
           rho <- lapply(kset, function(x) stats::acf(stats::na.omit(res[rep(kset, (m/kset) * N) == x]),
                                                      1, plot = F)$acf[2, 1, 1])
           expo <- lapply((m/kset), function(x) toeplitz(1:x) - 1)

           Gam <- Matrix::bdiag(Map(function(x, y) x^y, x = rho, y = expo))
           diagO <- diag(cov_mod(RES))
           Oh2 <- .sparseDiagonal(x = diagO)^0.5
           Omega <- Matrix::Matrix(Oh2 %*% Gam %*% Oh2)
         },
         shr = {
           Omega <- shr_mod(RES)
         },
         sam = {
           Omega <- cov_mod(RES)
         },
         omega = {
           if (is.null(Omega)) {
             stop("Please, put in option Omega your covariance matrix", call. = FALSE)
           }
           Omega <- Omega
         }
  )

  b_pos <- c(rep(0, kt - m), rep(1, m))

  if(!is.null(v)){
    keep <- "recf"
    rec_sol <- recoV(
      basef = BASEF, W = Omega, Ht = Zt, sol = sol, nn = nn, keep = keep, S = R, type = type,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type, v = v
    )
  }else if(type == "S"){
    rec_sol <- recoS(
      basef = BASEF, W = Omega, S = R, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
    )
  }else{
    rec_sol <- recoM(
      basef = BASEF, W = Omega, Ht = Zt, sol = sol, nn = nn, keep = keep, S = R,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
    )
  }

  if (keep == "list") {
    rec_sol$nn_check <- sum(rec_sol$recf < 0)
    rec_sol$rec_check <- all(rec_sol$recf %*% t(Zt) < 1e-6)

    rec_sol$recf <- as.vector(t(Dh) %*% as.vector(t(rec_sol$recf)))
    names(rec_sol)[names(rec_sol) == "W"] <- "Omega"
    rec_sol$Omega <- as.matrix(rec_sol$Omega)
    dimnames(rec_sol$Omega) <- NULL

    names_all_list <- c("recf", "W", "Omega", "nn_check", "rec_check", "varf", "M", "G", "S", "info")
    names_list <- names(rec_sol)
    rec_sol <- rec_sol[names_list[order(match(names_list, names_all_list))]]

    rec_sol$recf <- stats::setNames(rec_sol$recf, paste("k", rep(kset, h * (m/kset)), "h",
                                                        do.call("c", as.list(sapply(
                                                          (m/kset) * h,
                                                          function(x) seq(1:x)
                                                        ))),
                                                        sep = ""
    ))
    return(rec_sol)
  } else {
    if (length(rec_sol) == 1) {
      rec_sol$recf <- as.vector(t(Dh) %*% as.vector(t(rec_sol$recf)))
      rec_sol$recf <- stats::setNames(rec_sol$recf, paste("k", rep(kset, h * (m/kset)), "h",
                                                          do.call("c", as.list(sapply(
                                                            (m/kset) * h,
                                                            function(x) seq(1:x)
                                                          ))),
                                                          sep = ""
      ))
      return(rec_sol$recf)
    } else {
      rec_sol$recf <- as.vector(t(Dh) %*% as.vector(t(rec_sol$recf)))
      rec_sol$recf <- stats::setNames(rec_sol$recf, paste("k", rep(kset, h * (m/kset)), "h",
                                                          do.call("c", as.list(sapply(
                                                            (m/kset) * h,
                                                            function(x) seq(1:x)
                                                          ))),
                                                          sep = ""
      ))
      return(rec_sol)
    }
  }
}

#' @export
thfrec.list <- function(basef, m, ..., Omega){
  if(missing(m)){
    stop("The argument m is not specified", call. = FALSE)
  }

  # Prepare basef
  tools <- thf_tools(m=m)
  m <- max(tools$kset)
  h <- length(basef) / tools$kt
  if(h == 1){
    stop("You don't need thfrech for h = 1, please use thfrec()", call. = FALSE)
  }

  Dh <- Dmat(h = h, m = tools$kset, n = 1)
  baseh <- matrix(Dh%*%basef, nrow = h, byrow = T)

  if(h != length(Omega) | !is.list(Omega)){
    stop("Omega must be a list with ", h, " matrices", call. = FALSE)
  }
  # Reconciliation
  obj <- Map(function(b, Om){
    thfrec(basef = b,         # vector of base forecasts
           comb = "omega",    # custom combination
           Omega = Om,        # custom Omega
           m = tools$kset,
           ...)
  },
  # list of base forecasts divide by forecasts horizons
  b = split(baseh, 1:NROW(baseh)),
  # list of Omega divide by forecasts horizons
  Om = Omega)

  out <- list()

  out$recf <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "recf")))
  out$recf <- as.vector(t(Dh) %*% as.vector(t(out$recf)))
  out$recf <- stats::setNames(out$recf, paste("k", rep(tools$kset, h * (m/tools$kset)), "h",
                                              do.call("c", as.list(sapply(
                                                (m/tools$kset) * h,
                                                function(x) seq(1:x)
                                              ))),
                                              sep = ""
  ))

  out$varf <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "varf")))
  if(all(is.na(out$varf))){
    out$varf <- NULL
  }else{
    out$varf <- as.vector(t(Dh) %*% as.vector(t(out$varf)))
    out$varf <- stats::setNames(out$varf, paste("k", rep(tools$kset, h * (m/tools$kset)), "h",
                                                do.call("c", as.list(sapply(
                                                  (m/tools$kset) * h,
                                                  function(x) seq(1:x)
                                                ))),
                                                sep = ""
    ))
    out$varf <- out$varf[!is.na(out$varf)]
  }

  out$info <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "info")))
  if(all(is.na(out$info))){
    out$info <- NULL
  }else{
    rownames(out$info) <-  paste("h",1:NROW(baseh), sep="")
  }

  return(out)
}
