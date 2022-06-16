#' @title Optimal combination cross-temporal forecast reconciliation
#'
#' @description
#' \loadmathjax
#' Optimal (in least squares sense) combination cross-temporal forecast
#' reconciliation. The reconciled forecasts are calculated either through a
#' projection approach (Byron, 1978), or the equivalent structural approach
#' by Hyndman et al. (2011).
#'
#' @usage octrec(basef, m, C, comb, res, Ut, nb, mse = TRUE, corpcor = FALSE,
#'        type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
#'        nn_type = "osqp", settings = list(), bounds = NULL, W = NULL,
#'        Omega = NULL)
#'
#' @param basef  (\mjseqn{n \times h(k^\ast+m)}) matrix of base forecasts to be
#' reconciled, \mjseqn{\widehat{\mathbf{Y}}}; \mjseqn{n} is the total number of variables,
#' \mjseqn{m} is the highest time frequency, \mjseqn{k^\ast} is the sum of (a
#' subset of) (\mjseqn{p-1}) factors of \mjseqn{m}, excluding \mjseqn{m}, and
#' \mjseqn{h} is the forecast horizon for the lowest frequency time series.
#' Each row identifies a time series, and the forecasts are ordered as
#' [lowest_freq' ...  highest_freq']'.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#' @param comb Type of the reconciliation. It corresponds to a specific
#' (\mjseqn{n(k\ast + m) \times n(k^\ast + m)}) covariance matrix, where
#' \mjseqn{k^\ast} is the sum of (a subset of) (\mjseqn{p-1}) factors of
#' \mjseqn{m} (\mjseqn{m} is not considered) and \mjseqn{n} is the number
#' of variables:
#' \itemize{
#'   \item \bold{ols} (Identity);
#'   \item \bold{struc} (Cross-temporal structural variances);
#'   \item \bold{cs_struc} (Cross-sectional structural variances and temporal independence);
#'   \item \bold{t_struc} (Cross-sectional independence and temporal structural variances);
#'   \item \bold{wlsh} (Hierarchy variances);
#'   \item \bold{wlsv} (Series variances);
#'   \item \bold{bdshr} (Shrunk cross-covariance matrix, cross-sectional framework);
#'   \item \bold{bdsam} (Sample cross-covariance matrix, cross-sectional framework);
#'   \item \bold{acov} (Series auto-covariance matrix);
#'   \item \bold{Sshr} (Series shrunk cross-covariance matrix);
#'   \item \bold{Ssam} (Series cross-covariance matrix);
#'   \item \bold{shr}  (Shrunk cross-covariance matrix);
#'   \item \bold{sam} (Sample cross-covariance matrix);
#'   \item \bold{w} use your personal matrix W in param \code{W};
#'   \item \bold{omega} use your personal matrix Omega in param \code{Omega}.
#' }
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \mjseqn{(\textbf{U}'\textbf{y} = \mathbf{0})} spanning the null space valid
#' for the reconciled forecasts. It can be used instead of parameter
#' \code{C}, but \code{nb} (\mjseqn{n = n_a + n_b}) is needed if
#' \mjseqn{\textbf{U}' \neq [\textbf{I} \ -\textbf{C}]}. If the hierarchy
#' admits a structural representation, \mjseqn{\textbf{U}'} has dimension
#' (\mjseqn{n_a \times n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb}
#' and \code{Ut} are not used.
#' @param res (\mjseqn{n \times N(k^\ast + m)}) matrix containing the residuals at
#' all the temporal frequencies ordered [lowest_freq' ...  highest_freq']'
#' (columns) for each variable (row), needed to estimate the covariance matrix
#' when \code{comb =} \code{\{"sam",} \code{"wlsv",} \code{"wlsh",}
#' \code{"acov",} \code{"Ssam",} \code{"Sshr",} \code{"Sshr1",} \code{"shr"\}}.
#' @param W,Omega This option permits to directly enter the covariance matrix:
#' \enumerate{
#'   \item \code{W} must be a p.d. (\mjseqn{n(k^\ast + m) \times n(k^\ast + m)})
#'   matrix or a list of \mjseqn{h} matrix (one for each forecast horizon);
#'   \item \code{Omega} must be a p.d. (\mjseqn{n(k^\ast + m) \times n(k^\ast + m)})
#'   matrix or a list of \code{h} matrix (one for each forecast horizon);
#'   \item if \code{comb} is different from "\code{w}" or "\code{omega}",
#'   \code{W} or \code{Omega} is not used.
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
#' structural approach with summing matrix S.
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
#' @param bounds (\mjseqn{n(k^\ast + m) \times 2}) matrix of the bounds on the
#' variables: the first column is the lower bound, and the second column is the
#' upper bound.
#' @param v vector index of the fixed base forecast (\mjseqn{\mbox{min}(v) > 0}
#' and \mjseqn{\mbox{max}(v) < n(k^\ast + m)}).
#'
#' @details
#' Considering contemporaneous and temporal dimensions in the
#' same framework requires to extend and adapt the notations
#' used in \link[FoReco]{htsrec} and \link[FoReco]{thfrec}.
#' To do that, we define the matrix containing the base forecasts
#' at any considered temporal frequency as
#' \mjsdeqn{
#' \widehat{\textbf{Y}}_{n \times h(k^\ast+m))} =
#' \left[\begin{array}{ccccc}
#' \widehat{\textbf{A}}^{[m]} & \widehat{\textbf{A}}^{[k_{p-1}]} & \cdots &
#' \widehat{\textbf{A}}^{[k_2]} & \widehat{\textbf{A}}^{[1]} \cr
#' \widehat{\textbf{B}}^{[m]} & \widehat{\textbf{B}}^{[k_{p-1}]} & \cdots &
#' \widehat{\textbf{B}}^{[k_2]} & \widehat{\textbf{B}}^{[1]}
#' \end{array}
#' \right] \qquad k \in {\cal K},}
#' where \mjseqn{\cal K} is a subset of \mjseqn{p} factors of \mjseqn{m} and,
#' \mjseqn{\widehat{\textbf{B}}^{[k]}} and \mjseqn{\widehat{\textbf{A}}^{[k]}}
#' are the matrices containing the \mjseqn{k}-order temporal aggregates of the
#' bts and uts, of dimension (\mjseqn{n_b \times h m/k}) and
#' (\mjseqn{n_a \times h m/k}), respectively.
#'
#' Let us consider the multivariate regression model
#' \mjsdeqn{\widehat{\mathbf{Y}} = \mathbf{Y} + \mathbf{E} ,}
#' where the involved matrices have each dimension
#' \mjseqn{[n \times (k^\ast+m)]} and contain, respectively, the base
#' (\mjseqn{\widehat{\mathbf{Y}}}) and the target forecasts
#' (\mjseqn{\mathbf{Y}}), and the coherency errors (\mjseqn{\mathbf{E}}) for
#' the \mjseqn{n} component variables of the linearly constrained time series
#' of interest. For each variable, \mjseqn{k^\ast + m} base forecasts are
#' available, pertaining to all aggregation levels of the temporal hierarchy
#' for a complete cycle of high-frequency observation, \mjseqn{m}. Consider
#' now two vectorized versions of model, by transforming the matrices either
#' in original form:
#' \mjsdeqn{\mbox{vec}\left(\widehat{\mathbf{Y}}\right) =
#' \mbox{vec}\left(\mathbf{Y}\right) + \mathbf{\varepsilon} \; \mbox{ with }
#' \; \mathbf{\varepsilon} = \mbox{vec}\left(\mathbf{E}\right)}
#' or in transposed form:
#' \mjsdeqn{\mbox{vec}\left(\widehat{\mathbf{Y}}'\right) =
#' \mbox{vec}\left(\mathbf{Y}'\right) + \mathbf{\eta} \; \mbox{ with }
#' \; \mathbf{\eta} = \mbox{vec}\left(\mathbf{E}'\right).}
#' Denote with \mjseqn{\mathbf{P}} the \mjseqn{[n(k^\ast+m) \times n(k^\ast+m)]}
#' commutation matrix such that
#' \mjseqn{\mathbf{P}\mbox{vec}(\mathbf{Y}) = \mbox{vec}(\mathbf{Y}')},
#' \mjseqn{\mathbf{P}\mbox{vec}(\widehat{\mathbf{Y}}) = \mbox{vec}(\widehat{\mathbf{Y}}')}
#' and \mjseqn{\mathbf{P}\mathbf{\varepsilon} = {\bf \eta}}.
#' Let \mjseqn{\mathbf{W} = \mathrm{E}[\mathbf{\varepsilon\varepsilon}']} be the
#' covariance matrix of vector \mjseqn{\mathbf{\varepsilon}}, and
#' \mjseqn{\mathbf{\Omega} = \mathrm{E}[\mathbf{\eta\eta}']} the covariance matrix of
#' vector \mjseqn{\mathbf{\eta}}. Clearly, \mjseqn{\mathbf{W}} and
#' \mjseqn{\mathbf{\Omega}} are different parameterizations of the same
#' statistical object for which the following relationships hold:
#' \mjsdeqn{\mathbf{\Omega} = \mathbf{P}\mathbf{W}\mathbf{P}',
#' \qquad \mathbf{W} = \mathbf{P}' \mathbf{\Omega}\mathbf{P} .}
#' In order to apply the general point forecast reconciliation according to the
#' projection approach (\code{type = "M"}) to a cross-temporal forecast
#' reconciliation problem, we may consider either two \emph{vec}-forms , e.g.
#' if we follow the first:
#' \mjsdeqn{
#' \tilde{\mathbf{y}}= \hat{\mathbf{y}} - \mathbf{\Omega}\mathbf{H}\left(
#' \mathbf{H}'\mathbf{\Omega}\mathbf{H}\right)^{-1}\mathbf{H}'\hat{\mathbf{y}} =
#' {\mathbf{M}}\hat{\mathbf{y}},}
#' where \mjseqn{\widehat{\mathbf{y}} = \mbox{vec}(\widehat{\mathbf{Y}}')} is the
#' row vectorization of the base forecasts matrix \mjseqn{\widehat{\mathbf{Y}}}
#' The alternative equivalent solution (\code{type = "S"}) (following the
#' structural reconciliation approach by Hyndman et al., 2011) is
#' \mjsdeqn{\widetilde{\mathbf{y}} = \widetilde{\mathbf{S}}\left(\widetilde{\mathbf{S}}'
#' \mathbf{\Omega}^{-1}\widetilde{\mathbf{S}}\right)^{-1}\widetilde{\mathbf{S}}'
#' \mathbf{\Omega}^{-1}\widehat{\mathbf{y}} = \widetilde{\mathbf{S}}\mathbf{G}\widehat{\mathbf{y}}.}
#' where \mjseqn{\widetilde{\mathbf{S}}} is the cross-temporal summing matrix.
#'
#' \strong{Bounds on the reconciled forecasts}
#'
#' When the reconciliation uses the optimization package osqp,
#' the user may impose bounds on the reconciled forecasts.
#' The parameter \code{bounds} permits to consider lower (\mjseqn{\mathbf{a}}) and
#' upper (\mjseqn{\mathbf{b}}) bounds like \mjseqn{\mathbf{a} \leq
#' \widetilde{\mathbf{y}} \leq \mathbf{b}}, where \mjseqn{\widehat{\mathbf{y}} =
#' \mbox{vec}(\widehat{\mathbf{Y}}')}, such that:
#' \mjsdeqn{ \begin{array}{c}
#' a_1 \leq \widetilde{y}_1 \leq b_1 \cr
#' \dots \cr
#' a_{n(k^\ast + m)} \leq \widetilde{y}_{n(k^\ast + m)} \leq b_{n(k^\ast + m)} \cr
#' \end{array} \Rightarrow
#' \mbox{bounds} = [\mathbf{a} \; \mathbf{b}] =
#' \left[\begin{array}{cc}
#' a_1 & b_1 \cr
#' \vdots & \vdots \cr
#' a_{n(k^\ast + m)} & b_{n(k^\ast + m)} \cr
#' \end{array}\right],}
#' where \mjseqn{a_i \in [- \infty, + \infty]} and \mjseqn{b_i \in [- \infty, + \infty]}.
#' If \mjseqn{y_i} is unbounded, the i-th row of \code{bounds} would be equal
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
#' returns only the (\mjseqn{n \times h(k^\ast + m)}) reconciled forecasts
#' matrix, otherwise (\code{keep="all"}) it returns a list that mainly depends
#' on what type of representation (\code{type}) and solution technique
#' (\code{sol}) have been used:
#' \item{\code{recf}}{(\mjseqn{n \times h(k^\ast + m)}) reconciled forecasts matrix, \mjseqn{\widetilde{\textbf{Y}}}.}
#' \item{\code{Omega}}{Covariance matrix used for reconciled forecasts (\mjseqn{\mbox{vec}(\widehat{\textbf{Y}}')} representation).}
#' \item{\code{W}}{Covariance matrix used for reconciled forecasts (\mjseqn{\mbox{vec}(\widehat{\textbf{Y}})} representation).}
#' \item{\code{nn_check}}{Number of negative values (if zero, there are no values below zero).}
#' \item{\code{rec_check}}{Logical value: \code{rec_check = TRUE} when the constraints have been fulfilled,}
#' \item{\code{varf} (\code{type="direct"})}{(\mjseqn{n \times (k^\ast + m)}) reconciled forecasts variance matrix for \mjseqn{h=1}, \mjseqn{\mbox{diag}(\mathbf{MW}}).}
#' \item{\code{M} (\code{type="direct"})}{Projection matrix (projection approach).}
#' \item{\code{G} (\code{type="S"} and \code{type="direct"})}{Projection matrix (structural approach, \mjseqn{\mathbf{M}=\mathbf{S}\mathbf{G}}).}
#' \item{\code{S} (\code{type="S"} and \code{type="direct"})}{Cross-temporal summing matrix (\mjseqn{\widetilde{\textbf{F}}\mbox{vec}(\widehat{\textbf{Y}}')} representation).}
#' \item{\code{info} (\code{type="osqp"})}{matrix with some useful indicators (columns)
#' for each forecast horizon \mjseqn{h} (rows): run time (\code{run_time}), number of iteration,
#' norm of primal residual (\code{pri_res}), status of osqp's solution (\code{status}) and
#' polish's status (\code{status_polish}).}
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social accounts matrices,
#' \emph{Journal of the Royal Statistical Society A}, 141, 3, 359-367.
#'
#' Di Fonzo, T., and Girolimetto, D. (2021), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, in press.
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
#' Quadratic Programming Solver using the `OSQP' Library, R package version 0.6.0.3
#' (October 10, 2019), \href{https://CRAN.R-project.org/package=osqp}{https://CRAN.R-project.org/package=osqp}.
#'
#' @family reconciliation procedures
#' @examples
#' data(FoReco_data)
#' obj <- octrec(FoReco_data$base, m = 12, C = FoReco_data$C,
#'               comb = "bdshr", res = FoReco_data$res)
#'
#' @export
#'
#' @import Matrix osqp methods
#'
octrec <- function(basef, m, C, comb, res, Ut, nb, mse = TRUE, corpcor = FALSE,
                   type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
                   nn_type = "osqp", settings = list(), bounds = NULL, W = NULL, Omega = NULL) {

  if(missing(comb)){
    stop("The argument comb is not specified.", call. = FALSE)
  }else if(comb != "omega" | comb != "w"){
    Omega <-  NULL
    W <-  NULL
  }

  if(is.list(W)){
    UseMethod("octrec", W)
  }else{
    UseMethod("octrec", Omega)
  }
}

#' @export
octrec.default <- function(basef, m, C, comb, res, Ut, nb, mse = TRUE,
                           corpcor = FALSE, type = "M", sol = "direct",
                           keep = "list", v = NULL, nn = FALSE,
                           nn_type = "osqp", settings = osqpSettings(),
                           bounds = NULL, W = NULL, Omega = NULL) {
  if (missing(m)) {
    stop("The argument m is not specified", call. = FALSE)
  }

  comb <- match.arg(comb, c(
    "ols", "struc", "sam", "wlsv", "wlsh", "shr",
    "acov", "Ssam", "Sshr", "bdshr", "bdsam", "w", "omega",
    "cs_struc", "t_struc"))

  type <- match.arg(type, c("M", "S"))
  keep <- match.arg(keep, c("list", "recf"))

  nn_type <- match.arg(nn_type, c("osqp", "KAnn", "fbpp", "sntz"))

  if (missing(basef)) {
    stop("The argument basef is not specified", call. = FALSE)
  }

  if (missing(C)) {
    if (missing(Ut)) {
      stop("Please, give C or Ut.", call. = FALSE)
    } else if(missing(nb)){
      ctf <- ctf_tools(Ut = Ut, m = m, h = 1, sparse = TRUE)
    } else {
      ctf <- ctf_tools(Ut = Ut, nb = nb, m = m, h = 1, sparse = TRUE)
    }
  }else{
    ctf <- ctf_tools(C = C, m = m,  h = 1, sparse = TRUE)
  }

  n <- ctf$hts$n
  na <- ctf$hts$na
  nb <- ctf$hts$nb
  #C <- ctf$hts$C
  Scs <- ctf$hts$S
  #Ut <- ctf$hts$Ut

  kset <- ctf$thf$kset
  m <- ctf$thf$m
  p <- ctf$thf$p
  kt <- ctf$thf$kt
  ks <- ctf$thf$ks
  R <- ctf$thf$R

  # matrix
  #Zt <- tools$Zt

  if (NROW(basef) != n) {
    stop("Incorrect dimension of Ut or basef (they don't have same columns).", call. = FALSE)
  }

  Ht <- ctf$ctf$Ht
  S <- ctf$ctf$Fmat

  # Base forecast
  if (NCOL(basef) %% kt != 0) {
    stop("basef has a number of row not in line with the frequency of the series", call. = FALSE)
  }

  if(nn){
    if(nn_type == "fbpp" | nn_type == "KAnn"){
      type = "M"
    }
  }

  if(is.null(S)){
    if (type == "S") {
      stop("Type = S needs the cross-sectional C matrix. \nPlease consider using gecoma() to get the right C when working with Ut.", call. = FALSE)
    }

    if(nn_type == "sntz"){
      stop("nn_type = sntz needs the cross-sectional C matrix.\nPlease consider using gecoma() to get the right C when working with Ut.", call. = FALSE)
    }

    if (comb == "struc" | comb == "cs_struc") {
      stop("comb = ", comb, " needs the cross-sectional C matrix.\nPlease consider using gecoma() to get the right C when working with Ut.", call. = FALSE)
    }
  }

  h <- NCOL(basef) / kt
  Dh <- Dmat(h = h, m = kset, n = n)
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

    DN <- Dmat(h = N, m = kset, n = n)
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
    shr_mod <- function(x, ...) Matrix(unclass(corpcor::cov.shrink(x, verbose = FALSE, lambda.var=0, ...)))
  } else {
    shr_mod <- function(x, ...) shrink_estim(x, minT = mse)[[1]]
  }

  switch(comb,
         ols = {
           Omega <- .sparseDiagonal(n * kt)
         },
         struc = {
           Omega <- .sparseDiagonal(x = rowSums(S))
         },
         "cs_struc" = {
           Omega <- .sparseDiagonal(x = rep(rowSums(Scs), each = kt))
         },
         "t_struc" = {
           Omega <- .sparseDiagonal(x = rep(rowSums(R), n))
         },
         sam = {
           Omega <- cov_mod(E)
         },
         wlsh = {
           Omega <- .sparseDiagonal(x = diag(cov_mod(E)))
         },
         wlsv = {
           var_freq <- apply(res, 1, function(z) sapply(kset, function(x) cov_mod(z[rep(kset, (m/kset) * N) == x])))
           Omega <- .sparseDiagonal(x = rep(as.vector(var_freq), rep((m/kset), n)))
         },
         shr = {
           Omega <- shr_mod(E)
         },
         acov = {
           mat1 <- bdiag(rep(lapply((m/kset), function(x) matrix(1, nrow = x, ncol = x)), n))
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
           if (is.null(W)) {
             stop("Please, put in option W your covariance matrix", call. = FALSE)
           }
           P <- commat(n, kt)
           Omega <- P %*% W %*% t(P)
         },
         omega = {
           if (is.null(Omega)) {
             stop("Please, put in option Omega your covariance matrix", call. = FALSE)
           }
           Omega <- Omega
         },
         bdshr = {
           blockW <- lapply(kset, function(x) shr_mod(t(res[, rep(kset, N * (m/kset)) == x])))
           blockW <- rep(blockW, (m/kset))
           P <- commat(n, kt)
           Omega <- P %*% bdiag(blockW) %*% t(P)
         },
         bdsam = {
           blockW <- lapply(kset, function(x) cov_mod(t(res[, rep(kset, N * (m/kset)) == x])))
           blockW <- rep(blockW, (m/kset))
           P <- commat(n, kt)
           Omega <- P %*% bdiag(blockW) %*% t(P)
         }
  )

  b_pos <- c(rep(0, na * kt), rep(rep(kset, (m/kset)), nb) == 1)

  if(!is.null(v)){
    keep <- "recf"
    rec_sol <- recoV(
      basef = Ybase, W = Omega, Ht = Ht, sol = sol, nn = nn, keep = keep, S = S, type = type,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type, v = v
    )
  }else if (type == "S") {
    rec_sol <- recoS(
      basef = Ybase, W = Omega, S = S, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
    )
  } else {
    rec_sol <- recoM(
      basef = Ybase, W = Omega, Ht = Ht, sol = sol, nn = nn, keep = keep, S = S,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
    )
  }

  if (keep == "list") {
    P <- commat(n, kt)
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
    colnames(rec_sol$recf) <- paste("k", rep(kset, h * (m/kset)), "h",
                                    do.call("c", as.list(sapply(
                                      (m/kset) * h,
                                      function(x) seq(1:x)
                                    ))),
                                    sep = ""
    )
    return(rec_sol)
  } else {
    if (length(rec_sol) == 1) {
      rec_sol$recf <- matrix(t(Dh) %*% as.vector(t(rec_sol$recf)), nrow = n, byrow = TRUE)
      rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
      colnames(rec_sol$recf) <- paste("k", rep(kset, h * (m/kset)), "h",
                                      do.call("c", as.list(sapply(
                                        (m/kset) * h,
                                        function(x) seq(1:x)
                                      ))),
                                      sep = ""
      )
      return(rec_sol$recf)
    } else {
      rec_sol$recf <- matrix(t(Dh) %*% as.vector(t(rec_sol$recf)), nrow = n, byrow = TRUE)
      rownames(rec_sol$recf) <- if (is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
      colnames(rec_sol$recf) <- paste("k", rep(kset, h * (m/kset)), "h",
                                      do.call("c", as.list(sapply(
                                        (m/kset) * h,
                                        function(x) seq(1:x)
                                      ))),
                                      sep = ""
      )
      return(rec_sol)
    }
  }
}

#' @export
octrec.list <- function(basef, m, ..., W, Omega){
  if(missing(m)){
    stop("The argument m is not specified", call. = FALSE)
  }

  if(!is.matrix(basef)){
    stop("basef must be a matrix", call. = FALSE)
  }

  # Prepare basef
  tools <- thf_tools(m=m)
  m <- max(tools$kset)
  h <- NCOL(basef) / tools$kt
  if(h == 1){
    stop("You don't need octrech for h = 1, please use octrec()", call. = FALSE)
  }

  n <- NROW(basef)

  Dh <- Dmat(h = h, m = tools$kset, n = n)
  Ybase <- matrix(Dh %*% as.vector(t(basef)), nrow = h, byrow = TRUE)
  baseh <- lapply(split(Ybase,1:h), matrix, nrow = 3, byrow = TRUE)

  # Reconciliation
  if(is.list(W)){
    obj <- Map(function(b, W){
      octrec.default(basef = b,         # vector of base forecasts
                     comb = "w",    # custom combination
                     W = W,
                     m = tools$kset,
                     ...)
    },
    # list of base forecasts divide by forecasts horizons
    b = baseh,
    # list of Omega divide by forecasts horizons
    W = W)
  }else{
    obj <- Map(function(b, Om){
      octrec.default(basef = b,         # vector of base forecasts
                     comb = "w",    # custom combination
                     Omega = Om,
                     m = tools$kset,
                     ...)
    },
    # list of base forecasts divide by forecasts horizons
    b = baseh,
    # list of Omega divide by forecasts horizons
    Om = Omega)
  }

  out <- list()
  recf <- lapply(obj, function(x) extract_data(x = x, name = "recf"))
  recf <- do.call(rbind, lapply(recf, function(x) as.vector(t(x))))
  out$recf <- matrix(t(Dh) %*% as.vector(t(recf)), nrow = n, byrow = TRUE)
  colnames(out$recf) <- paste("k", rep(tools$kset, h * (m/tools$kset)), "h",
                              do.call("c", as.list(sapply(
                                (m/tools$kset) * h,
                                function(x) seq(1:x)
                              ))),
                              sep = ""
  )
  rownames(out$recf) <- if(is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)

  varf <- lapply(obj, function(x) extract_data(x = x, name = "varf"))
  if(all(is.na(out$varf))){
    out$varf <- NULL
  }else{
    varf <- do.call(rbind, lapply(varf, function(x) as.vector(t(x))))
    out$varf <- matrix(t(Dh) %*% as.vector(t(varf)), nrow = n, byrow = TRUE)
    colnames(out$varf) <- paste("k", rep(tools$kset, h * (m/tools$kset)), "h",
                                do.call("c", as.list(sapply(
                                  (m/tools$kset) * h,
                                  function(x) seq(1:x)
                                ))),
                                sep = "")
    rownames(out$varf) <- if(is.null(rownames(basef))) paste("serie", 1:n, sep = "") else rownames(basef)
    out$varf <- out$varf[ , apply(out$varf, 2, function(x) !any(is.na(x)))]
  }

  out$info <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "info")))
  if(all(is.na(out$info))){
    out$info <- NULL
  }else{
    rownames(out$info) <-  paste("h",1:NROW(baseh), sep="")
    out$info <- out$info[apply(out$info, 1, function(x) !any(is.na(x))), , drop = FALSE]
  }

  return(out)
}
