#' @title Cross-sectional (contemporaneous) forecast reconciliation
#'
#' @description
#' Cross-sectional (contemporaneous) forecast reconciliation of a linearly constrained
#' (e.g., hierarchical/grouped) multiple time series.
#' The reconciled forecasts are calculated either through a
#' projection approach (Byron, 1978, see also van Erven and Cugliari, 2015, and
#' Wickramasuriya et al., 2019), or the equivalent structural approach by Hyndman
#' et al. (2011). Moreover, the classic bottom-up approach is available.
#'
#' @usage htsrec(basef, comb, C, res, Ut, nb, mse = TRUE, corpcor = FALSE,
#'        type = "M", sol = "direct", keep = "list",  v = NULL, nn = FALSE,
#'        nn_type = "osqp", settings = list(), bounds = NULL, W = NULL)
#'
#' @param basef (\mjseqn{h \times n}) matrix of base forecasts to be reconciled;
#' \mjseqn{h} is the forecast horizon and \mjseqn{n} is the total number of time series.
#' @param comb Type of the reconciliation. Except for Bottom-up, each
#' option corresponds to a specific (\mjseqn{n \times n}) covariance matrix:
#' \itemize{
#'     \item \bold{bu} (Bottom-up);
#'     \item \bold{ols} (Identity);
#'     \item \bold{struc} (Structural variances);
#'     \item \bold{wls} (Series variances) - uses res;
#'     \item \bold{shr} (Shrunk covariance matrix - MinT-shr) - uses res;
#'     \item \bold{sam} (Sample covariance matrix - MinT-sam) - uses res;
#'     \item \bold{w} use your personal matrix W in param \code{W}.
#'   }
#' @param C (\mjseqn{n_a \times n_b}) cross-sectional (contemporaneous) matrix
#' mapping the bottom level series into the higher level ones.
#' @param Ut Zero constraints cross-sectional (contemporaneous) kernel matrix
#' \mjseqn{(\mathbf{U}'\mathbf{y} = \mathbf{0})} spanning the null space valid
#' for the reconciled forecasts. It can be used instead of parameter
#' \code{C}, but \code{nb} (\mjseqn{n = n_a + n_b}) is needed if
#' \mjseqn{\mathbf{U}' \neq [\mathbf{I} \ -\mathbf{C}]}{}. If the hierarchy
#' admits a structural representation, \mjseqn{\mathbf{U}'} has dimension
#' (\mjseqn{n_a \times n}).
#' @param nb Number of bottom time series; if \code{C} is present, \code{nb}
#' and \code{Ut} are not used.
#' @param res (\mjseqn{N \times n}) in-sample residuals matrix needed when \code{comb =}
#' \code{\{"wls",} \code{"shr",} \code{"sam"\}}. The columns must be in
#' the same order as \code{basef}.
#' @param W This option permits to directly enter the covariance matrix:
#'   \enumerate{
#'     \item \code{W} must be a p.d. (\mjseqn{n \times n}) matrix or a list of \mjseqn{h} matrix (one for each forecast horizon);
#'     \item if \code{comb} is different from "\code{w}", \code{W} is not used.
#'   }
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
#' @param bounds (\mjseqn{n \times 2}) matrix containing the cross-sectional bounds:
#' the first column is the lower bound, and the second column is the upper bound.
#' @param v vector index of the fixed base forecast (\mjseqn{\mbox{min}(v) > 0}
#' and \mjseqn{\mbox{max}(v) < n}).
#'
#' @details
#' \loadmathjax
#' Let \mjseqn{\mathbf{y}} be a (\mjseqn{n \times 1}) vector of target point
#' forecasts which are wished to satisfy the system of linearly independent
#' constraints \mjsdeqn{\mathbf{U}'\mathbf{y} = \mathbf{0}_{(r \times 1)},}
#' where \mjseqn{\mathbf{U}'} is a (\mjseqn{r \times n}) matrix, with
#' rank\mjseqn{(\mathbf{U}') = r \leq n}, and \mjseqn{\mathbf{0}_{(r \times 1)}}
#' is a (\mjseqn{r \times 1}) null vector. Let \mjseqn{\widehat{\mathbf{y}}}
#' be a (\mjseqn{n \times 1}) vector of unbiased point forecasts, not
#' fulfilling the linear constraints (i.e., \mjseqn{\mathbf{U}'\widehat{\mathbf{y}}
#' \ne \mathbf{0}}).
#'
#' We consider a regression-based reconciliation method assuming
#' that \mjseqn{\widehat{\mathbf{y}}} is related to \mjseqn{\mathbf{y}} by
#' \mjsdeqn{\widehat{\mathbf{y}} = \mathbf{y} + \mathbf{\varepsilon},}
#' where \mjseqn{\mathbf{\varepsilon}} is a (\mjseqn{n \times 1}) vector
#' of zero mean disturbances, with known p.d. covariance matrix
#' \mjseqn{\mathbf{W}}. The reconciled forecasts \mjseqn{\widetilde{\mathbf{y}}}
#' are found by minimizing the generalized least squares (GLS) objective
#' function \mjseqn{\left(\widehat{\mathbf{y}} - \mathbf{y}\right)'\mathbf{W}^{-1}
#' \left(\widehat{\mathbf{y}} - \mathbf{y}\right)} constrained by
#' \mjseqn{\mathbf{U}'\mathbf{y} = \mathbf{0}_{(r \times 1)}}:
#'
#' \mjsdeqn{\widetilde{\mathbf{y}} = \mbox{argmin}_\mathbf{y} \left(\mathbf{y} -
#' \widehat{\mathbf{y}} \right)' \mathbf{W}^{-1} \left(\mathbf{y} -
#' \widehat{\mathbf{y}} \right), \quad \mbox{s.t. } \mathbf{U}'\mathbf{y} =
#' \mathbf{0}.}
#'
#' The solution is given by
#' \mjsdeqn{\widetilde{\mathbf{y}}= \widehat{\mathbf{y}} - \mathbf{W}\mathbf{U}
#' \left(\mathbf{U}’\mathbf{WU}\right)^{-1}\mathbf{U}'\widehat{\mathbf{y}}=
#' \mathbf{M}\widehat{\mathbf{y}},}
#' where \mjseqn{\mathbf{M} = \mathbf{I}_n - \mathbf{W}\mathbf{U}\left(
#' \mathbf{U}’\mathbf{WU}\right)^{-1}\mathbf{U}’}
#' is a (\mjseqn{n \times n}) projection matrix. This solution is used by
#' \code{\link[FoReco]{htsrec}} when \code{type = "M"}.
#'
#' Denoting with \mjseqn{\mathbf{d}_{\widehat{\mathbf{y}}} = \mathbf{0} -
#' \mathbf{U}'\widehat{\mathbf{y}}} the (\mjseqn{r \times 1}) vector containing
#' the \emph{coherency} errors of the base forecasts, we can re-state the solution as
#' \mjsdeqn{\widetilde{\mathbf{y}}= \widehat{\mathbf{y}} + \mathbf{WU} \left(
#' \mathbf{U}'\mathbf{WU}\right)^{-1}\mathbf{d}_{\widehat{y}},}
#' which makes it clear that the reconciliation formula simply adjusts the
#' vector \mjseqn{\widehat{\mathbf{y}}} with a linear
#' combination -- according to the smoothing matrix
#' \mjseqn{\mathbf{L} = \mathbf{WU} \left(\mathbf{U}’\mathbf{WU}\right)^{-1}} --
#' of the coherency errors of the base forecasts.
#'
#' In addition, if the error term \mjseqn{\mathbf{\varepsilon}} is gaussian, the reconciliation
#' error \mjseqn{\widetilde{\varepsilon} = \widetilde{\mathbf{y}} - \mathbf{y}} is
#' a zero-mean gaussian vector with covariance matrix
#' \mjsdeqn{E\left( \widetilde{\mathbf{y}} - \mathbf{y}\right) \left(
#' \widetilde{\mathbf{y}} - \mathbf{y}\right)' = \mathbf{W} - \mathbf{WU}
#' \left(\mathbf{U}'\mathbf{WU}\right)^{-1}\mathbf{U}' = \mathbf{MW}.}
#'
#' Hyndman et al. (2011, see also Wickramasuriya et al., 2019) propose an
#' equivalent, alternative formulation as for the reconciled estimates, obtained
#' by GLS estimation of the model
#' \mjsdeqn{\widehat{\mathbf{y}} = \mathbf{S}\mathbf{\beta} + \mathbf{\varepsilon},}
#' where \mjseqn{\mathbf{S}} is the \emph{structural summation matrix} describing
#' the aggregation relationships operating on \mjseqn{\mathbf{y}}, and
#' \mjseqn{\mathbf{\beta}} is a subset of \mjseqn{\mathbf{y}} containing the
#' target forecasts of the bottom level series, such that
#' \mjseqn{\mathbf{y} = \mathbf{S}\mathbf{\beta}}. Since the hypotheses on
#' \mjseqn{\mathbf{\varepsilon}} remain unchanged,
#' \mjsdeqn{\widetilde{\mathbf{\beta}} = \left(\mathbf{S}'\mathbf{W}^{-1}\mathbf{S}
#' \right)^{-1}\mathbf{S}'\mathbf{W}^{-1}\widehat{\mathbf{y}}}
#' is the best linear unbiased estimate of \mjseqn{\mathbf{\beta}}, and
#' the whole reconciled forecasts vector is given by
#' \mjsdeqn{\widetilde{\mathbf{y}} = \mathbf{S}\widetilde{\mathbf{\beta}} = \mathbf{SG}
#' \widehat{\mathbf{y}},}
#' where \mjseqn{\mathbf{G} = \left(\mathbf{S}'\mathbf{W}^{-1}
#' \mathbf{S}\right)^{-1}\mathbf{S}'\mathbf{W}^{-1}}, and \mjseqn{\mathbf{M}=\mathbf{SG}}.
#' This solution is used by \code{\link[FoReco]{htsrec}} when \code{type = "S"}.
#'
#' \strong{Bounds on the reconciled forecasts}
#'
#' The user may impose bounds on the reconciled forecasts.
#' The parameter \code{bounds} permits to consider lower (\mjseqn{\mathbf{a}}) and
#' upper (\mjseqn{\mathbf{b}}) bounds like \mjseqn{\mathbf{a} \leq
#' \widetilde{\mathbf{y}} \leq \mathbf{b}} such that:
#' \mjsdeqn{ \begin{array}{c}
#' a_1 \leq \widetilde{y}_1 \leq b_1 \cr
#' \dots \cr
#' a_n \leq \widetilde{y}_n \leq b_n \cr
#' \end{array} \Rightarrow
#' \mbox{bounds} = [\mathbf{a} \; \mathbf{b}] =
#' \left[\begin{array}{cc}
#' a_1 & b_1 \cr
#' \vdots & \vdots \cr
#' a_n & b_n \cr
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
#'   \item \code{nn_type = "sntz"} ("set-negative-to-zero")
#'   \item \code{nn_type = "osqp"} (Stellato et al., 2020)
#' }
#'
#' @return
#' If the parameter \code{keep} is equal to \code{"recf"}, then the function
#' returns only the (\mjseqn{h \times n}) reconciled forecasts matrix, otherwise (\code{keep="all"})
#' it returns a list that mainly depends on what type of representation (\code{type})
#' and solution technique (\code{sol}) have been used:
#' \item{\code{recf}}{(\mjseqn{h \times n}) reconciled forecasts matrix, \mjseqn{\widetilde{\mathbf{Y}}}.}
#' \item{\code{W}}{Covariance matrix used for forecast reconciliation, \mjseqn{\mathbf{W}}.}
#' \item{\code{nn_check}}{Number of negative values (if zero, there are no values below zero).}
#' \item{\code{rec_check}}{Logical value: \code{rec_check = TRUE} when the constraints have been fulfilled.}
#' \item{\code{varf} (\code{type="direct"})}{(\mjseqn{n \times 1}) reconciled forecasts variance vector for \mjseqn{h=1}, \mjseqn{\mbox{diag}(\mathbf{MW}}).}
#' \item{\code{M} (\code{type="direct"})}{Projection matrix, \mjseqn{\mathbf{M}} (projection approach).}
#' \item{\code{G} (\code{type="S"} and \code{type="direct"})}{Projection matrix, \mjseqn{\mathbf{G}} (structural approach, \mjseqn{\mathbf{M}=\mathbf{S}\mathbf{G}}).}
#' \item{\code{S} (\code{type="S"} and \code{type="direct"})}{Cross-sectional summing matrix, \mjseqn{\mathbf{S}}.}
#' \item{\code{info} (\code{type="osqp"})}{matrix with information in columns
#' for each forecast horizon \mjseqn{h} (rows): run time (\code{run_time}),
#' number of iteration (\code{iter}), norm of primal residual (\code{pri_res}),
#' status of osqp's solution (\code{status}) and polish's status
#' (\code{status_polish}). It will also be returned with \code{nn = TRUE} if
#' a solver (see \code{nn_type}) will be used.}
#'
#' Only if \code{comb = "bu"}, the function returns \code{recf}, \code{S} and \code{M}.
#'
#' @references
#' Byron, R.P. (1978), The estimation of large social accounts matrices,
#' \emph{Journal of the Royal Statistical Society A}, 141, 3, 359-367.
#'
#' Di Fonzo, T., and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39(1), 39-57.
#'
#' Di Fonzo, T., Marini, M. (2011), Simultaneous and two-step reconciliation of
#' systems of time series: methodological and practical issues,
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)},
#' 60, 2, 143-164
#'
#' Hyndman, R.J., Ahmed, R.A., Athanasopoulos, G., Shang, H.L.(2011),
#' Optimal combination forecasts for hierarchical time series,
#' \emph{Computational Statistics & Data Analysis}, 55, 9, 2579-2589.
#'
#' Kourentzes, N., Athanasopoulos, G. (2021),
#' Elucidate structure in intermittent demand series,
#' \emph{European Journal of Operational Research}, 288, 1, pp. 141–152.
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
#' van Erven, T., Cugliari, J. (2015), Game-theoretically Optimal Reconciliation
#' of Contemporaneous Hierarchical Time Series Forecasts, in Antoniadis, A.,
#' Poggi, J.M., Brossat, X. (eds.), \emph{Modeling and Stochastic Learning for
#' Forecasting in High Dimensions}, Berlin, Springer, 297-317.
#'
#' Wickramasuriya, S.L., Athanasopoulos, G., Hyndman, R.J. (2019), Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace minimization,
#' \emph{Journal of the American Statistical Association}, 114, 526, 804-819.
#'
#' @keywords bottom-up
#' @family reconciliation procedures
#'
#' @examples
#' data(FoReco_data)
#' # monthly base forecasts
#' mbase <- FoReco2matrix(FoReco_data$base, m = 12)$k1
#' # monthly residuals
#' mres <- FoReco2matrix(FoReco_data$res, m = 12)$k1
#' obj <- htsrec(mbase, C = FoReco_data$C, comb = "shr", res = mres)
#'
#' # FoReco is able to work also with covariance matrix that are not equal
#' # across all the forecast horizon. For example, we can consider the
#' # normalized squared differences (see Di Fonzo and Marini, 2011) where
#' # Wh = diag(|yh|):
#' Wh <- lapply(split(mbase, row(mbase)), function(x) diag(abs(x)))
#'
#' # Now we can introduce the list of the covariance matrix in htsrec throught
#' # the parameter "W" and setting comb = "w".
#' objh <- htsrec(mbase, C = FoReco_data$C, W = Wh, comb = "w")
#'
#' @export
#'
#' @import Matrix osqp methods
#'
htsrec <- function(basef, comb, C, res, Ut, nb, mse = TRUE, corpcor = FALSE,
                   type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
                   nn_type = "osqp", settings = list(), bounds = NULL, W = NULL){

  if(missing(comb)){
    stop("The argument comb is not specified.", call. = FALSE)
  }else if(comb != "w"){
    W <-  NULL
  }

  UseMethod("htsrec", W)
}

#' @export
htsrec.default <- function(basef, comb, C, res, Ut, nb, mse = TRUE, corpcor = FALSE,
                           type = "M", sol = "direct", keep = "list", v = NULL, nn = FALSE,
                           nn_type = "osqp", settings = list(), bounds = NULL, W) {

  comb <- match.arg(comb, c("bu", "ols", "struc", "wls", "shr", "sam", "w"))

  type <- match.arg(type, c("M", "S"))
  keep <- match.arg(keep, c("list", "recf"))
  nn_type <- match.arg(nn_type, c("osqp", "KAnn", "fbpp", "sntz"))

  # base forecasts condition
  if (missing(basef)) {
    stop("The argument basef is not specified.", call. = FALSE)
  }
  if (NCOL(basef) == 1) {
    basef <- t(basef)
  }

  # Using Ut or C
  if (missing(C)) {
    if (missing(Ut)) {
      stop("Please, give C or Ut.", call. = FALSE)
    } else if(missing(nb)){
      hts <- hts_tools(Ut = Ut, h = 1, sparse = TRUE)
    } else {
      hts <- hts_tools(Ut = Ut, nb = nb, h = 1, sparse = TRUE)
    }
  } else {
    hts <- hts_tools(C = C, h = 1, sparse = TRUE)
  }

  n <- hts$n
  na <- hts$na
  nb <- hts$nb
  C <- hts$C
  S <- hts$S
  Ut <- hts$Ut

  if(nn){
    if(nn_type == "fbpp" | nn_type == "KAnn"){
      type = "M"
    }
  }

  if(is.null(S)){
    if (type == "S") {
      stop("Type = S needs C matrix. \nPlease consider using ut2c() to get the right C when working with Ut.", call. = FALSE)
    }

    if(nn_type == "sntz"){
      stop("nn_type = sntz needs C matrix.\nPlease consider using ut2c() to get the right C when working with Ut.", call. = FALSE)
    }

    if (comb == "bu" | comb == "struc") {
      stop("Param comb equal to bu or struc needs C matrix.\nPlease consider using ut2c() to get the right C when working with Ut.", call. = FALSE)
    }
  }

  if ((NCOL(basef) != n & comb != "bu") | ((NCOL(basef) != nb & NCOL(basef) != n) & comb == "bu")) {
    stop("Incorrect dimension of Ut, C or basef.", call. = FALSE)
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
    shr_mod <- function(x, ...) Matrix(unclass(corpcor::cov.shrink(x, verbose = FALSE, lambda.var=0, ...)))
  } else {
    shr_mod <- function(x, ...) shrink_estim(x, minT = mse)[[1]]
  }

  switch(comb,
         bu = {
           if (NCOL(basef) != nb) {
             basef_bu <- basef[, (na + 1):n]
           }else{
             basef_bu <- basef
           }

           names_bu <- if (is.null(colnames(basef))) paste("serie", 1:n, sep = "") else colnames(basef)

           if(nn){
             basef_bu <- basef_bu * (basef_bu > 0)
           }

           outf <- basef_bu %*% t(S)

           rownames(outf) <- paste("h", 1:NROW(outf), sep = "")
           if(length(names_bu) == NCOL(outf)){
             colnames(outf) <- names_bu
           }else{
             colnames(outf) <- paste("serie", 1:n, sep = "")
           }


           if (keep == "list") {
             return(list(recf = as.matrix(outf), S = S, M = S %*% cbind(matrix(0, nb, na), diag(nb))))
           } else {
             return(as.matrix(outf))
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
           if (is.null(W)) {
             stop("Please, put in option W your covariance matrix", call. = FALSE)
           }
           W <- W
         }
  )

  b_pos <- c(rep(0, na), rep(1, nb))

  if(!is.null(v)){
    keep <- "recf"
    rec_sol <- recoV(
      basef = basef, W = W, Ht = Ut, sol = sol, nn = nn, keep = keep, S = S, type = type,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type, v = v
    )
  }else if(type == "S"){
    rec_sol <- recoS(
      basef = basef, W = W, S = S, sol = sol, nn = nn, keep = keep,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
    )
  }else{
    rec_sol <- recoM(
      basef = basef, W = W, Ht = Ut, sol = sol, nn = nn, keep = keep, S = S,
      settings = settings, b_pos = b_pos, bounds = bounds, nn_type = nn_type
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


#' @export
htsrec.list <- function(basef, comb = "w", ..., W){
  if(NROW(basef) == 1){
    stop("Please put in W a matrix, not a list (the forecasts horizon is 1)", call. = FALSE)
  }

  if(NROW(basef) != length(W) | !is.list(W)){
    stop("Wh must be a list with ", NROW(basef), " matrices", call. = FALSE)
  }

  if(comb != "w"){
    comb <- "w"
  }

  obj <- Map(function(b, W){
    htsrec(basef = b,         # vector of base forecasts
           comb = "w",        # custom combination
           W = W,             # custom W
           ...)
  },
  # list of base forecasts divide by forecasts horizons
  b = split(basef, 1:NROW(basef)),
  # list of W divide by forecasts horizons
  W = W)

  out <- list()
  out$recf <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "recf")))
  colnames(out$recf) <- colnames(basef)
  rownames(out$recf) <-  paste("h",1:NROW(basef), sep="")

  out$varf <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "varf")))
  if(all(is.na(out$varf))){
    out$varf <- NULL
  }else{
    colnames(out$varf) <- colnames(basef)
    rownames(out$varf) <-  paste("h",1:NROW(basef), sep="")
  }
  out$info <- do.call(rbind, lapply(obj, function(x) extract_data(x = x, name = "info")))
  if(all(is.na(out$info))){
    out$info <- NULL
  }else{
    rownames(out$info) <-  paste("h",1:NROW(basef), sep="")
  }

  return(out)
}


extract_data <- function(x, name){
  if(is.list(x)){
    if(is.null(x[[name]])){
      NA
    }else{
      x[[name]]
    }
  }else if(name == "recf"){
    x
  }else{
    NA
  }
}

