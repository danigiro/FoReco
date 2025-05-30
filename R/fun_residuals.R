#' One-step and multi-step residuals
#'
#' @description
#' These functions can be used to arrange residuals to reconcile temporal or
#' cross-temporal forecasts.
#'
#' [res2matrix] takes as input a set of temporal and cross-temporal residuals and
#' re-organizes them into a matrix where the rows correspond to different forecast
#' horizons, capturing the temporal dimension. Meanwhile, the columns are ordered
#' based on the specific arrangement as described in Di Fonzo and Girolimetto (2023).
#'
#' @inheritParams ctrec
#' @param res A (\eqn{n \times N(k^\ast+m)}) numeric matrix (cross-temporal framework)
#' or an (\eqn{N(k^\ast+m) \times 1}) numeric vector (temporal framework) representing
#' the in-sample residuals at all the temporal frequencies ordered from the lowest
#' frequency to the highest frequency (columns) for each variable (rows).
#'
#' @return [res2matrix] returns a (\eqn{N \times n(k^\ast + m)}) matrix, where \eqn{n = 1}
#' for the temporal framework.
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3, 1134-1151.
#' \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' @examples
#' h <- 10
#' agg_order <- 4
#' tmp <- tetools(agg_order)
#' kt <- tmp$dim["kt"]
#'
#' # Simulate vector (temporal case)
#' vec <- rnorm(kt*h)
#' out <- res2matrix(vec, agg_order) # matrix h x kt
#'
#' # Simulate (n x kt) matrix (cross-temporal case) with n = 3
#' mat <- rbind(rnorm(kt*h), rnorm(kt*h), rnorm(kt*h))
#' out <- res2matrix(mat, agg_order) # matrix h x (3*kt)
#'
#' @family Utilities
#' @rdname residuals
#'
#' @export
res2matrix <- function(res, agg_order){
  kset <- tetools(agg_order = agg_order)$set

  if(is.vector(res)){
    if(length(res) %% sum(max(kset)/kset) != 0){
      cli_abort("Incorrect {.arg res} length.", call = NULL)
    }

    N <- length(res) / sum(max(kset)/kset)
    vec2hmat(vec = res, h = N, kset = kset)
  }else{
    n <- NROW(res)
    if(NCOL(res) %% sum(max(kset)/kset) != 0){
      cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
    }

    mat2hmat(res, h = NCOL(res) / sum(max(kset)/kset), kset = kset, n = n)
  }
}


#' @description
#' [arrange_hres] takes as input a list of multi-step residuals and is
#' designed to organize them in accordance with their time order (Girolimetto et al.
#' 2023). When applied, this function ensures that the sequence of multi-step
#' residuals aligns with the chronological order in which they occurred.
#'
#' @param list_res A list of \eqn{H} multi-step residuals. Each element in the list
#' can be either a (\eqn{T \times 1}) vector (temporal framework) or a (\eqn{T \times n})
#' matrix (cross-temporal framework).
#'
#' @details
#' Let \eqn{Z_t}, \eqn{t=1,\dots,T}, be a univariate time series. We can define the multi-step
#' residuals such us
#' \deqn{\widehat{\varepsilon}_{h,t} = Z_{t+h} - \widehat{Z}_{t+h|t} \qquad h \le t \le T-h}
#' where \eqn{\widehat{Z}_{t+h|t}} is the \eqn{h}-step fitted value, calculated as the \eqn{h}-step ahead
#' forecast condition to the information up to time \eqn{t}. Given the list of errors at different steps
#' \deqn{\left([\widehat{\varepsilon}_{1,1}, \; \dots, \; \widehat{\varepsilon}_{1,T}], \dots, [\widehat{\varepsilon}_{H,1}, \; \dots, \; \widehat{\varepsilon}_{H,T}]\right),}
#' [arrange_hres] returns a \eqn{T}-vector with the residuals, organized in the following way:
#' \deqn{[\varepsilon_{1,1} \; \varepsilon_{2,2} \; \dots \; \varepsilon_{H,H} \; \varepsilon_{1,H+1} \; \dots \; \varepsilon_{H,T-H}]'}
#' A similar organisation can be apply to a multivariate time series.
#'
#' @return [arrange_hres] returns a (\eqn{N(k^\ast+m) \times 1}) vector (temporal framework)
#' or a (\eqn{n \times N(k^\ast+m)}) matrix  (cross-temporal framework) of multi-step residuals.
#'
#' @examples
#' # Input: 4 (forecast horizons) vectors with 4*10 elements
#' input <-  list(rnorm(4*10), rnorm(4*10), rnorm(4*10), rnorm(4*10))
#' # Output: 1 vector with 4*10 elements
#' out <- arrange_hres(input)
#'
#' # Matrix version
#' input <-  list(matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10),
#'                matrix(rnorm(4*10*3), 4*10), matrix(rnorm(4*10*3), 4*10))
#' out <- arrange_hres(input)
#'
#' @rdname residuals
#'
#' @export
arrange_hres <- function(list_res){
  if(!is.list(list_res)){
    cli_abort("{.arg list_res} is not a list.", call = NULL)
  }

  if(is.list(list_res) & length(list_res)<2){
    return(list_res[[1]])
  }

  out <- list_res[[1]]
  tsp(out) <- NULL
  H <- length(list_res)
  for(h in 2:H){
    outh <- list_res[[h]]
    id <- seq(h, by = H, length.out = NROW(out)/H)
    if(is.vector(out)){
      out[id] <- outh[id]
    }else{
      out[id,] <- outh[id,]
    }
  }
  return(out)
}


