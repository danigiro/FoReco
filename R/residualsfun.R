#' Re-arrange the multi-step residuals
#'
#' @param list_res a list of \mjseqn{H} multi-step residuals. Each element
#' of the list can be a vector (univariate time series) or a matrix (multivariate time series).
#'
#' @details
#' Let \mjseqn{Z_t}, \mjseqn{t=1,\dots,T}, be a univariate time series. We can define the multi-step
#' residuals such us
#' \mjsdeqn{\widehat{\varepsilon}_{h,t} = Z_{t+h} - \widehat{Z}_{t+h|t} \qquad h \le t \le T-h}
#' where \mjseqn{\widehat{Z}_{t+h|t}} is the \mjseqn{h}-step fitted value, calculated as the \mjseqn{h}-step ahead
#' forecast given the time \mjseqn{t}. Given the list of errors at different step
#' (\mjseqn{[\widehat{\varepsilon}_{1,1}, \; \dots, \; \widehat{\varepsilon}_{1,T}]}, ..., \mjseqn{[\widehat{\varepsilon}_{H,1}, \; \dots, \; \widehat{\varepsilon}_{H,T}]})
#' this function returns a \mjseqn{T}-vector with the residuals, organized in the following way:
#' \mjsdeqn{[\varepsilon_{1,1} \; \varepsilon_{2,2} \; \dots \; \varepsilon_{H,H} \; \varepsilon_{1,H+1} \; \dots \; \varepsilon_{H,T-H}]'}
#' Same idea can be apply for a multivariate time series.
#'
#' @return A vector or a matrix of multi-step residuals
#'
#' @family utilities
#' @importFrom stats tsp<-
#'
#' @export
arrange_hres <- function(list_res){
  if(!is.list(list_res)){
    warning("res has to be a list")
    return(list_res)
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


#' Arrange temporal and cross-temporal residuals in a matrix form
#'
#' @param res (\mjseqn{n \times N(k^\ast + m)}) matrix or (\mjseqn{N(k^\ast + m)})
#' vector containing the residuals at all the temporal frequencies
#' ordered [lowest_freq' ...  highest_freq']'.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#'
#' @return a (\mjseqn{N \times n(k^\ast + m)}) matrix
#' (if the input res is a vector then \mjseqn{n = 1})
#'
#' @family utilities
#'
#' @export
#'
#' @examples
#' data(FoReco_data)
#' mat <- residuals_matrix(FoReco_data$res, m = 12)
#'
residuals_matrix <- function(res, m){
  info <- thf_tools(m = m)
  kt <- info$kt
  kset <- info$kset

  if(is.vector(res)){
    if (length(res) %% kt != 0) {
      stop("res vector has a number of row not in line with frequency of the series", call. = FALSE)
    }

    N <- length(res) / kt
    DN <- Dmat(h = N, m = kset, n = 1)
    matrix(DN %*% res, nrow = N, byrow = T)
  }else{
    n <- NROW(res)

    if (NCOL(res) %% kt != 0) {
      stop("res has a number of columns not in line with frequency of the series", call. = FALSE)
    }

    N <- NCOL(res) / kt

    DN <- Dmat(h = N, m = kset, n = n)
    matrix(DN %*% as.vector(t(res)), nrow = N, byrow = TRUE)
  }
}
