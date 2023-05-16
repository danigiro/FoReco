# Bootstrap index
#
# @param N The sample size.
# @param boot_size The number of bootstrap replicates.
# @param m Highest available sampling frequency per seasonal cycle (max. order
# of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
# of \mjseqn{m}.
# @param h Block size of the bootstrap in the cross-sectional framework, or the forecast horizon
# for the most temporally aggregated series in the temporal/cross-temporal framework.
# @param seed An integer seed.
#
# @return The indeces' matrix in the cross-sectional framework, or a list of indeces' matrix in
# temporal and cross-temporal
#
boot_index <- function(N, boot_size, m, h, seed){
  if(missing(m)){
    if(missing(h)){
      stop("Provide h for cross-sectional index")
    }

    if(is.null(seed))
      warning("no seed")
    else
      set.seed(seed)

    index <- sample(1:(N-h+1), size = boot_size, replace = TRUE)
    kid <- sapply(index, function(x){
      base::seq(x, x+h-1, by = 1)
    })
  }else{
    if(missing(h)){
      h <- 1
    }
    info <- thf_tools(m = m)
    kset <- info$kset
    m <- info$m

    if(is.null(seed))
      warning("no seed")
    else
      set.seed(seed)

    index <- base::sample(1:(N-h+1), size = boot_size, replace = TRUE)
    kid <- lapply(kset, function(k) sapply(index, function(x){
      Mk <- m/k
      base::seq(Mk*(x-1)+1, Mk*(x+h-1), by = 1)
    }))
    kid <- lapply(kid, rbind)
    names(kid) <- paste0("k", kset)
  }
  return(kid)
}

#' Cross-sectional Joint Bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take into account
#' the correlation between different time series (Panagiotelis et al. 2023).
#'
#' @param fit A list of \mjseqn{n} base forecast models. It is important to note that the models
#' must have the \code{simulate()} function available and implemented as with
#' the package \pkg{forecast}, with the following mandatory parameters:
#' \emph{object}, \emph{innov}, \emph{future}, and \emph{nsim}.
#' @param boot_size The number of bootstrap replicates.
#' @param h Block size of the bootstrap, which is typically equivalent to the forecast horizon.
#' @param seed An integer seed.
#'
#' @return A list with two elements: the seed used to sample the errors and a 3-d array
#' (\mjseqn{boot\_size\times n \times h})
#'
#' @references
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. & Hyndman, R. J. (2023),
#' Probabilistic forecast reconciliation: Properties, evaluation and score optimisation,
#' \emph{European Journal of Operational Research} 306(2), 693–706.
#'
#' @family bootstrap
#' @importFrom stats residuals simulate
#'
#' @export
boot_cs <- function(fit, boot_size, h, seed = NULL){
  res <- sapply(fit, residuals)
  N <- NROW(res)

  if(is.null(seed))
    seed <- stats::rpois(1, 1000)

  index <- boot_index(N = N, boot_size = boot_size, h = h, seed = seed)
  fboot <- apply(index, 2, function(id) sapply(1:length(fit), function(x){
    unname(simulate(fit[[x]], innov = res[id, x], future = TRUE, nsim = length(res[id, x])))
  }), simplify = FALSE)
  fboot <- aperm(simplify2array(fboot), c(3, 2, 1))
  return(list(sample = fboot,
              seed = seed))
}

#' Temporal Joint Bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take into account
#' the correlation between different temporal aggregation order (Girolimetto et al. 2023).
#'
#' @param fit A list of \mjseqn{(k^\ast+m)} base forecast models ordered as
#' [lowest_freq' ...  highest_freq']'. It is important to note that the models
#' must have the \code{simulate()} function available and implemented as with
#' the package \pkg{forecast}, with the following mandatory parameters:
#' \emph{object}, \emph{innov}, \emph{future}, and \emph{nsim}.
#' @param boot_size The number of bootstrap replicates.
#' @param h Forecast horizon for the most temporally aggregated series.
#' @param seed An integer seed.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#'
#' @return A list with two elements: the seed used to sample the errors and
#' a (\mjseqn{boot\_size\times h(k^\ast+m)}) matrix
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., & Hyndman, R. J. (2023),
#' Cross-temporal Probabilistic Forecast Reconciliation,
#' \doi{10.48550/arXiv.2303.17277}.
#'
#' @family bootstrap
#' @importFrom stats residuals simulate
#'
#' @export
boot_te <- function(fit, boot_size, m, h = 1, seed = NULL){
  info <- thf_tools(m = m)
  m <- info$m
  kt <- info$kt
  kset <- info$kset

  if(length(fit)!=length(kset)){
    stop("fit != kset")
  }
  names(fit) <- paste0("k",kset)

  res_list <- lapply(fit, function(mod) residuals(mod))

  N <- NROW(res_list[[paste0("k",m)]])

  if(is.null(seed))
    seed <- stats::rpois(1, 1000)

  index <- boot_index(N = N, boot_size = boot_size, m = m, h = h, seed = seed)

  fboot <- lapply(1:boot_size, function(i){
    lapply(kset, function(k){
      id <- index[[paste0("k",k)]][,i]
      fit_i <- fit[[paste0("k",k)]]
      res_vec <- res_list[[paste0("k",k)]][id]
      simulate(fit_i, innov = res_vec, future = TRUE)
    })
  })
  fboot <- t(sapply(fboot, Reduce, f = "c"))
  return(list(sample = fboot,
              seed = seed))
}

#' Cross-temporal Joint Bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take into account
#' the correlation between different time series and temporal aggregation order at the same time
#' (Girolimetto et al. 2023).
#'
#' @param fit A list of \mjseqn{n} elements. Each elements is a list with the \mjseqn{(k^\ast+m)}
#' base forecast models ordered as [lowest_freq' ...  highest_freq']' of the cross-sectional
#' series. It is important to note that the models must have the \code{simulate()}
#' function available and implemented as with the package \pkg{forecast}, with
#' the following mandatory parameters: \emph{object}, \emph{innov}, \emph{future}, and \emph{nsim}.
#' @param boot_size The number of bootstrap replicates.
#' @param h Forecast horizon for the most temporally aggregated series.
#' @param seed An integer seed.
#' @param m Highest available sampling frequency per seasonal cycle (max. order
#' of temporal aggregation, \mjseqn{m}), or a subset of \mjseqn{p} factors
#' of \mjseqn{m}.
#'
#' @return A list with two elements: the seed used to sample the errors and
#' a (\mjseqn{boot\_size\times hn(k^\ast+m)}) matrix
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., & Hyndman, R. J. (2023),
#' Cross-temporal Probabilistic Forecast Reconciliation,
#' \doi{10.48550/arXiv.2303.17277}.
#'
#' @family bootstrap
#' @importFrom stats residuals simulate
#'
#' @export
boot_ct <- function(fit, boot_size, m, h = 1, seed = NULL){
  info <- thf_tools(m = m)
  m <- info$m
  kt <- info$kt
  kset <- info$kset

  if(length(fit)!=length(kset)){
    stop("fit != kset")
  }
  names(fit) <- paste0("k",kset)

  res_list <- lapply(fit, function(mod) sapply(mod, residuals))

  N <- NROW(res_list[[paste0("k",m)]])

  if(is.null(seed))
    seed <- stats::rpois(1, 1000)

  index <- boot_index(N = N, boot_size = boot_size, m = kset, h = h, seed = seed)

  fboot <- lapply(1:boot_size, function(i){
    lapply(kset, function(k){
      id <- index[[paste0("k",k)]][,i]
      fit_i <- fit[[paste0("k",k)]]
      res_mat <- res_list[[paste0("k",k)]][id,,drop = FALSE]
      out <- sapply(1:length(fit_i), function(x){
        simulate(fit_i[[x]], innov = res_mat[, x], future = TRUE)
      })
      if(is.vector(out)){
        out <- unname(rbind(out))
      }
      colnames(out) <- names(fit_i)
      out
    })
  })

  fboot <- lapply(fboot, Reduce, f = "rbind")
  return(list(sample = fboot,
              seed = seed))
}
