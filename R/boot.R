#' Cross-sectional joint block bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take
#' into account the correlation between different time series (Panagiotelis et
#' al. 2023).
#'
#' @usage
#' csboot(model_list, boot_size, block_size, seed = NULL, xreg = NULL, ...)
#'
#' @param model_list A list of all the \eqn{n} base forecasts models. A
#'   \code{simulate()} function for each model has to be available and
#'   implemented according to the package
#'   \href{https://CRAN.R-project.org/package=forecast}{\pkg{forecast}},
#'   with the following mandatory parameters: \emph{object}, \emph{innov},
#'   \emph{future}, and \emph{nsim}.
#' @param boot_size The number of bootstrap replicates.
#' @param block_size Block size of the bootstrap, which is typically equivalent
#'   to the forecast horizon.
#' @param xreg An optional 3-d numeric array of dimensions
#'   (\eqn{\text{block\_size} \times n \times N_{xreg}})
#'   containing the new values of \code{xreg} to be used for forecasting.
#'   It can contain \code{NA}s.
#' @param seed An integer seed.
#' @param ... Additional arguments for the \code{simulate()} function.
#'
#' @return A list with two elements: the seed used to sample the errors and a
#'   3-d array (\eqn{\text{block\_size} \times n \times \text{boot\_size}}).
#'
#' @references
#' Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J.
#' (2023), Probabilistic forecast reconciliation: Properties, evaluation and
#' score optimisation, \emph{European Journal of Operational Research}, 306(2),
#' 693–706. \doi{http://dx.doi.org/10.1016/j.ejor.2022.07.040}
#'
#' @family Bootstrap samples
#' @family Framework: cross-sectional
#'
#' @export
csboot <- function(
  model_list,
  boot_size,
  block_size,
  seed = NULL,
  xreg = NULL,
  ...
) {
  res <- sapply(model_list, residuals)
  N <- NROW(res)
  n <- NCOL(res)

  if (!is.null(xreg)) {
    if (!is.array(xreg) || length(dim(xreg)) != 3) {
      cli_abort(
        "{.arg xreg} must be a 3-d array of dimensions (block_size x n x N_xreg).",
        call = NULL
      )
    }
    if (dim(xreg)[1] != block_size) {
      cli_abort(
        "First dimension of {.arg xreg} must equal {.arg block_size} ({block_size}).",
        call = NULL
      )
    }
    if (dim(xreg)[2] != n) {
      cli_abort(
        "Second dimension of {.arg xreg} must equal the number of series ({n}).",
        call = NULL
      )
    }
  }

  if (is.null(seed)) {
    seed <- stats::rpois(1, 1000)
  }

  index <- boot_index(
    N = N,
    boot_size = boot_size,
    block_size = block_size,
    seed = seed
  )

  fboot <- apply(
    index,
    2,
    function(id) {
      sapply(1:length(model_list), function(x) {
        if (!is.null(xreg)) {
          if (all(is.na(xreg[, x, ]))) {
            unname(simulate(
              model_list[[x]],
              innov = res[id, x],
              future = TRUE,
              nsim = length(res[id, x]),
              ...
            ))
          } else {
            unname(simulate(
              model_list[[x]],
              innov = res[id, x],
              future = TRUE,
              nsim = length(res[id, x]),
              xreg = xreg[, x, ],
              ...
            ))
          }
        } else {
          unname(simulate(
            model_list[[x]],
            innov = res[id, x],
            future = TRUE,
            nsim = length(res[id, x]),
            ...
          ))
        }
      })
    },
    simplify = FALSE
  )

  fboot <- simplify2array(fboot)
  dimnames(fboot) <- list(
    paste0("h-", 1:dim(fboot)[1]),
    paste0("s-", 1:dim(fboot)[2]),
    paste0("b-", 1:dim(fboot)[3])
  )
  if (!is.null(names(model_list))) {
    dimnames(fboot)[[2]] <- names(model_list)
  }
  return(list(sample = fboot, seed = seed))
}

#' Temporal joint block bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take
#' into account the correlation between different temporal aggregation orders
#' (Girolimetto et al. 2023).
#'
#' @usage
#' teboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL,
#'        xreg = NULL, ...)
#'
#' @param model_list A list of all the \eqn{p} base forecasts models
#'   ordered from the lowest frequency (most temporally aggregated) to the
#'   highest frequency. A \code{simulate()} function for each model has to be
#'   available and implemented according to the package
#'   \href{https://CRAN.R-project.org/package=forecast}{\pkg{forecast}},
#'   with the following mandatory parameters: \emph{object}, \emph{innov},
#'   \emph{future}, and \emph{nsim}.
#' @param block_size Block size of the bootstrap, which is typically equivalent
#'   to the forecast horizon for the most temporally aggregated series.
#' @param xreg A (\eqn{\text{block\_size}(k^\ast+m) \times N_{xreg}}) optional
#'   numeric matrix containing the new values of \code{xreg} to be used for
#'   forecasting ordered from the lowest frequency to the highest frequency.
#'   It can contains \code{NA}s.
#' @inheritParams terec
#' @inheritParams csboot
#'
#' @return A list with two elements: the seed used to sample the errors and
#'   a (\eqn{\text{boot\_size}\times (k^\ast+m)\text{block\_size}}) matrix.
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3,
#' 1134-1151. \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' @family Bootstrap samples
#' @family Framework: temporal
#'
#' @export
teboot <- function(
  model_list,
  boot_size,
  agg_order,
  block_size = 1,
  seed = NULL,
  xreg = NULL,
  ...
) {
  info <- tetools(agg_order = agg_order)

  if (length(model_list) != length(info$set)) {
    cli_abort("Incorrect {.arg model_list} dimension.", call = NULL)
  }
  names(model_list) <- paste0("k", info$set)

  res_list <- lapply(model_list, function(mod) residuals(mod))

  N <- NROW(res_list[[paste0("k", info$dim[["m"]])]])

  if (is.null(seed)) {
    seed <- stats::rpois(1, 1000)
  }

  if (!is.null(xreg)) {
    kt <- length(model_list)

    if (NROW(xreg) %% kt != 0) {
      cli_abort("Incorrect {.arg xreg} length.", call = NULL)
    } else {
      xreg <- FoReco2matrix(t(xreg), agg_order = agg_order)
    }
  }

  index <- boot_index(
    N = N,
    boot_size = boot_size,
    agg_order = info$set,
    block_size = block_size,
    seed = seed
  )

  fboot <- lapply(1:boot_size, function(i) {
    lapply(info$set, function(k) {
      cat(k, " ")
      id <- index[[paste0("k", k)]][, i]
      fit_i <- model_list[[paste0("k", k)]]
      res_vec <- res_list[[paste0("k", k)]][id]
      if (!is.null(xreg)) {
        xreg_i <- xreg[[paste0("k-", k)]]
        if (all(is.na(xreg_i))) {
          simulate(fit_i, innov = res_vec, future = TRUE, ...)
        } else {
          simulate(fit_i, innov = res_vec, future = TRUE, xreg = xreg_i, ...)
        }
      } else {
        simulate(fit_i, innov = res_vec, future = TRUE, ...)
      }
    })
  })
  fboot <- t(sapply(fboot, Reduce, f = "c"))
  return(list(sample = fboot, seed = seed))
}

#' Cross-temporal joint block bootstrap
#'
#' Joint block bootstrap for generating probabilistic base forecasts that take
#' into account the correlation between variables at different temporal
#' aggregation orders (Girolimetto et al. 2023).
#'
#' @usage
#' ctboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL,
#'        xreg = NULL, ...)
#'
#' @param model_list A list of \eqn{p} elements, one for each temporal level
#'   ordered from the lowest frequency (most temporally aggregated) to
#'   the highest frequency. Each elements is a list with the \eqn{n} base
#'   forecasts models for each cross-sectional series. A \code{simulate()}
#'   function for each model has to be available and implemented according
#'   to the package
#'   \href{https://CRAN.R-project.org/package=forecast}{\pkg{forecast}},
#'   with the following mandatory parameters: \emph{object}, \emph{innov},
#'   \emph{future}, and \emph{nsim}.
#' @param xreg An optional 3-d numeric array of dimensions
#'   (\eqn{n \times \text{boot\_size}(k^\ast+m) \times N_{xreg}})
#'   containing the new values of \code{xreg} to be used for forecasting
#'   ordered from the lowest frequency to the highest frequency (columns) for
#'   each variable (rows). It can contain \code{NA}s.
#' @inheritParams ctrec
#' @inheritParams teboot
#'
#' @return A list with two elements: the seed used to sample the errors and a
#'   list with \eqn{\text{boot\_size}} matrix of size
#'   (\eqn{n\times(k^\ast+m)\text{block\_size}}) matrix.
#'
#' @references
#' Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J. (2024),
#' Cross-temporal probabilistic forecast reconciliation: Methodological and
#' practical issues. \emph{International Journal of Forecasting}, 40, 3,
#' 1134-1151. \doi{10.1016/j.ijforecast.2023.10.003}
#'
#' @family Bootstrap samples
#' @family Framework: cross-temporal
#'
#' @export
ctboot <- function(
  model_list,
  boot_size,
  agg_order,
  block_size = 1,
  seed = NULL,
  xreg = NULL,
  ...
) {
  info <- tetools(agg_order = agg_order)

  if (length(model_list) != length(info$set)) {
    cli_abort("Incorrect {.arg model_list} dimension.", call = NULL)
  }
  names(model_list) <- paste0("k", info$set)

  res_list <- lapply(model_list, function(mod) sapply(mod, residuals))

  N <- NROW(res_list[[paste0("k", info$dim[["m"]])]])

  if (is.null(seed)) {
    seed <- stats::rpois(1, 1000)
  }

  if (!is.null(xreg)) {
    n <- length(model_list[[1]])

    if (!is.array(xreg) || length(dim(xreg)) != 3) {
      cli_abort(
        "{.arg xreg} must be a 3-d array of dimensions (n x boot_size(k*+m) x N_xreg).",
        call = NULL
      )
    }
    if (dim(xreg)[1] != n) {
      cli_abort(
        "First dimension of {.arg xreg} must equal the number of series ({n}).",
        call = NULL
      )
    }
    if (dim(xreg)[2] %% info$dim[["kt"]] != 0) {
      cli_abort(
        "Second dimension of {.arg xreg} must be a multiple of {kt}.",
        call = NULL
      )
    }

    # Convert each slice along the third dimension and store as a list of
    # FoReco2matrix results, then recombine into a 3-d structure
    xreg <- lapply(seq_len(dim(xreg)[3]), function(j) {
      FoReco2matrix(xreg[,, j], agg_order = agg_order)
    })
    # xreg is now a list of length N_xreg, each element is a list keyed by "k..."
    # Rearrange to: a list keyed by "k...", each element is an array (rows x cols x N_xreg)
    xreg <- lapply(setNames(nm = paste0("k-", info$set)), function(k) {
      slices <- lapply(xreg, function(j) j[[k]])
      simplify2array(slices) # rows x cols x N_xreg
    })
  }

  index <- boot_index(
    N = N,
    boot_size = boot_size,
    agg_order = info$set,
    block_size = block_size,
    seed = seed
  )

  fboot <- lapply(1:boot_size, function(i) {
    lapply(info$set, function(k) {
      id <- index[[paste0("k", k)]][, i]
      fit_i <- model_list[[paste0("k", k)]]
      res_mat <- res_list[[paste0("k", k)]][id, , drop = FALSE]

      if (!is.null(xreg)) {
        xreg_i <- xreg[[paste0("k-", k)]] # rows x cols x N_xreg
      }

      out <- sapply(1:length(fit_i), function(x) {
        if (!is.null(xreg)) {
          xreg_x <- xreg_i[, x, ] # rows x N_xreg matrix for series x
          if (all(is.na(xreg_x))) {
            simulate(fit_i[[x]], innov = res_mat[, x], future = TRUE, ...)
          } else {
            simulate(
              fit_i[[x]],
              innov = res_mat[, x],
              future = TRUE,
              xreg = xreg_x,
              ...
            )
          }
        } else {
          simulate(fit_i[[x]], innov = res_mat[, x], future = TRUE, ...)
        }
      })
      if (is.vector(out)) {
        out <- unname(rbind(out))
      }
      colnames(out) <- names(fit_i)
      t(out)
    })
  })

  fboot <- lapply(fboot, Reduce, f = "cbind")
  return(list(sample = fboot, seed = seed))
}


boot_index <- function(N, boot_size, agg_order, block_size, seed) {
  if (missing(agg_order)) {
    if (missing(block_size)) {
      cli_abort(
        "Argument {.arg block_size} is missing, with no default.",
        call = NULL
      )
    }

    if (is.null(seed)) {
      cli_warn("No {.arg seed} provided.", call = NULL)
    } else {
      set.seed(seed)
    }

    index <- sample(1:(N - block_size + 1), size = boot_size, replace = TRUE)
    kid <- sapply(index, function(x) {
      base::seq(x, x + block_size - 1, by = 1)
    })
    kid <- rbind(kid)
  } else {
    if (missing(block_size)) {
      block_size <- 1
    }

    if (length(agg_order) > 1) {
      kset <- agg_order
      m <- max(agg_order)
    } else {
      info <- tetools(agg_order = agg_order)
      kset <- info$set
      m <- info$dim[["m"]]
    }

    if (is.null(seed)) {
      cli_warn("No {.arg seed} provided.", call = NULL)
    } else {
      set.seed(seed)
    }

    index <- base::sample(
      1:(N - block_size + 1),
      size = boot_size,
      replace = TRUE
    )
    kid <- lapply(kset, function(k) {
      sapply(index, function(x) {
        Mk <- m / k
        base::seq(Mk * (x - 1) + 1, Mk * (x + block_size - 1), by = 1)
      })
    })
    kid <- lapply(kid, rbind)
    names(kid) <- paste0("k", kset)
  }
  return(kid)
}
