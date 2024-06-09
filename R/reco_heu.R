#' Iterative cross-temporal reconciliation
#'
#' This function performs the iterative procedure described in Di Fonzo and Girolimetto (2023),
#' which produces cross-temporally reconciled forecasts by alternating forecast
#' reconciliation along one single dimension (either cross-sectional or temporal)
#' at each iteration step.
#'
#' @usage
#' iterec(base, cslist, telist, res = NULL, itmax = 100, tol = 1e-5,
#'        type = "tcs", norm = "inf", verbose = TRUE)
#'
#' @inheritParams ctrec
#' @param cslist A list of elements for the cross-sectional reconciliation.
#' See [csrec] for a complete list (excluded \code{base} and \code{res}).
#' @param telist A list of elements for the temporal reconciliation.
#' See [terec] for a complete list (excluded \code{base} and \code{res}).
#' @param itmax Max number of iteration (\code{100}, \emph{default}).
#' @param tol Convergence tolerance (\code{1e-5}, \emph{default}).
#' @param type A string specifying the uni-dimensional reconciliation order:
#' temporal and then cross-sectional ("\code{tcs}") or cross-sectional and
#' then temporal ("\code{cst}").
#' @param norm Norm used to calculate the temporal and the cross-sectional
#' incoherence: infinity norm ("\code{inf}", \emph{default}), one norm ("\code{one}"), and
#' 2-norm ("\code{two}").
#' @param verbose If \code{TRUE}, reconciliation information are printed.
#'
#' @inherit ctrec return
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' @examples
#' set.seed(123)
#' # (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
#' base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#' # (3 x 70) in-sample residuals matrix (simulated)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#'
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' m <- 4 # from quarterly to annual temporal aggregation
#'
#' rite <- iterec(base = base,
#'                cslist = list(agg_mat = A, comb = "shr"),
#'                telist = list(agg_order = m, comb = "wlsv"),
#'                res = res)
#'
#' @family Framework: cross-temporal
#' @export
iterec <- function(base, cslist, telist, res = NULL, itmax = 100, tol = 1e-5,
                   type = "tcs", norm = "inf", verbose = TRUE){

  type <- match.arg(type, c("tcs", "cst"))
  norm <- match.arg(norm, c("one", "inf", "two"))

  # Check if 'base' is provided
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  # Check cslist and telist
  if(missing(cslist) | missing(telist)){
    cli_abort("Argument {.arg cslist} and/or {.arg telist} are missing,
              with no default.", call = NULL)
  }

  # Check telist: agg_order
  if(any(names(telist)=="agg_order")){
    kset <- tetools(agg_order = telist$agg_order)$set
    h <- NCOL(base)/sum(kset)

    # Check dimension of base and res
    if(NCOL(base) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
    }

    if(!is.null(res) && NCOL(res) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
    }

    te_const_mat <- tetools(agg_order = telist$agg_order, fh = h, tew = telist$tew)$cons_mat
  }else{
    cli_abort("Argument {.arg agg_order} is missing in {.arg telist},
              with no default.", call = NULL)
  }

  # Check cslist: agg_mat or const_mat
  if(any(names(cslist) %in% c("agg_mat", "cons_mat"))){
    if(any(names(cslist)=="agg_mat")){
      cs_const_mat <- cstools(agg_mat = cslist$agg_mat)$cons_mat
    }else{
      cs_const_mat <- cstools(cons_mat = cslist$cons_mat)$cons_mat
    }
    n <- NCOL(cs_const_mat)
    # Check dimension of base and res
    if(NROW(base) != n){
      cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
    }

    if(!is.null(res) && NROW(res) != n){
      cli_abort("Incorrect {.arg res} rows dimension.", call = NULL)
    }
  }else{
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing in {.arg cslist},
              with no default.", call = NULL)
  }

  # Check comb
  if(!any(names(telist)=="comb") & verbose){
    cli_alert_warning(c("Argument {.arg comb} is missing in {.arg telist}. ",
                        "Default option: {.strong {.code ols}}."))
    telist$comb <- "ols"
  }

  if(!any(names(cslist)=="comb") & verbose){
    cli_alert_warning(c("Argument {.arg comb} is missing in {.arg cslist}. ",
                        "Default option: {.strong {.code ols}}."))
    cslist$comb <- "ols"
  }

  if(type == "tcs"){
    step1_func <- testep
    step2_func <- csstep
    rorder <- c("te", "cs")
  }else{
    step1_func <- csstep
    step2_func <- testep
    rorder <- c("cs", "te")
  }
  names_label <- c("cs" = "Cross-sectional",
                   "te" = "Temporal")
  flag <- -1
  tmp <- base
  dmat <- array(NA, dim = c(2, 3, itmax+1, 2),
                dimnames = list(c("cs", "te"),
                                c("inf", "one", "two"),
                                NULL,
                                paste0("Step ", c(1,2), " - ",rorder)))
  dmat[, , 1,2] <- dmat[, , 1,1] <- discrepancy(base, cs_const_mat = cs_const_mat,
                                                te_const_mat = te_const_mat)

  pwidth <- max(nchar(round(dmat[, norm, 1,1], digits = 2)), 15)
  if(verbose){
    cli_rule("{.strong Iterative heuristic cross-temporal forecast reconciliation}")
    cli_text("{.emph Legend: i\u00a0=\u00a0iteration; s\u00a0=\u00a0step. Norm\u00a0=\u00a0\"{norm}\".}")
    cat("\n")
    #cli_rule()
    #cat("-------------------------------------------------------------------\n",
    # "Iterative heuristic cross-temporal forecast reconciliation\n",
    # "Legend: i = iteration; s = step; te = temporal; \ncs = cross-sectional;
    #  ct = cross-temporal.\n",
    # "-------------------------------------------------------------------\n", sep = "")

    c1r1 <- gsub(" ", "\u00a0", format("i.s", width = nchar(itmax)+2, justify = "right"))
    c2r1 <- gsub(" ", "\u00a0", format(names_label[rorder[1]], width = pwidth, justify = "right"))
    c3r1 <- gsub(" ", "\u00a0", format(names_label[rorder[2]], width = pwidth, justify = "right"))
    cli_text("{.strong {c1r1}} | {.strong {c2r1}} | {.strong {c3r1}} |\n")

    #cat(rule(width = 8+(nchar(itmax)+2)+(pwidth*2)), "\n")
    cat(#format("i.s", width = nchar(itmax)+2, justify = "right"), " | ",
      #format(rorder[1], width = pwidth, justify = "right"), " | ",
      #format(rorder[2], width = pwidth, justify = "right"), " |\n",
      formatC(0, width = nchar(itmax)+2, format = "f", digits = 0), " | ",
      paste(fnumber(dmat[rorder, norm, 1,1], pwidth), collapse = " | "), " |\n", sep = "")
  }

  csc_list <- cscov_list(cslist = cslist, kset = kset, res = res, n = NROW(base))
  tec_list <- tecov_list(telist = telist, n = NROW(base), res = res)

  for(i in 1:itmax){
    tmp <- step1_func(base = tmp, telist = telist, cslist = cslist, kset = kset,
                      csc_list = csc_list, tec_list = tec_list)
    dmat[, , i+1,1] <- discrepancy(tmp, cs_const_mat = cs_const_mat, te_const_mat = te_const_mat)
    if(verbose){
      cat(formatC(i+0.1, width = nchar(itmax)+2, format = "f", digits = 1), " | ",
          paste(fnumber(dmat[rorder, norm, i+1,1], pwidth), collapse = " | "), " |\n", sep = "")
    }

    tmp <- step2_func(base = tmp, telist = telist, cslist = cslist, kset = kset,
                      csc_list = csc_list, tec_list = tec_list)
    dmat[, , i+1,2] <- discrepancy(tmp, cs_const_mat = cs_const_mat, te_const_mat = te_const_mat)
    if(verbose){
      cat(formatC(i+0.2, width = nchar(itmax)+2, format = "f", digits = 1), " | ",
          paste(fnumber(dmat[rorder, norm, i+1,2], pwidth), collapse = " | "), " |\n", sep = "")
    }

    if(all(dmat[, norm, i+1,2] < tol)){
      flag <- 0
      break
    }
  }
  dmat <- dmat[,,1:(i+1),, drop = FALSE]
  flag <- check_results(flag = flag, dmat = dmat, i = i,
                        itmax = itmax, norm = norm, verbose = verbose)
  if(verbose){
    cli_rule()
    #cat("-------------------------------------------------------------------\n")
  }

  out <- as.matrix(tmp)
  colnames(out) <- namesTE(kset = kset, h = h)
  rownames(out) <- namesCS(n = NROW(out), names_vec = rownames(base))
  attr(out, "FoReco") <- list2env(list(discrepancy = dmat,
                                       flag = flag, ite = i, norm = norm,
                                       framework = "Cross-temporal",
                                       forecast_horizon = h,
                                       comb = c("cs" = cslist$comb,
                                                "te" = telist$comb)[rorder],
                                       te_set = kset,
                                       cs_n = n,
                                       rfun = "iterec"))
  return(out)
}

#' Heuristic cross-temporal reconciliation
#'
#' [tcsrec] replicates the procedure by Kourentzes and Athanasopoulos (2019):
#' (i) for each time series the forecasts at any temporal aggregation order are
#' reconciled using temporal hierarchies; (ii) time-by-time cross-sectional
#' reconciliation is performed; and (iii) the projection matrices obtained at
#' step (ii) are then averaged and used to cross-sectionally reconcile the
#' forecasts obtained at step (i). In [cstrec], the order of application of the
#' two reconciliation steps (temporal first, then cross-sectional), is inverted
#' compared to [tcsrec] (Di Fonzo and Girolimetto, 2023).
#'
#' @inheritParams iterec
#' @param avg If \code{avg = "KA"} (\emph{default}), the final projection
#' matrix \eqn{\mathbf{M}} is the one proposed by Kourentzes and
#' Athanasopoulos (2019), otherwise it is calculated as simple average of
#' all the involved projection matrices at step 2 of the procedure (see
#' Di Fonzo and Girolimetto, 2023).
#'
#' @usage
#' # First-temporal-then-cross-sectional forecast reconciliation
#' tcsrec(base, cslist, telist, res = NULL, avg = "KA")
#'
#' @section Warning:
#' The two-step heuristic reconciliation allows considering
#' non negativity constraints only in the first step. This means that non-negativity
#' is not guaranteed in the final reconciled values.
#'
#' @inherit ctrec return
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Cross-temporal forecast reconciliation:
#' Optimal combination method and heuristic alternatives, \emph{International Journal
#' of Forecasting}, 39, 1, 39-57. \doi{10.1016/j.ijforecast.2021.08.004}
#'
#' Kourentzes, N. and Athanasopoulos, G. (2019), Cross-temporal coherent forecasts
#' for Australian tourism, \emph{Annals of Tourism Research}, 75, 393-409.
#' \doi{10.1016/j.annals.2019.02.001}
#'
#' @examples
#' set.seed(123)
#' # (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
#' base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#' # (3 x 70) in-sample residuals matrix (simulated)
#' res <- rbind(rnorm(70), rnorm(70), rnorm(70))
#'
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' m <- 4 # from quarterly to annual temporal aggregation
#'
#' rtcs <- tcsrec(base = base,
#'                cslist = list(agg_mat = A, comb = "shr"),
#'                telist = list(agg_order = m, comb = "wlsv"),
#'                res = res)
#'
#' rcst <- tcsrec(base = base,
#'                cslist = list(agg_mat = A, comb = "shr"),
#'                telist = list(agg_order = m, comb = "wlsv"),
#'                res = res)
#'
#' @rdname heuristic-reco
#'
#' @family Framework: cross-temporal
#' @export
tcsrec <- function(base, cslist, telist, res = NULL, avg = "KA"){

  #avg <- match.arg(avg, c("KA", "wmean"))
  # Check if 'base' is provided
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  # Check cslist and telist
  if(missing(cslist) | missing(telist)){
    cli_abort("Argument {.arg cslist} and/or {.arg telist} are missing,
              with no default.", call = NULL)
  }


  # Check telist: agg_order
  if(any(names(telist)=="agg_order")){
    kset <- tetools(agg_order = telist$agg_order)$set
    h <- NCOL(base)/sum(kset)

    # Check dimension of base and res
    if(NCOL(base) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
    }

    if(!is.null(res) && NCOL(res) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
    }

  }else{
    cli_abort("Argument {.arg agg_order} is missing in {.arg telist},
              with no default.", call = NULL)
  }

  # Check cslist: agg_mat or const_mat
  if(any(!(names(cslist) %in% c("agg_mat", "cons_mat")))){
    if(any(names(cslist)=="agg_mat")){
      n <- cstools(agg_mat = cslist$agg_mat)$dim[["n"]]
    }else{
      n <- cstools(cons_mat = cslist$cons_mat)$dim[["n"]]
    }

    # Check dimension of base and res
    if(NROW(base) != n){
      cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
    }

    if(!is.null(res) && NROW(res) != n){
      cli_abort("Incorrect {.arg res} rows dimension.", call = NULL)
    }
  }else{
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing in {.arg cslist},
              with no default.", call = NULL)
  }

  # Check comb
  if(!any(names(telist)=="comb")){
    telist$comb <- "ols"
  }

  if(!any(names(cslist)=="comb")){
    cslist$comb <- "ols"
  }

  step1 <- testepM(base = base, telist = telist, res = res)
  mM <- csmeanM(cslist = cslist, kset = kset, res = res, avg = avg)

  out <- as.matrix(mM %*% step1)
  colnames(out) <- namesTE(kset = kset, h = h)
  rownames(out) <- namesCS(n = NROW(out), names_vec = rownames(base))
  attr(out, "FoReco") <- list2env(list(proj_mean = mM,
                                       framework = "Cross-temporal",
                                       forecast_horizon = h,
                                       comb = c("te" = telist$comb,
                                                "cs" = cslist$comb),
                                       te_set = kset,
                                       cs_n = n,
                                       rfun = "tcsrec"))
  return(out)
}

#' @rdname heuristic-reco
#'
#' @usage
#' # First-cross-sectional-then-temporal forecast reconciliation
#' cstrec(base, cslist, telist, res = NULL)
#'
#' @export
cstrec <- function(base, cslist, telist, res = NULL){

  # Check if 'base' is provided
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  # Check cslist and telist
  if(missing(cslist) | missing(telist)){
    cli_abort("Argument {.arg cslist} and/or {.arg telist} are missing,
              with no default.", call = NULL)
  }


  # Check telist: agg_order
  if(any(names(telist)=="agg_order")){
    kset <- tetools(agg_order = telist$agg_order)$set
    h <- NCOL(base)/sum(kset)

    # Check dimension of base and res
    if(NCOL(base) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
    }

    if(!is.null(res) && NCOL(res) %% sum(kset) != 0){
      cli_abort("Incorrect {.arg res} columns dimension.", call = NULL)
    }

  }else{
    cli_abort("Argument {.arg agg_order} is missing in {.arg telist},
              with no default.", call = NULL)
  }

  # Check cslist: agg_mat or const_mat
  if(any(!(names(cslist) %in% c("agg_mat", "cons_mat")))){
    if(any(names(cslist)=="agg_mat")){
      n <- cstools(agg_mat = cslist$agg_mat)$dim[["n"]]
    }else{
      n <- cstools(cons_mat = cslist$cons_mat)$dim[["n"]]
    }

    # Check dimension of base and res
    if(NROW(base) != n){
      cli_abort("Incorrect {.arg base} rows dimension.", call = NULL)
    }

    if(!is.null(res) && NROW(res) != n){
      cli_abort("Incorrect {.arg res} rows dimension.", call = NULL)
    }
  }else{
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing in {.arg cslist},
              with no default.", call = NULL)
  }

  # Check comb
  if(!any(names(telist)=="comb")){
    telist$comb <- "ols"
  }

  if(!any(names(cslist)=="comb")){
    cslist$comb <- "ols"
  }

  step1 <- csstepM(base = base, cslist = cslist, kset = kset, res = res)
  mM <- temeanM(telist = telist, res = res, n = NROW(base))

  h_pos <- rep(rep(1:h, length(kset)), rep((max(kset)/kset), each = h))
  out <- matrix(NA, NROW(base), NCOL(base))
  for(i in 1:h){
    out[, i==h_pos] <- as.vector(step1[, i==h_pos, drop = FALSE]%*%t(mM))
  }
  colnames(out) <- namesTE(kset = kset, h = h)
  rownames(out) <- namesCS(n = NROW(out), names_vec = rownames(base))
  attr(out, "FoReco") <- list2env(list(proj_mean = mM,
                                       framework = "Cross-temporal",
                                       forecast_horizon = h,
                                       comb = c("cs" = cslist$comb,
                                                "te" = telist$comb),
                                       te_set = kset,
                                       cs_n = n,
                                       rfun = "cstrec"))
  return(out)
}

# Projection matrix mean (cross-sectional framework)
csmeanM <- function(cslist, res = NULL, kset, avg = "KA", ...){
  if(is.null(res)){
    do.call("csprojmat", cslist)
  }else{
    input_list <- NULL
    input_list[as.character(kset)] <- list(cslist)
    if(!is.null(res)){
      res <- mat2list(res, kset)
    }

    M <- lapply(as.character(kset), function(k){
      input_list[[k]]$res <- res[[k]]
      do.call("csprojmat", input_list[[k]])
    })

    if(avg == "KA"){
      return(Reduce("+", M)/length(M))
    }else{
      Mw <- mapply(function(a, A) a * A, A = M, a = split(kset, 1:length(kset)), SIMPLIFY = FALSE)
      return(Reduce("+", Mw)/sum(kset))
    }
  }
}

# Projection matrix mean (temporal framework)
temeanM <- function(telist, res = NULL, n, ...){
  if(is.null(res)){
    do.call("csprojmat", telist)
  }else{
    input_list <- NULL
    input_list[1:n] <- list(telist)
    if(!is.null(res)){
      res <- split(res, 1:NROW(res))
    }

    M <- sapply(1:n, function(i){
      input_list[[i]]$res <- res[[i]]
      do.call("teprojmat", input_list[[i]])
    })
    return(Reduce("+", M)/length(M))
  }
}

# Cross-sectional covariance matrices: list
# cscov_list <- function(cslist, res = NULL, kset, n, ...){
#   input_list <- NULL
#   cslist$n <- n
#   input_list[as.character(kset)] <- list(cslist)
#   if(!is.null(res)){
#     res <- mat2list(res, kset)
#   }
#
#   out <- lapply(as.character(kset), function(k){
#     input_list[[k]]$res <- res[[k]]
#     do.call("cscov", input_list[[k]])
#   })
#   names(out) <- as.character(kset)
#   return(out)
# }
cscov_list <- function(cslist, res = NULL, kset, n, ...){
  input_list <- NULL
  cslist$n <- n
  input_list[as.character(kset)] <- list(cslist)
  if(is.null(cslist$comb) || cslist$comb %in% c("ols", "str")){
    out <- rep(cslist$comb, length(kset))
  }else{
    res <- mat2list(res, kset)
    out <- lapply(as.character(kset), function(k){
      input_list[[k]]$res <- res[[k]]
      do.call("cscov", input_list[[k]])
    })
  }

  names(out) <- as.character(kset)
  return(out)
}

# Temporal covariance matrices: list
tecov_list <- function(telist, res = NULL, n, ...){
  input_list <- NULL
  input_list[1:n] <- list(telist)
  if(is.null(telist$comb) || telist$comb %in% c("ols", "str")){
    out <- rep(telist$comb, n)
  }else{
    res <- split(res, 1:NROW(res))

    out <- lapply(1:n, function(i){
      input_list[[i]]$res <- res[[i]]
      do.call("tecov", input_list[[i]])
    })
  }

  return(out)
}

# Cross-sectional reconciliation step
csstep <- function(base, cslist, csc_list, kset, ...){
  if(is.null(cslist$comb) || cslist$comb %in% c("ols", "str")){
    cslist$base <- t(base)
    out <- do.call("csrec", cslist)
    return(unname(t(out)))
  }else{
    base <- mat2list(base, kset)
    input_list <- NULL
    input_list[as.character(kset)] <- list(cslist)

    out <- sapply(as.character(kset), function(k){
      input_list[[k]]$base <- base[[k]]
      input_list[[k]]$comb <- csc_list[[k]]
      do.call("csrec", input_list[[k]])
    })
    return(unname(t(do.call("rbind", out))))
  }
}

csstepM <- function(base, cslist, res = NULL, kset, ...){
  base <- mat2list(base, kset)
  input_list <- NULL
  input_list[as.character(kset)] <- list(cslist)
  if(!is.null(res)){
    res <- mat2list(res, kset)
  }

  out <- sapply(as.character(kset), function(k){
    input_list[[k]]$base <- base[[k]]
    input_list[[k]]$res <- res[[k]]
    do.call("csrec", input_list[[k]])
  })
  unname(t(do.call("rbind", out)))
}

# Temporal reconciliation step
testep <- function(base, telist, tec_list, kset, ...){
  if(is.null(telist$comb) || telist$comb %in% c("ols", "str")){
    h <- NCOL(base)/sum(max(kset)/kset)
    n <- NROW(base)
    id <- rep(rep(rev(kset), h*max(kset)/kset), n)
    telist$base <- as.vector(t(base))[order(id)]
    out <- do.call("terec", telist)[order(order(id))]
    out <- matrix(out, ncol = n)
  }else{
    base <- split(base, 1:NROW(base))
    input_list <- NULL
    input_list[1:NROW(base)] <- list(telist)

    out <- sapply(1:NROW(base), function(i){
      input_list[[i]]$base <- base[[i]]
      input_list[[i]]$comb <- tec_list[[i]]
      do.call("terec", input_list[[i]])
    })
  }
  unname(t(out))
}

testepM <- function(base, telist, res = NULL, ...){
  base <- split(base, 1:NROW(base))
  input_list <- NULL
  input_list[1:NROW(base)] <- list(telist)
  if(!is.null(res)){
    res <- split(res, 1:NROW(res))
  }

  out <- sapply(1:NROW(base), function(i){
    input_list[[i]]$base <- base[[i]]
    input_list[[i]]$res <- res[[i]]
    do.call("terec", input_list[[i]])
  })
  unname(t(out))
}

# Check iterative convergence
check_results <- function(flag, dmat, i, itmax, norm, verbose){
  if(flag == 0){ # Convergence achieved
    if(verbose){
      cat("\n")
      cli_alert_success("Convergence achieved at iteration {i}.")
      #cat("-------------------------------------------------------------------\n")
      #cat("Convergence achieved at iteration", i, "\n")
    }
    check_sorted <- all(apply(dmat[,norm,,], c(1,3),
                              function(x) !is.unsorted(rev(na.omit(x[-1])))))
    if(!check_sorted){
      # incoherence has increased in the next iteration (at least one time)
      if(verbose){
        cli_alert_warning("Incoherence has increased in the next iteration (at least one time)")
      }
      flag <- 1
      return(flag)
    }else{
      # Perfect convergence
      return(flag)
    }
  }else{
    # Convergence not achieved (maximum iteration limit reached)
    if(verbose){
      cat("\n")
      cli_alert_danger("Convergence NOT achieved. Maximum number of iterations reached ({itmax}).")
      #cat("-------------------------------------------------------------------\n")
      #cat("Convergence NOT achieved: Maximum number of iterations reached (",
      #    itmax, ")\n", sep = "")
    }
    return(flag)
  }
}

# Discrepancy function
discrepancy <- function(mat, cs_const_mat, te_const_mat){
  dmat <- matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("cs", "te"),
                                                         c("inf", "one", "two")))
  dmat[,1] <- c(max(abs(cs_const_mat%*%mat)),
                max(abs(mat%*%t(te_const_mat))))
  dmat[,2] <- c(sum(abs(cs_const_mat%*%mat)),
                sum(abs(mat%*%t(te_const_mat))))
  dmat[,3] <- c(sum((cs_const_mat%*%mat)^2),
                sum((mat%*%t(te_const_mat))^2))
  dmat[dmat<sqrt(.Machine$double.eps)] <- 0
  dmat
}
