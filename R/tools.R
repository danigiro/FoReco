#' Cross-sectional reconciliation tools
#'
#' @description
#' Some useful tools for the cross-sectional forecast reconciliation of a
#' linearly constrained (e.g., hierarchical/grouped) multiple time series.
#'
#' @inheritParams csrec
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @returns A list with four elements:
#' \item{dim}{A vector containing information about the number of series for the complete
#' system (\code{n}), for upper levels (\code{na}) and bottom level (\code{nb}).}
#' \item{agg_mat}{The cross-sectional aggregation matrix.}
#' \item{strc_mat}{The cross-sectional structural matrix.}
#' \item{cons_mat}{The cross-sectional zero constraints matrix.}
#'
#' @examples
#' # Cross-sectional framework
#' # One level hierarchy A = [1 1]
#' A <- matrix(1, 1, 2)
#' obj <- cstools(agg_mat = A)
#'
#' @family Framework: cross-sectional
#' @family Utilities
#' @export
cstools <- function(agg_mat, cons_mat, sparse = TRUE){

  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  }else if(!missing(agg_mat)){
    if(!(is.matrix(agg_mat) | is(agg_mat, "Matrix"))){
      cli_abort("{.arg agg_mat} is not a matrix.", call = NULL)
    }

  }else{
    if(!(is.matrix(cons_mat) | is(cons_mat, "Matrix"))){
      cli_abort("{.arg cons_mat} is not a matrix.", call = NULL)
    }

    n <- NCOL(cons_mat)
    na <- NROW(cons_mat)
    flag <- all(cons_mat[, 1:na, drop = FALSE] == .sparseDiagonal(na))
    if(flag){
      agg_mat <- -cons_mat[, -c(1:na), drop = FALSE]
      if(!is(agg_mat, "Matrix")){
        agg_mat <- Matrix(agg_mat, sparse = TRUE)
      }
    }else{
      if(sparse){
        return(list(dim = c("n" = n),
                    cons_mat = Matrix(cons_mat, sparse = TRUE)))
      }else{
        return(list(dim = c("n" = n),
                    cons_mat = as.matrix(cons_mat)))
      }
    }
  }

  if(!is(agg_mat, "Matrix")){
    agg_mat <- Matrix(agg_mat, sparse = TRUE)
  }

  if(any(rowSums(abs(agg_mat)) == 0)){
    cli_alert_info("Removed 0s rows in {.arg agg_mat}.")
    agg_mat <- agg_mat[rowSums(abs(agg_mat)) != 0,]
  }

  nb <- NCOL(agg_mat)
  na <- NROW(agg_mat)
  n <- na + nb
  strc_mat <- rbind(agg_mat, .sparseDiagonal(nb))
  if(!is.null(rownames(agg_mat)) & !is.null(colnames(agg_mat))){
    rownames(strc_mat) <- c(rownames(agg_mat), colnames(agg_mat))
  }
  cons_mat <- cbind(.sparseDiagonal(na), - agg_mat)
  if(!is.null(rownames(agg_mat)) & !is.null(colnames(agg_mat))){
    colnames(cons_mat) <- c(rownames(agg_mat), colnames(agg_mat))
  }

  return(sparse2dense(list(dim = c("n" = n, "na" = na, "nb" = nb),
                           agg_mat = agg_mat,
                           strc_mat = strc_mat,
                           cons_mat = cons_mat), sparse = sparse))
}

#' Temporal reconciliation tools
#'
#' @description
#' Some useful tools for forecast reconciliation through temporal hierarchies.
#'
#' @inheritParams terec
#' @param fh Forecast horizon for the lowest frequency (most temporally aggregated)
#' time series (\emph{default} is \code{1}).
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @returns A list with five elements:
#' \item{dim}{A vector containing information about the maximum aggregation order
#' (\code{m}), the number of factor (\code{p}), the partial (\code{ks}) and total
#' sum (\code{kt}) of factors.}
#' \item{set}{The vector of the temporal aggregation orders (in decreasing order).}
#' \item{agg_mat}{The temporal linear combination or aggregation matrix.}
#' \item{strc_mat}{The temporal structural matrix.}
#' \item{cons_mat}{The temporal zero constraints matrix.}
#'
#' @examples
#' # Temporal framework (quarterly data)
#' obj <- tetools(agg_order = 4, sparse = FALSE)
#'
#' @family Framework: temporal
#' @family Utilities
#' @export
tetools <- function(agg_order, fh = 1, tew = "sum", sparse = TRUE){

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else if(length(agg_order)>1){
    kset <- sort(agg_order, decreasing = TRUE)
  }else{
    kset <- rev(all_factors(agg_order))
  }
  m <- max(kset)

  if(m < 2){
    cli_abort("Argument {.arg agg_order} is not > 1.", call = NULL)
  }

  if(min(kset)!=1){
    kset <- sort(c(kset, 1), decreasing = TRUE)
  }

  if(!all(kset %in% rev(all_factors(m)))){
    cli_abort("({paste(kset, collapse = ', ')}) is not a subset
              of ({paste(rev(all_factors(m)), collapse = ', ')}).", call = NULL)
  }

  p <- length(kset)
  ks <- sum(m/kset[-p])
  kt <- sum(m/kset)

  if(is.null(tew)){
    tew <- "sum"
  }
  # weights
  if(is.character(tew[[1]][1])){
    if(tew == "sum"){
      weights <- lapply(kset[-p], function(x) rep(1, x))
    }else if(tew == "avg"){
      weights <- lapply(kset[-p], function(x) rep(1/x, x))
    }else if(tew == "last"){
      weights <- lapply(kset[-p], function(x) rep(c(0,1), c(x-1, 1)))
    }else if(tew == "first"){
      weights <- lapply(kset[-p], function(x) rep(c(1,0), c(1, x-1)))
    }else{
      cli_abort("{.code {tew}} tew is not implemented yet.", call = NULL)
    }
  }else if(is.list(tew)){
    if(length(tew) != (p-1)){
      cli_abort("{.arg weights} is not a list with {p-1} elements.", call = NULL)
    }

    if(any(sapply(tew, length) != kset[-p])){
      cli_abort("The {p-1} elements of {.arg weights} must be a vector of length {kset[-p]}, respectively.", call = NULL)
    }
    weights <- tew
  }

  agg_mat <- do.call("rbind", Map(
    function(x, y) kronecker(x, y), lapply(m/kset[-p], function(x) .sparseDiagonal(x * fh)),
    lapply(weights, rbind)
  ))
  cons_mat <- cbind(.sparseDiagonal(fh * ks), - agg_mat)
  strc_mat <- rbind(agg_mat, .sparseDiagonal(m * fh))

  return(sparse2dense(list(dim = c("m" = m,
                                    "p" = p,
                                    "ks" = ks,
                                    "kt" = kt),
                            set = kset,
                            agg_mat = agg_mat,
                            strc_mat = strc_mat,
                            cons_mat = cons_mat), sparse = sparse))
}

#' Cross-temporal reconciliation tools
#'
#' @description
#' Some useful tools for the cross-temporal forecast reconciliation of a linearly constrained
#' (e.g., hierarchical/grouped) multiple time series.
#'
#' @inheritParams ctrec
#' @param fh Forecast horizon for the lowest frequency (most temporally aggregated)
#' time series (\emph{default} is \code{1}).
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @returns A list with four elements:
#' \item{dim}{A vector containing information about the number of series for the
#' complete system (\code{n}), for upper levels (\code{na}) and bottom level
#' (\code{nb}), the maximum aggregation order  (\code{m}), the number of factor
#' (\code{p}), the partial (\code{ks}) and total sum (\code{kt}) of factors.}
#' \item{set}{The vector of the temporal aggregation orders (in decreasing order).}
#' \item{agg_mat}{The cross-temporal aggregation matrix.}
#' \item{strc_mat}{The cross-temporal structural matrix.}
#' \item{cons_mat}{The cross-temporal zero constraints matrix.}
#'
#' @examples
#' # Cross-temporal framework
#' A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
#' m <- 4 # from quarterly to annual temporal aggregation
#' cttools(agg_mat = A, agg_order = m)
#'
#' @family Framework: cross-temporal
#' @family Utilities
#' @export
cttools <- function(agg_mat, cons_mat, agg_order, tew = "sum", fh = 1, sparse = TRUE){

  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  }else if(!missing(agg_mat)){
    csobj <- cstools(agg_mat = agg_mat)
  }else{
    csobj <- cstools(cons_mat = cons_mat)
  }

  Ccs <- csobj$cons_mat
  Acs <- csobj$agg_mat
  Scs <- csobj$strc_mat
  n <- csobj$dim[["n"]]

  if(missing(agg_order)){
    cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
  }else{
    teobj <- tetools(agg_order = agg_order, fh = fh, tew = tew)
  }
  Cte <- teobj$cons_mat
  Ate <- teobj$agg_mat
  Ste <- teobj$strc_mat
  m <- teobj$dim[["m"]]
  ks <- teobj$dim[["ks"]]
  kt <- teobj$dim[["kt"]]

  #Cbr <- rbind(kronecker(Ccs, .sparseDiagonal(h * kt)), kronecker(.sparseDiagonal(n), Cte))
  P <- commat(fh * kt, n)

  Cup <- cbind(Matrix(0, fh * NROW(Ccs) * m, n * fh * ks),
               kronecker(.sparseDiagonal(fh * m), Ccs)) %*% P
  cons_mat <- rbind(Cup, kronecker(.sparseDiagonal(n), Cte))

  if(!is.null(Acs)){
    agg_mat <- rbind(kronecker(Acs, Ste),kronecker(.sparseDiagonal(NCOL(Acs)), Ate))
    strc_mat <- kronecker(Scs, Ste)
    return(sparse2dense(list(dim = c(csobj$dim, teobj$dim),
                             set = teobj$set,
                             agg_mat = agg_mat,
                             strc_mat = strc_mat,
                             cons_mat = cons_mat), sparse = sparse))
  }else{
    return(sparse2dense(list(dim = c(csobj$dim, teobj$dim),
                             set = teobj$set,
                             cons_mat = cons_mat), sparse = sparse))
  }
}
