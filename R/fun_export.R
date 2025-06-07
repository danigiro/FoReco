#' @title Commutation matrix
#'
#' @description
#' This function returns the (\eqn{r c \times r c})
#' commutation matrix \eqn{\mathbf{P}} such that
#' \eqn{\mathbf{P} \mbox{vec}(\mathbf{Y}) = \mbox{vec}(\mathbf{Y}'),}
#' where \eqn{\mathbf{Y}} is a (\eqn{r \times c}) matrix (Magnus and Neudecker, 2019).
#'
#' @usage
#' # Commutation matrix
#' commat(r, c)
#'
#' @param r Number of rows of \eqn{\mathbf{Y}}.
#' @param c Number of columns of \eqn{\mathbf{Y}}.
#'
#' @returns A sparse (\eqn{r c \times r c}) matrix \eqn{\mathbf{P}} ([commat]),
#' or the vector of indexes for the rows of commutation matrix \eqn{\mathbf{P}}
#' ([commat_index])
#'
#' @references
#' Magnus, J.R. and Neudecker, H. (2019), Matrix Differential Calculus with Applications
#' in Statistics and Econometrics, third edition, New York, Wiley, pp. 54-55.
#'
#' @family Utilities
#'
#' @examples
#' Y <- matrix(rnorm(30), 5, 6)
#' P <- commat(5, 6)
#' P %*% as.vector(Y) == as.vector(t(Y)) # check
#'
#' @export
commat <- function(r, c){
  I <- seq(1:(r*c))[order(rep(1:r, c))]
  #I <- Matrix(1:(r * c), c, r, byrow = T) # initialize a matrix of indices of size (c x r)
  #I <- as.vector(I) # vectorize the required indices
  P <- Diagonal(r * c) # Initialize an identity matrix
  P <- P[I, ] # Re-arrange the rows of the identity matrix
  return(P)
}


#' @rdname commat
#'
#' @usage
#' # Vector of indexes for the rows of commutation matrix
#' commat_index(r, c)
#'
#' @export
commat_index <- function(r, c){
  return(seq(1:(r*c))[order(rep(1:r, c))])
}

#' @title Shrinkage of the covariance matrix
#'
#' @description
#' Shrinkage of the covariance matrix according to \enc{Schäfer}{Schafer} and Strimmer (2005).
#'
#' @param x A numeric matrix containing the in-sample residuals.
#' @param mse If \code{TRUE} (\emph{default}), the residuals used to compute the covariance
#' matrix are not mean-corrected.
#'
#' @returns A shrunk covariance matrix.
#'
#' @references
#' \enc{Schäfer}{Schafer}, J.L. and Strimmer, K. (2005), A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics,
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 4, 1
#'
#' @family Utilities
#'
#' @export
shrink_estim <- function(x, mse = TRUE){
  if(is.matrix(x) == TRUE && is.numeric(x) == FALSE){
    cli_abort("{.arg x} is not a numeric matrix.", call = NULL)
  }
  x <- remove_na(x)
  n <- nrow(x)

  # Full
  covm <- sample_estim(x = x, mse = mse)

  # Target
  tar <- Diagonal(x = diag(covm))

  x <- na.omit(x)
  if(NROW(x)>3){
    # Lambda
    xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
    xs[is.nan(xs)] <- xs[is.na(xs)] <- 0
    #xs <- xs[stats::complete.cases(xs), ]
    vS <- (1 / (n * (n - 1))) * (crossprod(xs^2) - ((1 / n) * (crossprod(xs))^2))
    diag(vS) <- 0
    corm <- covcor(covm)
    corm[is.nan(corm)] <- 0
    diag(corm) <- diag(corm)-1
    corm <- corm^2
    lambda <- sum(vS) / sum(corm)
    if(is.nan(lambda)){
      lambda <- 1
    }
    lambda <- max(min(lambda, 1), 0)
  }else{
    lambda <- 1
  }

  # Shrinkage
  shrink_cov <- lambda * tar + (1 - lambda) * covm
  shrink_cov <- drop0(shrink_cov)
  attr(shrink_cov, "lambda") <- lambda

  return(shrink_cov)
}

#' @title Shrinkage of the covariance matrix using the Oracle approximation
#'
#' @description
#' Shrinkage of the covariance matrix according to the Oracle Approximating Shrinkage (OAS)
#' of Chen et al. (2009) and Ando and Xiao (2023).
#'
#' @param x A numeric matrix containing the in-sample residuals.
#' @param mse If \code{TRUE} (\emph{default}), the residuals used to compute the covariance
#' matrix are not mean-corrected.
#'
#' @returns A shrunk covariance matrix.
#'
#' @references
#' Ando, S., and Xiao, M. (2023), High-dimensional covariance matrix estimation:
#' shrinkage toward a diagonal target. \emph{IMF Working Papers}, 2023(257), A001.
#'
#' Chen, Y., Wiesel, A., and Hero, A. O. (2009), Shrinkage estimation of high dimensional
#' covariance matrices, \emph{2009 IEEE international conference on acoustics, speech and
#' signal processing}, 2937–2940. IEEE.
#'
#' @family Utilities
#'
#' @export
shrink_oasd <- function(x, mse = TRUE){
  if(is.matrix(x) == TRUE && is.numeric(x) == FALSE){
    cli_abort("{.arg x} is not a numeric matrix.", call = NULL)
  }
  x <- remove_na(x)
  n <- nrow(x)

  # Full
  covm <- sample_estim(x = x, mse = mse)

  # Target
  tar <- diag(diag(covm))

  # lambda
  tr_ss <- sum(diag(covm%*%covm))
  var <- diag(covm)
  tr_s <- sum(var^2)
  num <- tr_ss - tr_s
  den <- tr_ss + sum(var)^2 - 2*tr_s
  phi <- num/den
  lambda <- min(1, 1/((n+1)*phi))

  # Shrinkage
  shrink_cov <- lambda * tar + (1 - lambda) * covm
  shrink_cov <- drop0(shrink_cov)
  attr(shrink_cov, "lambda") <- lambda

  return(shrink_cov)
}

#' Aggregation matrix of a (possibly) unbalanced hierarchy in balanced form
#'
#' A hierarchy with \eqn{L} upper levels is said to be balanced if each variable at level
#' \eqn{l} has at least one child at level \eqn{l+1}. When this doesn't hold, the hierarchy
#' is unbalanced. This function transforms an aggregation matrix of an unbalanced hierarchy
#' into an aggregation matrix of a balanced one. This function is used to reconcile forecasts
#' with [cslcc], which operates exclusively with balanced hierarchies.
#'
#' @inheritParams cslcc
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#'
#' @return A list containing four elements:
#' \item{bam}{The balanced aggregation matrix.}
#' \item{agg_mat}{The input matrix.}
#' \item{nodes}{A (\eqn{L \times 1}) numeric vector indicating the number of variables
#' in each of the \eqn{L} upper levels of the balanced hierarchy.}
#' \item{id}{The identification number of each variable in the balanced hierarchy.
#' It may contains duplicated values.}
#'
#' @family Utilities
#'
#' @export
#'
#' @examples
#' #    Unbalanced     ->      Balanced
#' #        T                     T
#' #    |-------|             |-------|
#' #    A       |             A       B
#' #  |---|     |           |---|     |
#' # AA   AB    B          AA   AB    BA
#' A <- matrix(c(1, 1, 1,
#'               1, 1, 0), 2, byrow = TRUE)
#' obj <- balance_hierarchy(agg_mat = A, nodes = c(1, 1))
#' obj$bam
balance_hierarchy <- function(agg_mat, nodes = "auto", sparse = TRUE){
  if(missing(agg_mat)){
    cli_abort("Argument {.arg agg_mat} is missing, with no default.", call = NULL)
  }

  tmp <- cstools(agg_mat = agg_mat)
  n <- tmp$dim[["n"]]
  na <- tmp$dim[["na"]]
  nb <- tmp$dim[["nb"]]
  agg_mat <- tmp$agg_mat

  # check hierarchy agg_mat (only 0 and 1)
  if(!all(unique(as.vector(agg_mat)) %in% c(1,0))){
    cli_abort("A hierarchical aggregation matrix has only 1s and 0s. Check {.arg agg_mat}.", call = NULL)
  }

  if(is.character(nodes[1])){
    out <- find_nodes(agg_mat)
    return(sparse2dense(out, sparse = sparse))
  }

  # Check nodes
  if(sum(nodes)!=na){
    cli_abort("The sum of {.arg nodes} must be equal to {na}, number of upper level time series.", call = NULL)
  }

  lev_na <- rep(1:length(nodes), nodes)
  Ident <- .sparseDiagonal(nb)

  tmp <- lapply(1:length(nodes),
                function(id){
                  idna <- which(lev_na == id)
                  lev_mat <- agg_mat[idna,,drop = FALSE]
                  if(!all(colSums(lev_mat) %in% c(1,0))){
                    cli_warn("In level {id} at least one bottom time
                                     series is present in two upper time series.
                                     Check {.arg agg_mat}.", call = NULL)
                  }
                  id_unb <- which(colSums(lev_mat)==0)
                  id_base <- c(idna, na+id_unb)
                  bal_mat <- rbind(lev_mat, Ident[id_unb,,drop = FALSE])
                  new_nodes <- NROW(bal_mat)
                  return(list(bam = bal_mat, id_base = id_base, new_nodes = new_nodes))
                })
  bam <- do.call(rbind, lapply(tmp, "[[", "bam"))
  id_base <- do.call(c, lapply(tmp, "[[", "id_base"))
  new_nodes <- do.call(c, lapply(tmp, "[[", "new_nodes"))

  out <- list(bam = bam,
              agg_mat = agg_mat,
              nodes = new_nodes,
              id = c(unname(id_base), (na+1):n))
  return(sparse2dense(out, sparse = sparse))
}

find_nodes <- function(agg_mat){

  strc_mat <- rbind(agg_mat, .sparseDiagonal(NCOL(agg_mat)))
  excl_list <- lapply(split(t(strc_mat), 1:NCOL(strc_mat)), function(x) which(x == 1))

  lop <- sapply(1:NROW(excl_list), function(id){
    tmp <- sort(unique(unlist(lapply(excl_list, function(lx) if(id %in% lx) lx))))
    tmp[tmp != id]
  })

  levels <- NULL
  k <- 1
  i <- 1
  while (i <= NROW(agg_mat)) {
    levid <- i
    tmp_i <- levid:NROW(strc_mat)
    if (sum(agg_mat[levid, ]) == NCOL(agg_mat)) {
      levels[[k]] <- levid
    }else {
      while (sum(strc_mat[levid, ]) != NCOL(agg_mat)) {
        del <- sort(unique(unlist(unique(lop[levid]))))
        check <- tmp_i[!(tmp_i %in% del) & !(tmp_i %in% levid)]
        levid <- c(levid, check[1])
      }
      levels[[k]] <- levid
    }
    tmp_i <- sort(setdiff(tmp_i, levels[[k]]))
    if (length(tmp_i) > 0) {
      i <- tmp_i[1]
    }else {
      break
    }
    k <- k + 1
  }

  return(list(bam = strc_mat[unlist(levels), , drop = FALSE],
              agg_mat = agg_mat,
              nodes = sapply(levels, length),
              id = c(unlist(levels), (NROW(agg_mat)+1):NROW(strc_mat))))
}



#' Aggregation matrix of a balanced hierarchy in (possibly) unbalanced form
#'
#' A hierarchy with \eqn{L} upper levels is said to be balanced if each variable at level
#' \eqn{l} has at least one child at level \eqn{l+1}. When this doesn't hold, the
#' hierarchy is unbalanced.
#' This function transforms an aggregation matrix of a balanced hierarchy
#' into an aggregation matrix of an unbalanced one, by removing possible duplicated series.
#'
#' @inheritParams balance_hierarchy
#' @param more_info If \code{TRUE}, it returns only the aggregation matrix
#' of the unbalanced hierarchy. \emph{Default} is \code{FALSE}.
#'
#' @returns A list containing four
#' elements (\code{more_info = TRUE}):
#' \item{ubm}{The aggregation matrix of the unbalanced hierarchy.}
#' \item{agg_mat}{The input matrix.}
#' \item{idrm}{The identification number of the duplicated variables (row numbers of
#' the aggregation matrix \code{agg_mat}).}
#' \item{id}{The identification number of each variable in the balanced hierarchy.
#' It may contains duplicated values.}
#'
#' @family Utilities
#'
#' @export
#'
#' @examples
#' #     Balanced     ->     Unbalanced
#' #        T                    T
#' #    |-------|            |-------|
#' #    A       B            A       |
#' #  |---|     |          |---|     |
#' # AA   AB    BA        AA   AB    BA
#' A <- matrix(c(1, 1, 1,
#'               1, 1, 0,
#'               0, 0, 1), 3, byrow = TRUE)
#' obj <- unbalance_hierarchy(agg_mat = A)
#' obj
unbalance_hierarchy <- function(agg_mat, more_info = FALSE, sparse = TRUE){
  agg_mat <- cstools(agg_mat = agg_mat)$agg_mat
  if(!more_info){
    ubm <- Matrix(unique(as.matrix(agg_mat)), sparse = TRUE)
    out <- ubm[-which(rowSums(abs(ubm)) == 1), , drop = FALSE]
  }else{
    idnb <- which(rowSums(abs(agg_mat)) == 1)
    idna <- anyDuplicated(as.matrix(agg_mat))
    idna <- idna[idna!=0]
    idrm <- sort(unique(c(idna,idnb)))
    ubm <- agg_mat[-idrm, ,drop = FALSE]

    iddu <- apply(agg_mat[idrm,, drop = FALSE], 1,
                  function(x){
                    if(any(sum(abs(x)) == 1)){
                      NROW(ubm) + which(x == 1)
                    }else{
                      which(apply(ubm, 1, function(y)
                        all(x == y)))[1]
                    }
                  })
    id <- 1:NROW(agg_mat)
    id[id %in% idrm] <- iddu
    out <- list(ubm = ubm,
                agg_mat = agg_mat,
                idrm = idrm,
                id = c(id, (NROW(ubm)+1):sum(dim(ubm))))
  }
  return(sparse2dense(out, sparse = sparse))
}

#' Informations on the reconciliation process
#'
#' @description
#' This function extracts reconciliation information from the output of any reconciled
#' function implemented by \pkg{FoReco}.
#'
#' @param x An output from any reconciliation function implemented by \pkg{FoReco}.
#' @param verbose If \code{TRUE} (\emph{defaults}), reconciliation information are printed.
#'
#' @returns A list containing the following reconciliation process informations:
#'   \item{rfun}{the reconciliation function.}
#'   \item{cs_n}{the cross-sectional number of variables.}
#'   \item{te_set}{the set of temporal aggregation orders.}
#'   \item{forecast_horizon}{the forecast horizon
#'   (in temporal and cross-temporal frameworks, for the most temporally aggregated series).}
#'   \item{framework}{the reconciliation framework (cross-sectional, temporal or cross-temporal).}
#'   \item{info}{non-negative reconciled forecast convergence information.}
#'   \item{lcc}{list of level conditional reconciled forecasts (+ BU) for
#'   [cslcc], [telcc] and [ctlcc].}
#'   \item{nn}{if \code{TRUE}, all the forecasts are not negative.}
#'   \item{comb}{the covariance approximation.}
#'
#' @family Utilities
#' @export
#'
recoinfo <- function(x, verbose = TRUE){
  if(is.null(attr(x,"FoReco"))){
    cli_warn(c("!"="No information available."), call = NULL)
    invisible(NULL)
  }else{
    out <- as.list(attr(x,"FoReco"))
    out$nn <- all(x>=0)
    if(verbose){
      if(out$rfun %in% c("cslcc", "telcc", "ctlcc")){
        title <- "Level Conditional Coherent "
      }else if(out$rfun %in% c("csrec", "terec", "ctrec")){
        title <- "Optimal "
      }else if(out$rfun %in% c("iterec", "tcsrec", "cstrec")){
        title <- "Heuristic "
      }else if(out$rfun %in% c("ctbu", "csbu", "tebu")){
        title <- "Bottom-up "
      }else if(out$rfun %in% c("cttd", "cstd", "tetd")){
        title <- "Top-down "
      }else if(out$rfun %in% c("ctmo", "csmo", "temo")){
        title <- "Middle-out "
      }else{
        title <- " "
      }
      cli_alert_success("{.emph {title}}{.strong {out$framework}} Forecast Reconciliation")
      if(!is.null(out$rfun)) cli_alert_info("{.pkg FoReco} function: {.strong {.code {out$rfun}}}")

      if(!is.null(out$comb)){
        if(length(out$comb)>1){
          tmp <- paste(names(out$comb), out$comb, sep = "-")
        }else{
          tmp <- out$comb
        }
        cli_alert_info("Covariance approximation: {.strong {.code {tmp}}}")
      }
      if(!is.null(out$nn)) cli_alert_info("Non-negative forecasts: {.strong {.code {out$nn}}}")
    }

    invisible(out)
  }
}

#' Reconciled forecasts to matrix/vector
#'
#' @description
#' This function splits the temporal vectors and the cross-temporal matrices in a
#' list according to the temporal aggregation order
#'
#' @param x An output from any reconciliation function implemented by \pkg{FoReco}.
#' @inheritParams ctrec
#' @param keep_names If \code{FALSE} (\emph{default}), the rownames names of the output matrices are removed.
#'
#' @returns A list of matrices or vectors distinct by temporal aggregation order.
#'
#' @family Utilities
#' @export
#' @examples
#' set.seed(123)
#' # (3 x 7) base forecasts matrix (simulated), Z = X + Y and m = 4
#' base <- rbind(rnorm(7, rep(c(20, 10, 5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))),
#'               rnorm(7, rep(c(10, 5, 2.5), c(1, 2, 4))))
#'
#' reco <- ctrec(base = base, agg_mat = t(c(1,1)), agg_order = 4, comb = "ols")
#' matrix_list <- FoReco2matrix(reco)
FoReco2matrix <- function(x, agg_order, keep_names = FALSE){
  if(!is.null(attr(x, "FoReco"))){
    fr <- recoinfo(x, verbose = FALSE)
    frame <- fr$framework
    set <- fr$te_set
    h <- fr$forecast_horizon
  }else if(!missing(agg_order)){
    frame <- ifelse(NCOL(x)==1, "Temporal", "Cross-temporal")
    set <- tetools(agg_order = agg_order)$set
    h <- ifelse(NCOL(x)==1, length(x), NCOL(x))/sum(max(set)/set)
  }else{
    frame <- "Cross-sectional"
  }

  if(frame == "Cross-sectional"){
    attr(x, "FoReco") <- NULL
    return(list("k-1"= x))
  }else{
    id <- rep(set, h*max(set)/set)

    if(NCOL(x)==1){
      out <- split(x, factor(id, set))
      if(!keep_names){
        out <- lapply(out, unname)
      }
    }else{
      out <- lapply(setNames(set, set), function(k){
        mat <- t(x[, id == k, drop = FALSE])
        if(!keep_names)
          rownames(mat) <- NULL
        mat
      })
    }

    names(out) <- paste0("k-", names(out))
    return(out)
  }
}

#' Cross-sectional aggregation matrix of a dataframe
#'
#' @description
#' This function allows the user to easily build the (\eqn{n_a \times n_b})
#' cross-sectional aggregation matrix starting from a data frame.
#'
#' @usage
#' df2aggmat(formula, data, sep = "_", sparse = TRUE, top_label = "Total",
#'           verbose = TRUE)
#'
#' @param formula  Specification of the hierarchical structure: grouped hierarchies are specified
#' using \code{~ g1 * g2} and nested hierarchies are specified using \code{~ parent / child}.
#' Mixtures of the two formulations are also possible, like \code{~ g1 * (grandparent / parent / child)}.
#' @param data A dataset in which each column contains the values of the variables in the formula
#' and each row identifies a bottom level time series.
#' @param sep Character to separate the names of the aggregated series, (\emph{default} is "\code{_}").
#' @param sparse Option to return sparse matrices (\emph{default} is \code{TRUE}).
#' @param top_label Label of the top level variable (\emph{default} is "\code{Total}").
#' @param verbose If \code{TRUE} (\emph{default}), hierarchy informations are printed.
#'
#' @return A (\code{na x nb}) matrix.
#'
#' @family Utilities
#'
#' @examples
#' ## Balanced hierarchy
#' #         T
#' #    |--------|
#' #    A        B
#' #  |---|   |--|--|
#' # AA   AB  BA BB BC
#' # Names of the bottom level variables
#' data_bts <- data.frame(X1 = c("A", "A", "B", "B", "B"),
#'                        X2 = c("A", "B", "A", "B", "C"),
#'                        stringsAsFactors = FALSE)
#' # Cross-sectional aggregation matrix
#' agg_mat <- df2aggmat(~ X1 / X2, data_bts, sep = "", verbose = FALSE)
#'
#' ## Unbalanced hierarchy
#' #                 T
#' #       |---------|---------|
#' #       A         B         C
#' #     |---|     |---|     |---|
#' #    AA   AB   BA   BB   CA   CB
#' #  |----|         |----|
#' # AAA  AAB       BBA  BBB
#' # Names of the bottom level variables
#' data_bts <- data.frame(X1 = c("A", "A", "A", "B", "B", "B", "C", "C"),
#'                        X2 = c("A", "A", "B", "A", "B", "B", "A", "B"),
#'                        X3 = c("A", "B", NA, NA, "A", "B", NA, NA),
#'                        stringsAsFactors = FALSE)
#' # Cross-sectional aggregation matrix
#' agg_mat <- df2aggmat(~ X1 / X2 / X3, data_bts, sep = "", verbose = FALSE)
#'
#' ## Group of two hierarchies
#' #     T          T         T | A  | B  | C
#' #  |--|--|  X  |---|  ->  ---+----+----+----
#' #  A  B  C     M   F       M | AM | BM | CM
#' #                          F | AF | BF | CF
#' # Names of the bottom level variables
#' data_bts <- data.frame(X1 = c("A", "A", "B", "B", "C", "C"),
#'                        Y1 = c("M", "F", "M", "F", "M", "F"),
#'                        stringsAsFactors = FALSE)
#' # Cross-sectional aggregation matrix
#' agg_mat <- df2aggmat(~ Y1 * X1, data_bts, sep = "", verbose = FALSE)
#'
#' @export
df2aggmat <- function(formula, data, sep = "_", sparse = TRUE,
                      top_label = "Total", verbose = TRUE) {
  if (missing(data)) {
    cli_abort("Argument {.arg data} is missing, with no default.", call = NULL)
  }

  if (NCOL(data) == 1) {
    out <- matrix(1, 1, NROW(unique(data)))
    rownames(out) <- top_label
    colnames(out) <- unique(data[, 1, , drop = TRUE])
    return(out)
  }

  if (missing(formula)) {
    formula <- as.formula(paste("~", paste(colnames(data), collapse = "*"), sep = ""))
    cli_alert_warning("Argument {.arg data} is missing, default is {.arg formula = ~ {formula[2]}}")
  }

  tm <- terms(formula)
  lev <- attr(tm, "factors")
  lev_vars <- rownames(lev)
  lev <- Map(function(x) lev_vars[x != 0], split(lev, col(lev)))
  lev <- unname(lev)
  lev <- lev[lapply(lev, length)<length(lev_vars)]

  if (!all(lev_vars %in% colnames(data))) {
    cli_abort("Colnames of {.arg data} don't contain all the variables in formula.", call = NULL)
  }
  data <- data[, which(colnames(data) %in% lev_vars)]

  out <- lapply(lev, function(x) {
    id <- which(!colnames(data) %in% x)
    datax <- data
    datax[, id] <- NA
    datax <- unique(datax)
    data_id <- data[, -id, drop = FALSE]
    datax <- cbind(datax, t(apply(
      datax[, -id, drop = FALSE], 1,
      function(x) {
        data_s <- data.frame(matrix(x,
                                    nrow = NROW(data_id),
                                    ncol = length(x), byrow = TRUE
        ))
        as.numeric(apply(data_id == data_s, 1, all))
      }
    )))
    return(datax)
  })
  out <- do.call("rbind", out)
  namerows <- apply(out[, 1:NCOL(data)], 1, function(x) paste(stats::na.omit(x), collapse = sep))
  namecols <- apply(data, 1, function(x) paste(stats::na.omit(x), collapse = sep))
  out <- as.matrix(out[, -c(1:NCOL(data)), drop = FALSE])
  out <- sparse2dense(out, sparse = sparse)
  rownames(out) <- namerows
  colnames(out) <- namecols

  out <- out[rowSums(out) > 1 & rowSums(out) < NCOL(out), , drop = FALSE]
  out <- rbind(Total = 1, out)
  rownames(out)[1] <- top_label
  if(verbose){
    cli_rule("{.strong Cross-sectional information}")
    cli_bullets(c("# levels: {length(lev)+1}",
                  "{.emph # total time series} ({.strong n}): {NROW(out) + NCOL(out)}",
                  "{.emph # upper time series} ({.strong na}): {NROW(out)}",
                  "{.emph # bottom time series} ({.strong nb}): {NCOL(out)}"))
  }
  return(out)
}



#' Non-overlapping temporal aggregation of a time series
#'
#' @description
#' Non-overlapping temporal aggregation of a time series according
#' to a specific aggregation order.
#'
#' @param agg_order A numeric vector with the aggregation orders to consider.
#' @param y Univariate or multivariate time series: a vector/matrix or a \code{ts} object.
#' @param align A string or a vector specifying the alignment of \code{y}. Options include:
#' "\code{end}" (end of the series, \emph{default}), "\code{start}" (start of the series),
#' an integer (or a vector of integers) indicating the starting period of the temporally
#' aggregated series.
#' @param rm_na If \code{TRUE} the missing values are removed.
#' @inheritParams tetools
#'
#' @return A list of vectors or \code{ts} objects.
#'
#' @family Utilities
#'
#' @export
#'
#' @examples
#' # Monthly time series (input vector)
#' y <- ts(rnorm(24), start = 2020, frequency = 12)
#' # Quarterly time series
#' x1 <- aggts(y, 3)
#' # Monthly, quarterly and annual time series
#' x2 <- aggts(y, c(1, 3, 12))
#' # All temporally aggregated time series
#' x3 <- aggts(y)
#'
#' # Ragged data
#' y2 <- ts(rnorm(11), start = c(2020, 3), frequency = 4)
#' # Annual time series: start in 2021
#' x4 <- aggts(y2, 4, align = 3)
#' # Semi-annual (start in 2nd semester of 2020) and annual (start in 2021) time series
#' x5 <- aggts(y2, c(2, 4), align = c(1, 3))
#'
aggts <- function(y, agg_order, tew = "sum", align = "end", rm_na = FALSE){

  if(is.ts(y)){
    tspy <- tsp(y)
    y <- as.matrix(y)
    kset <- rev(all_factors(tspy[3]))
  }else{
    if(is.vector(y)){
      y <- cbind(y)
    }
    kset <- tspy <- NULL
  }

  n <- NROW(y)

  if(missing(agg_order)){
    if(!is.null(kset)){
      agg_order <- kset
    }else{
      cli_abort("Argument {.arg agg_order} is missing, with no default.", call = NULL)
    }
  }

  if(is.character(align)){
    align <- match.arg(align, c("end","start"))
    if(align=='end'){
      start <- n%%agg_order + 1L
    }else if(align=='start'){
      start <- rep(1L, length(agg_order))
    }
  }else{
    if(length(align)==1){
      start <- rep(align, length(agg_order))
    }else if(length(align)==length(agg_order)){
      start <- align
    }else{
      cli_abort("Argument {.arg align} can be a character ('end' or 'start'), a
                number or a vector with the same lenght of {.arg agg_order}", call = NULL)
    }
  }

  start <- setNames(start, as.character(agg_order))
  agg_order <- sort(as.integer(agg_order))
  nk <- setNames(trunc(n/agg_order), as.character(agg_order))
  out <- lapply(agg_order, function(k){
    out <- apply(y, 2, function(col){
      tmp <- matrix(col[start[as.character(k)] - 1L + seq_len(k*nk[as.character(k)])],
                    ncol=nk[as.character(k)])
      if(tew == "sum"){
        tmp <- colSums(tmp)
      }else if(tew == "avg"){
        tmp <- colMeans(tmp)
      }else if(tew == "last"){
        tmp <- tmp[NROW(tmp),]
      }else if(tew == "first"){
        tmp <- tmp[1,]
      }else{
        cli_abort("{.code {tew}} is not implemented yet.", call = NULL)
      }

      if(is.character(align)){
        if(align=='end' & n%%k != 0){
          tmp <- c(NA, tmp)
        }else if(align=='start' & n%%k != 0){
          tmp <- c(tmp, NA)
        }
      }else{
        if(start[as.character(k)] != 1 ){
          tmp <- c(NA, tmp)
        }

        if((n-start[as.character(k)]+1)%%k != 0){
          tmp <- c(tmp, NA)
        }
      }
      return(tmp)
    }, simplify = FALSE)
    out <- do.call(cbind, out)

    if(NCOL(out)==1){
      out <- out[,]
    }

    if(!is.null(tspy)){
      out <- ts(out, frequency=tspy[3]/k,
                start=floor(tspy[1]))
    }

    if(rm_na){
      out <- na.omit(out)
    }
    out
  })

  if(length(out) == 1){
    return(out[[1]])
  }else{
    names(out) <- paste0("k-", agg_order)
    return(out)
  }
}

#' Set bounds for bounded forecast reconciliation
#'
#' This function defines the bounds matrix considering cross-sectional,
#' temporal, or cross-temporal frameworks. The output matrix can be used as
#' input for the \code{bounds} parameter in functions such as [csrec], [terec],
#' or [ctrec], to perform bounded reconciliations.
#'
#' @usage set_bounds(n, k, h, lb = -Inf, ub = Inf, approach = "osqp", bounds = NULL)
#'
#' @param n A (\eqn{b \times 1}) vector representing the \eqn{i}th cross-sectional
#' series (\eqn{i = 1, \dots, n}), where \eqn{b} is the number of bounds to be set.
#' @param k A (\eqn{b \times 1}) vector specifying the temporal aggregation orders
#' (\eqn{k = m, \dots, 1}).
#' @param h A (\eqn{b \times 1}) vector representing the forecast horizons
#' (\eqn{j = 1, \dots, m/k}).
#' @param lb,ub A (\eqn{b \times 1}) vector of lower and upper bounds.
#' @param approach A string specifying the algorithm to compute bounded reconciled forecasts:
#'   \itemize{
#'   \item "\code{osqp}": quadratic programming optimization
#'   (\href{https://osqp.org/}{\pkg{osqp}} solver).
#'   \item "\code{sftb}": heuristic "set-forecasts-to-bounds", which adjusts the reconciled
#'   forecasts to be within specified bounds without further optimization.
#'   }
#' @param bounds A matrix of previous bounds to be added. If not specified,
#' new bounds will be computed.
#'
#' @return A numeric matrix representing the computed bounds, which can be:
#' \itemize{
#'   \item Cross-sectional (\eqn{b \times 3}) matrix for cross-sectional reconciliation ([csrec]).
#'   \item Temporal (\eqn{b \times 4}) matrix for temporal reconciliation ([terec]).
#'   \item Cross-temporal (\eqn{b \times 5}) matrix for cross-temporal reconciliation ([ctrec]).
#' }
#'
#' @family Utilities
#'
#' @examples
#' # Example 1
#' # Two cross-sectional series (i = 2,3),
#' # with each series required to be between 0 and 1.
#' n <- c(2, 3)
#' lb <- c(0, 0)
#' ub <- c(1,1)
#' bounds_mat <- set_bounds(n = c(2, 3),
#'                          lb = rep(0, 2), # or lb = 0
#'                          ub = rep(1, 2)) # or ub = 1
#'
#' # Example 2
#' # All the monthly values are between 0 and 1.
#' bounds_mat <- set_bounds(k = rep(1, 12),  # or k = 1
#'                          h = 1:12,
#'                          lb = rep(0, 12), # or lb = 0
#'                          ub = rep(1, 12)) # or ub = 1
#'
#' # Example 3
#' # For two cross-sectional series (i = 2,3),
#' # all the monthly values are between 0 and 1.
#' bounds_mat <- set_bounds(n = rep(c(2, 3), each = 12),
#'                          k = 1,
#'                          h = rep(1:12, 2),
#'                          lb = 0, # or lb = 0
#'                          ub = 1) # or ub = 1
#'
#' @export
set_bounds <- function(n, k, h, lb = -Inf, ub = Inf, approach = "osqp", bounds = NULL){
  bounds_old <- bounds

  if(!missing(n) & !missing(k)){
    if(missing(h)){
      cli_abort("Argument {.arg h} is missing, with no default.", call = NULL)
    }

    max_size <- max(length(k), length(h), length(lb), length(ub), length(n))

    if(length(n) != 1 & length(n) != max_size){
      cli_abort("{.arg n} length must be either 1 or {max_size}")
    }

    if(length(k) != 1 & length(k) != max_size){
      cli_abort("{.arg k} length must be either 1 or {max_size}")
    }

    if(length(h) != 1 & length(h) != max_size){
      cli_abort("{.arg h} length must be either 1 or {max_size}")
    }

    if(length(lb) != 1 & length(lb) != max_size){
      cli_abort("{.arg lb} length must be either 1 or {max_size}")
    }

    if(length(ub) != 1 & length(ub) != max_size){
      cli_abort("{.arg up} length must be either 1 or {max_size}")
    }

    bounds <- cbind(n, k, h, lb, ub)

  }else if(!missing(n)){

    max_size = max(length(lb), length(ub), length(n))

    if(length(n) != 1 & length(n) != max_size){
      cli_abort("{.arg n} length must be either 1 or {max_size}")
    }

    if(length(lb) != 1 & length(lb) != max_size){
      cli_abort("{.arg lb} length must be either 1 or {max_size}")
    }

    if(length(ub) != 1 & length(ub) != max_size){
      cli_abort("{.arg up} length must be either 1 or {max_size}")
    }

    bounds <- cbind(n, lb, ub)

  }else if(!missing(k)){
    max_size = max(length(k), length(h), length(lb), length(ub))
    if(length(k) != 1 & length(k) != max_size){
      cli_abort("{.arg k} length must be either 1 or {max_size}")
    }

    if(length(h) != 1 & length(h) != max_size){
      cli_abort("{.arg h} length must be either 1 or {max_size}")
    }

    if(length(lb) != 1 & length(lb) != max_size){
      cli_abort("{.arg lb} length must be either 1 or {max_size}")
    }

    if(length(ub) != 1 & length(ub) != max_size){
      cli_abort("{.arg up} length must be either 1 or {max_size}")
    }

    bounds <- cbind(k, h, lb, ub)
  }else{
    cli_abort("No arguments provide.")
  }

  if(!is.null(bounds_old)){
    bounds <- unique(rbind(bounds, bounds_old))
    app_old <- attr(bounds_old, "approach")

    if(!is.null(app_old)){
      attr(bounds, "approach") <- app_old
    }else{
      attr(bounds, "approach") <- approach
    }
  }else{
    attr(bounds, "approach") <- approach
  }
  return(bounds)
}
