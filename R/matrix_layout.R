#' Convert between horizon-stacked and cross-temporal layouts
#'
#' These functions convert matrix between the two canonical layouts used in
#' cross-temporal reconciliation.
#' Let \eqn{m} be the maximum temporal aggregation order and \eqn{k^\ast} the sum
#' of a subset of the \eqn{(p-1)} proper factors of \eqn{m} (excluding \eqn{m});
#' let \eqn{h} be the forecast horizon for the lowest frequency series (e.g.,
#' most aggregated temporal forecast horizon) and \eqn{n} the number of variables:
#' \itemize{
#'   \item \emph{Horizon-stacked layout (cross-temporal version)}: a
#'   \eqn{h \times n(k^\ast + m)} matrix where rows are the most aggregated temporal
#'   forecast horizons, and the values in each row are ordered from the lowest frequency
#'   (most temporally aggregated) to the highest frequency grouped by variable.
#'   \item \emph{Cross-temporal layout}: a \eqn{n \times h(k^\ast + m)} matrix where
#'   rows are variables, and horizons for each temporal block appear consecutively.
#'   rows are variables, and the values in each row are ordered from the lowest frequency
#'   (most temporally aggregated) to the highest frequency.
#' }
#' Then, [as_ctmatrix] converts a \eqn{(h \times n(k^\ast+m))}
#' horizon-stacked to a \eqn{(n \times h(k^\ast+m))} cross-temporal matrix;
#' [as_hstack_ctlayout] performs the inverse transform.
#'
#' @param hmat A \eqn{h \times n(k^\ast+m)} numeric matrix in \emph{horizon-stacked} layout
#' (cross-temporal version).
#' @param ctmat A \eqn{n \times h(k^\ast+m)} numeric matrix in \emph{cross-temporal} layout.
#' @inheritParams ctcov
#' @param row_names Optional character vector of length \code{n} with row names
#'   for the \emph{cross-temporal} output of \code{as_ctmatrix()}.
#'   If \code{NULL} (\emph{default}) no custom names are assigned.
#'
#' @return [as_ctmatrix] returns a \eqn{n \times h(k^\ast+m)} numeric
#' matrix in \emph{cross-temporal} layout.
#'
#' @examples
#' h <- 2   # horizons
#' n <- 3   # variables
#' m <- 4   # temporal aggregation order
#' kt <- tetools(m)$dim["kt"]
#'
#' # Build a horizon-stacked matrix: h rows, n * k_t columns
#' input_ct <- matrix(seq_len(h * n * kt), nrow = n, byrow = TRUE)
#'
#' hmat <- as_hstack_ctlayout(input_ct, agg_order = m)
#' ctmat <- as_ctmatrix(hmat, agg_order = m, n = n)
#' # all.equal(ctmat, input_ct, check.attributes = FALSE)
#'
#' @rdname ctmatrix_layouts
#' @export
as_ctmatrix <- function(hmat, agg_order, n, row_names = NULL){
  tmp <- tetools(agg_order = agg_order)

  if(is.vector(hmat)){
    hmat <- rbind(hmat)
  }

  if(NCOL(hmat) != tmp$dim["kt"]*n){
    cli_abort(c("x" = "Incorrect {.arg hmat} columns dimension.",
                "i" = "Expected {.val {tmp$dim['kt'] * n}} columns, but found {.val {NCOL(hmat)}}.",
                "i" = "Please check the consistency of {.arg n}, {.arg agg_order}, and {.arg hmat}."),
              call = NULL)
  }

  out <- hmat2mat(hmat = hmat,
                  h = NROW(hmat),
                  kset = tmp$set,
                  n = n)
  if(is.null(row_names) | length(row_names) != n){
    rownames(out) <- paste0("s-", 1:NROW(out))
  }else{
    rownames(out) <- row_names
  }
  return(out)
}


#' @return [as_hstack_ctlayout] returns a \eqn{h \times n(k^\ast+m)} numeric
#' matrix in \emph{horizon-stacked} layout (cross-temporal version).
#'
#' @rdname ctmatrix_layouts
#' @family Utilities
#' @export
as_hstack_ctlayout <- function(ctmat, agg_order){
  tmp <- tetools(agg_order = agg_order)
  n <- NROW(ctmat)
  if(NCOL(ctmat) %% tmp$dim[["kt"]] != 0){
    cli_abort(c("x" = "Incorrect {.arg ctmat} columns dimension.",
                "i" = "Please check the consistency of {.arg agg_order} and {.arg ctmat}."),
              call = NULL)
  }
  h <- NCOL(ctmat) / tmp$dim[["kt"]]
  out <- mat2hmat(mat = ctmat,
                  h = h,
                  kset = tmp$set,
                  n = n)
  rownames(out) <- paste0("tao-", 1:NROW(out))
  if(is.null(rownames(ctmat))){
    row_names <- paste0("s-", 1:n)
  }else{
    row_names <- rownames(ctmat)
  }
  colnames(out) <- as.vector(sapply(row_names, function(s)
    paste0(s, " ", namesTE(kset = tmp$set, h = 1))
  ))
  return(out)
}


#' Convert between horizon-stacked and temporal layouts
#'
#' These functions convert matrix between the two canonical layouts used in
#' temporal reconciliation.
#' Let \eqn{m} be the maximum temporal aggregation order and \eqn{k^\ast} the sum
#' of a subset of the \eqn{(p-1)} proper factors of \eqn{m} (excluding \eqn{m});
#' let \eqn{h} be the forecast horizon for the lowest frequency series (e.g.,
#' most aggregated temporal forecast horizon):
#' \itemize{
#'   \item \emph{Horizon-stacked layout (temporal version)}: a
#'   \eqn{h \times (k^\ast + m)} matrix where rows are the most aggregated temporal
#'   forecast horizons, and the values in each row are ordered from the lowest frequency
#'   (most temporally aggregated) to the highest frequency.
#'   \item \emph{Temporal layout}: a (\eqn{h(k^\ast + m) \times 1}) numeric vector where
#'   values are ordered from the lowest frequency (most temporally aggregated) to the
#'   highest frequency.
#' }
#' Then, [as_tevector] converts a \eqn{(h \times (k^\ast+m))}
#' horizon-stacked matrix to a (\eqn{h(k^\ast + m) \times 1}) temporal vector;
#' [as_hstack_telayout] performs the inverse transform.
#'
#' @param hmat A \eqn{h \times (k^\ast+m)} numeric matrix in \emph{horizon-stacked} layout
#' (temporal version).
#' @param tevec A (\eqn{h(k^\ast + m) \times 1}) numeric vector in \emph{temporal} layout.
#' @inheritParams ctcov
#'
#' @return [as_tevector] returns a (\eqn{h(k^\ast + m) \times 1}) numeric vector
#' in \emph{temporal} layout.
#'
#' @examples
#' h <- 2   # horizons
#' m <- 4   # temporal aggregation order
#' kt <- tetools(m)$dim["kt"]
#'
#' # Build a horizon-stacked matrix: h rows, n * k_t columns
#' input_te <- seq_len(h * kt)
#'
#' hmat <- as_hstack_telayout(input_te, agg_order = m)
#' tevec <- as_tevector(hmat, agg_order = m)
#' # all.equal(tevec, input_te, check.attributes = FALSE)
#'
#' @rdname tematrix_layouts
#' @family Utilities
#' @export
as_tevector <- function(hmat, agg_order){
  tmp <- tetools(agg_order = agg_order)

  if(NCOL(hmat) != tmp$dim["kt"]){
    cli_abort(c("x" = "Incorrect {.arg hmat} columns dimension.",
                "i" = "Expected {.val {tmp$dim['kt']}} columns, but found {.val {NCOL(hmat)}}.",
                "i" = "Please check the consistency of {.arg agg_order} and {.arg hmat}."),
              call = NULL)
  }

  out <- hmat2vec(hmat = hmat,
                  h = NROW(hmat),
                  kset = tmp$set)
  return(out)
}


#' @return [as_hstack_telayout] returns a \eqn{h \times (k^\ast+m)} numeric
#' matrix in \emph{horizon-stacked} layout (temporal version).
#'
#' @rdname tematrix_layouts
#' @export
as_hstack_telayout <- function(tevec, agg_order){
  tmp <- tetools(agg_order = agg_order)
  if(length(tevec) %% tmp$dim[["kt"]] != 0){
    cli_abort(c("x" = "Incorrect {.arg tevec} length.",
                "i" = "Please check the consistency of {.arg agg_order} and {.arg tevec}."),
              call = NULL)
  }
  h <- length(tevec) / tmp$dim[["kt"]]
  out <- vec2hmat(vec = tevec,
                  h = h,
                  kset = tmp$set)
  rownames(out) <- paste0("tao-", 1:NROW(out))
  colnames(out) <- namesTE(kset = tmp$set, h = 1)
  return(out)
}

