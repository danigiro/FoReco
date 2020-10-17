#' Cross-sectional (contemporaneous) aggregation matrix
#'
#'
#' This function allows the user to easily build the (\code{na x nb}) cross-sectional
#' (contemporaneous) matrix mapping the \code{nb} bottom level series into the \code{na} higher level
#' ones. (Experimental version)
#'
#' @param formula  Specification of the hierarchical structure: grouped hierarchies are specified
#' using \code{~ g1 * g2} and nested hierarchies are specified using \code{~ parent / child}.
#' Mixtures of the two formulations are also possible, like \code{~ g1 * (grandparent / parent / child)}.
#' @param data A dataset in which each column contains the values of the variables in the formula
#' and each row identifies a bottom level time series.
#' @param sep Character to separate the names of the aggregated series (\emph{default} is \code{"_"}).
#' @param sparse Option to return sparse matrix (\emph{default} is \code{TRUE}).
#' @param top_label Label of the top level variable (\emph{default} is \code{"Total"}).
#'
#' @return A (\code{na x nb}) matrix.
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
#' C <- Cmatrix(~ X1 / X2, data_bts, sep = "")
#'
#' ## Unbalanced hierarchy (1)
#' #             T
#' #    |--------|------|
#' #    A        B      C
#' #  |---|   |--|--|
#' # AA   AB  BA BB BC
#' # Names of the bottom level variables
#' data_bts <- data.frame(X1 = c("A", "A", "B", "B", "B", "C"),
#'                        X2 = c("A", "B", "A", "B", "C", NA),
#'                        stringsAsFactors = FALSE)
#' # Cross-sectional aggregation matrix
#' C <- Cmatrix(~ X1 / X2, data_bts, sep = "")
#'
#' ## Unbalanced hierarchy (2)
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
#' C <- Cmatrix(~ X1 / X2 / X3, data_bts, sep = "")
#'
#' ## Grouped hierarchy
#' #         C               S
#' #    |--------|      |--------|
#' #    A        B      M        F
#' #  |---|    |---|
#' # AA   AB  BA   BB
#' # Names of the bottom level variables
#' data_bts <- data.frame(X1 = c("A", "A", "B", "B", "A", "A", "B", "B"),
#'                        X2 = c("A", "B", "A", "B", "A", "B", "A", "B"),
#'                        Y1 = c("M", "M", "M", "M", "F", "F", "F", "F"),
#'                        stringsAsFactors = FALSE)
#' # Cross-sectional aggregation matrix
#' C <- Cmatrix(~ Y1 * (X1 / X2), data_bts, sep = "")
#'
#' @export
#' @import Matrix
Cmatrix <- function(formula, data, sep = "_", sparse = TRUE, top_label = "Total") {
  if (missing(data)) {
    stop("The data parameter is required")
  }

  if (NCOL(data) == 1) {
    out <- matrix(1, 1, NROW(unique(data)))
    rownames(out) <- "Total"
    colnames(out) <- unique(data[, 1, , drop = TRUE])
    return(out)
  }

  if (missing(formula)) {
    message(sprintf(
      "Formula is missing, defaulting to `formula = ~ %s`",
      paste(colnames(data), collapse = "*")
    ))
    formula <- stats::as.formula(paste("~", paste(colnames(data), collapse = "*"), sep = ""))
  }

  tm <- stats::terms(formula)
  lev <- attr(tm, "factors")
  lev_vars <- rownames(lev)
  lev <- Map(function(x) lev_vars[x != 0], split(lev, col(lev)))
  lev <- unname(lev)
  lev <- lev[lapply(lev, length)<length(lev_vars)]

  if (!all(lev_vars %in% colnames(data))) {
    stop("Please, data must data must contain all the variables in formula.")
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
  if (sparse) {
    out <- Matrix(out, sparse = TRUE)
  }
  rownames(out) <- namerows
  colnames(out) <- namecols

  out <- out[rowSums(out) > 1 & rowSums(out) < NCOL(out), , drop = FALSE]
  out <- rbind(Total = 1, out)
  rownames(out)[1] <- top_label
  message("------ Cross-sectional information ------")
  message("  Number of total time series (n): ", NROW(out) + NCOL(out))
  message(" Number of upper time series (na): ", NROW(out))
  message("Number of bottom time series (nb): ", NCOL(out))
  message("                 Number of levels: ", length(lev)+1)
  return(out)
}
