#' FoReco Reconciliation Class
#'
#' @description
#' The `foreco` class represents reconciled forecasts produced by the FoReco
#' package. It extends a numeric matrix, vector, or distributional object with
#' additional attributes that store metadata about the reconciliation procedure
#' (framework, function used, forecast type, and other reconciliation-specific
#' information). The class provides dedicated methods for printing,
#' summarising, extracting components and visualising reconciled forecasts.
#'
#' @param reco A numeric matrix/vector (when `rtype = "point"`) or a
#'   distributional object (when `rtype = "probabilistic"`).
#' @param x,object An object of class `foreco`.
#' @param framework A character string identifying the reconciliation framework.
#'   Must be one of `"cross-sectional"`, `"temporal"`, or `"cross-temporal"`.
#' @param rfun A character scalar with the name of the FoReco function that
#'   produced the reconciled forecasts (e.g. `"csrec"`, `"terec"`, `"ctrec"`).
#' @param rtype A character string indicating the type of reconciled forecasts.
#'   Must be one of `"point"` or `"probabilistic"`.
#' @param rinfo An optional named list with additional reconciliation
#'   information (e.g. covariance approximation `comb`, machine-learning
#'   approach `ml`, non-negativity flag `nn`, cross-sectional size `cs_n`,
#'   temporal aggregation set `te_set`, forecast horizon `forecast_horizon`,
#'   level-conditional reconciled forecasts `lcc`).
#' @param nninfo An optional matrix with information about the non-negativity
#'   procedure applied during reconciliation. Stored as the `info` field of
#'   the `"FoReco"` attribute.
#' @param keep_forecasts Logical; if `TRUE` (the default), the reconciled
#'   forecasts are stored in the `rf` element of the returned
#'   `summary_foreco` object (and therefore printed at the end of the summary).
#'   Set to `FALSE` to obtain a lighter summary object that omits them.
#' @param n_row,n_col Integers giving the maximum number of rows and columns to
#'   display when printing. If `NULL` (the default) all rows/columns are shown.
#'   When the matrix is truncated, a summary line reports how many rows and
#'   columns have been omitted.
#' @param cs Optional integer vector selecting the cross-sectional series to
#'   keep. If `NULL` (the default) all series are returned.
#' @param te Optional vector (numeric or character) selecting the temporal
#'   aggregation orders to keep, matched against the elements of `te_set`. If
#'   `NULL` (the default) all orders are returned.
#' @param alpha Nominal coverage of the prediction interval drawn by
#'   `plot.foreco()` for probabilistic forecasts. Defaults to `0.95`.
#' @param keep_names Logical. If `TRUE`, the row/column names of the reconciled
#'   forecasts are preserved in the output of `components.foreco()`. Defaults
#'   to `FALSE`.
#' @param temporal_names Optional character vector of labels for the temporal
#'   aggregation orders returned by `components.foreco()`. Its length must
#'   match the number of returned orders, otherwise a warning is emitted and
#'   the default `"k-..."` labels are used.
#' @param simplify Logical. If \code{TRUE} and the result consists of a single
#'   component, the underlying object (a matrix or vector) is returned directly
#'   instead of being wrapped in a named list of length one. If \code{FALSE}
#'   (default), the output is always a named list.
#' @param ... Additional arguments passed to the underlying methods
#'   (e.g. `print()`).
#'
#' @details
#' `new_foreco_class()` is the low-level constructor. It is exported so
#' that companion packages can produce objects that integrate with FoReco's
#' `print()`, `summary()`, `plot()` and `components()` methods.
#'
#' `plot.foreco()` draws the reconciled forecasts as line/point plots. For
#' probabilistic forecasts (`rtype = "probabilistic"`) it also overlays a
#' shaded central `alpha * 100%` prediction interval, built from the
#' `(1 - alpha)/2` and `1 - (1 - alpha)/2` quantiles of the distributional
#' object; the median is shown as a dashed line and the interval limits as
#' dotted lines.
#'
#' @return
#' A `foreco` object extending the reconciled forecasts/distributions with
#' reconciliation metadata.
#'
#' `components.foreco()` returns a named list of reconciled forecasts split
#' by temporal aggregation order. For the cross-sectional framework the list
#' has a single element `"k-1"`.
#'
#' @examples
#' set.seed(123)
#' # Aggregation matrix for Z = X + Y
#' A   <- t(c(1, 1))
#' bts <- matrix(rnorm(6, mean = 10), 3, 2)
#' reco <- csbu(base = bts, agg_mat = A)
#'
#' # Print and summarise the reconciled forecasts
#' print(reco)
#' print(reco, n_row = 2, n_col = 2)
#' summary(reco)
#' summary(reco, keep_forecasts = FALSE)
#'
#' # Extract reconciled forecasts by temporal aggregation order
#' components(reco)
#'
#' # Remove the foreco class
#' drop_foreco_class(reco)
#'
#' @name foreco-class
#' @aliases foreco summary.foreco print.foreco print.summary_foreco plot.foreco
#'   components.foreco drop_foreco_class
NULL

#' @rdname foreco-class
#' @export
new_foreco_class <- function(
  reco,
  framework,
  rfun,
  rtype,
  rinfo = NULL,
  nninfo = NULL
) {
  framework <- match.arg(
    framework,
    c("cross-sectional", "temporal", "cross-temporal")
  )

  rtype <- match.arg(
    rtype,
    c("point", "probabilistic")
  )

  if (!is.character(rfun) || length(rfun) != 1L) {
    cli_abort("{.arg rfun} must be a character", call = NULL)
  }

  if (!is.null(rinfo)) {
    if (!is.list(rinfo)) {
      cli_abort("{.arg rinfo} must be a list", call = NULL)
    }
  }

  if (!is.null(nninfo)) {
    if (!is.matrix(nninfo)) {
      cli_abort("{.arg nninfo} must be a matrix.", call = NULL)
    }
  }

  reserved <- c("framework", "rfun", "rtype", "info")
  if (any(reserved %in% names(rinfo))) {
    cli_abort(
      "{.arg rinfo} must not contain the reserved names {.arg framework}, {.arg rfun} and {.arg rtype}",
      call = NULL
    )
  }
  attr(reco, "FoReco") <- c(
    list(framework = framework, rfun = rfun, rtype = rtype, info = nninfo),
    rinfo
  )
  class(reco) <- c("foreco", class(reco))
  return(reco)
}

#' @rdname foreco-class
#' @method summary foreco
#' @export
summary.foreco <- function(object, keep_forecasts = TRUE, ...) {
  out <- attr(object, "FoReco")
  if (keep_forecasts) {
    out$object <- object
  }
  if (out$rtype == "point") {
    if (is.matrix(object)) {
      info_rf <- paste0("(", paste0(dim(object), collapse = " x "), ") matrix")
    } else if (is.vector(object)) {
      info_rf <- paste0("(", length(object), " x 1) vector")
    } else {
      info_rf <- "unknown object"
    }
  } else {
    info_rf <- "distributional object"
  }
  out$info_rf <- info_rf
  class(out) <- "summary_foreco"
  return(out)
}

#' @rdname foreco-class
#' @method print summary_foreco
#' @export
print.summary_foreco <- function(
  x,
  n_row = 4L,
  n_col = 6L,
  ...
) {
  #cli_rule(right = "FoReco reconciliation summary")
  frm <- paste0(
    toupper(substr(x$framework, 1, 1)),
    substr(x$framework, 2, nchar(x$framework))
  )
  cli_alert_success(
    "{.strong {frm}} {.emph {x$rtype}} forecast reconciliation"
  )

  method <- c()
  if (!is.null(x$rfun)) {
    method <- c(method, "Function used: {.strong {.code {x$rfun}}}")
  }

  if (!is.null(x$comb)) {
    if (!all(is.character(x$comb))) {
      tmp <- "custom"
    } else if (length(x$comb) > 1) {
      tmp <- paste(names(x$comb), x$comb, sep = "-")
    } else {
      tmp <- x$comb
    }
    method <- c(
      method,
      "Covariance approximation approach: {.strong {.code {tmp}}}"
    )
  }

  if (!is.null(x$ml)) {
    method <- c(method, "Machine Learning approach: {.strong {.code {x$ml}}}")
  }

  if (!is.null(x$info_rf)) {
    method <- c(method, "Output: {x$info_rf}")
  }

  str_frm <- c()
  if (!is.null(x$cs_n)) {
    str_frm <- c(str_frm, "Number of cross-sectional series: {x$cs_n}")
  }

  if (!is.null(x$te_set)) {
    str_frm <- c(
      str_frm,
      "Temporal orders (k): {x$te_set}",
      "Forecast horizons (h) per k: {x$forecast_horizon * max(x$te_set)/x$te_set}"
    )
  } else {
    str_frm <- c(str_frm, "Forecast horizons (h): {x$forecast_horizon}")
  }

  if (!is.null(x$nn)) {
    str_frm <- c(
      str_frm,
      "Non-negative forecasts (check): {.strong {.code {x$nn}}}"
    )
  }

  cli_h3("Method")
  cli_ul(method)
  cli_h3("Structure")
  cli_ul(str_frm)

  if (!is.null(x$info)) {
    cli_h3("Non-negative reconciliation diagnostics")
    print(head(x$info, n = 5))
    if (NROW(x$info) > 5) {
      cli_alert_info(
        "Showing the first 5 rows of the non-negativity diagnostics info matrix."
      )
    }
  }

  if (!is.null(x$object)) {
    cli_h3("Reconciled forecasts")
    print_foreco(
      x$object,
      n_row = n_row,
      n_col = n_col,
      .caller = "print",
      .name = deparse(substitute(x))
    )
  }

  invisible(x)
}

style_comment <- cli::make_ansi_style(
  grDevices::grey(0.6),
  grey = TRUE,
  colors = 256
)

#' @rdname foreco-class
#' @method print foreco
#' @export
print.foreco <- function(x, n_row = NULL, n_col = NULL, ...) {
  print_foreco(
    x = x,
    n_row = n_row,
    n_col = n_col,
    .name = deparse(substitute(x)),
    ...
  )
}

#' @rdname foreco-class
#' @method plot foreco
#' @export
plot.foreco <- function(x, cs = NULL, te = 1, alpha = 0.95, ...) {
  fr <- summary(x)
  if (fr$framework != "cross-sectional" & is.null(te)) {
    cli_warn(
      c(
        "No temporal order specified via {.arg te}.",
        "i" = "Plotting the most disaggregated series (k = 1) by default.",
        "i" = "Pass a value to {.arg te} to choose a different order."
      ),
      call = NULL
    )
  }
  if (fr$rtype == "point") {
    if (fr$framework == "cross-sectional") {
      mat <- .drop_foreco(x)
      if (!is.null(cs)) {
        plot_x <- mat[, cs, drop = FALSE]
      } else {
        plot_x <- mat
      }
      ylab <- "Reconciled forecasts"
    } else {
      lcomp <- components(x, cs = cs, te = te)
      plot_x <- lcomp[[length(lcomp)]]
      ylab <- paste0(
        "Reconciled forecasts with k = ",
        sub("k-", "", names(lcomp)[length(lcomp)])
      )
    }
    matplot(
      plot_x,
      type = "l",
      ylab = ylab,
      xlab = "Forecasts horizons"
    )
    matpoints(plot_x, pch = 19)
  } else {
    if (
      !is.numeric(alpha) ||
        length(alpha) != 1L ||
        !is.finite(alpha) ||
        alpha <= 0 ||
        alpha >= 1
    ) {
      cli_abort(
        "{.arg alpha} must be a single number in {.code (0, 1)}, got {.val {alpha}}.",
        call = NULL
      )
    }
    ci_lo <- (1 - alpha) / 2
    ci_hi <- 1 - ci_lo
    if (fr$framework == "cross-sectional") {
      mat <- median(x)
      up_mat <- quantile(x, p = ci_hi)
      lw_mat <- quantile(x, p = ci_lo)
      if (!is.null(cs)) {
        mat <- mat[, cs, drop = FALSE]
        up_mat <- up_mat[, cs, drop = FALSE]
        lw_mat <- lw_mat[, cs, drop = FALSE]
      }
      ylab <- c("Reconciled forecasts")
    } else if (fr$framework == "cross-temporal") {
      mat <- as_ctmatrix(median(x), agg_order = fr$te_set, n = fr$cs_n)
      mat <- distr_to_point(mat, attr(x, "FoReco"))
      mat <- components(mat, cs = cs, te = te)
      ylab <- paste0(
        "Reconciled forecasts with k = ",
        sub("k-", "", names(mat)[length(mat)])
      )
      mat <- mat[[length(mat)]]

      up_mat <- as_ctmatrix(
        quantile(x, p = ci_hi),
        agg_order = fr$te_set,
        n = fr$cs_n
      )
      up_mat <- distr_to_point(up_mat, attr(x, "FoReco"))
      up_mat <- components(up_mat, cs = cs, te = te)
      up_mat <- up_mat[[length(up_mat)]]

      lw_mat <- as_ctmatrix(
        quantile(x, ci_lo),
        agg_order = fr$te_set,
        n = fr$cs_n
      )
      lw_mat <- distr_to_point(lw_mat, attr(x, "FoReco"))
      lw_mat <- components(lw_mat, cs = cs, te = te)
      lw_mat <- lw_mat[[length(lw_mat)]]
    } else {
      mat <- as_tevector(median(x), agg_order = fr$te_set)
      mat <- distr_to_point(mat, attr(x, "FoReco"))
      mat <- components(mat, te = te)
      ylab <- paste0(
        "Reconciled forecasts with k = ",
        sub("k-", "", names(mat)[length(mat)])
      )
      mat <- mat[[length(mat)]]

      up_mat <- as_tevector(quantile(x, p = ci_hi), agg_order = fr$te_set)
      up_mat <- distr_to_point(up_mat, attr(x, "FoReco"))
      up_mat <- components(up_mat, te = te)
      up_mat <- up_mat[[length(up_mat)]]

      lw_mat <- as_tevector(quantile(x, ci_lo), agg_order = fr$te_set)
      lw_mat <- distr_to_point(lw_mat, attr(x, "FoReco"))
      lw_mat <- components(lw_mat, te = te)
      lw_mat <- lw_mat[[length(lw_mat)]]
    }
    xh <- 1:NROW(mat)
    xx <- as.vector(rbind(
      matrix(xh, NROW(mat), NCOL(mat)),
      matrix(rev(xh), NROW(mat), NCOL(mat)),
      NA
    ))
    if (is.matrix(lw_mat)) {
      yy <- as.vector(rbind(up_mat, apply(lw_mat, 2, rev), NA))
    } else {
      yy <- as.vector(c(up_mat, rev(lw_mat), NA))
    }

    matplot(
      xh,
      mat,
      type = "n",
      ylim = range(c(lw_mat, up_mat), na.rm = TRUE),
      ylab = ylab,
      xlab = "Forecasts horizons"
    )
    polygon(xx, yy, col = "grey80", border = NA)
    matlines(xh, mat, lty = 2)
    matlines(xh, up_mat, lty = 3)
    matlines(xh, lw_mat, lty = 3)
    matpoints(xh, mat, pch = 19)
  }
}

#' @importFrom generics components
#' @export
generics::components

#' @rdname foreco-class
#' @method components foreco
#' @export
components.foreco <- function(
  object,
  cs = NULL,
  te = NULL,
  keep_names = FALSE,
  temporal_names = NULL,
  simplify = FALSE,
  ...
) {
  fr <- summary(object)

  frame <- fr$framework
  set <- fr$te_set
  h <- fr$forecast_horizon

  if (fr$rtype == "point") {
    x <- .drop_foreco(object)
  } else {
    cli_abort("{.fn components} is not supported for probabilistic forecasts.")
  }

  if (frame == "cross-sectional") {
    attr(x, "FoReco") <- NULL

    if (!is.null(cs)) {
      x <- x[, cs, drop = FALSE]
    }

    if (simplify) {
      return(x)
    } else {
      return(list("k-1" = x))
    }
  } else {
    id <- rep(set, h * max(set) / set)

    if (NCOL(x) == 1) {
      out <- split(x, factor(id, set))
      if (!keep_names) {
        out <- lapply(out, unname)
      }
    } else {
      out <- lapply(setNames(set, set), function(k) {
        mat <- t(x[, id == k, drop = FALSE])
        if (!keep_names) {
          rownames(mat) <- NULL
        }
        if (!is.null(cs)) {
          mat <- mat[, cs, drop = FALSE]
        }

        return(mat)
      })
    }

    if (!is.null(te)) {
      extract_k <- which(names(out) %in% te)
      if (length(extract_k) == 0) {
        cli_abort(
          c(
            "{.arg te} does not match any temporal aggregation order.",
            "i" = "Available orders: {.val {names(out)}}.",
            "x" = "Got: {.val {as.character(te)}}."
          ),
          call = NULL
        )
      } else {
        names(out) <- paste0("k-", names(out))
        out <- out[extract_k]
      }
    } else {
      names(out) <- paste0("k-", names(out))
    }

    if (!is.null(temporal_names)) {
      if (length(temporal_names) == length(out)) {
        names(out) <- paste0(temporal_names, " (", names(out), ")")
      } else {
        cli_warn(
          c(
            paste0(
              "Length of {.arg temporal_names} ({length(temporal_names)})",
              " does not match the number of returned temporal aggregation",
              " orders ({length(out)})."
            ),
            "i" = "Default labels {.val {names(out)}} are used instead."
          ),
          call = NULL
        )
      }
    }

    if (simplify && length(out) == 1) {
      return(out[[1]])
    } else {
      return(out)
    }
  }
}

#' @method as.matrix foreco
#' @export
as.matrix.foreco <- function(x, ...) {
  as.matrix(.drop_foreco(x))
}

#' @method as.data.frame foreco
#' @export
as.data.frame.foreco <- function(x, ...) {
  as.data.frame(.drop_foreco(x), ...)
}

#' @method Ops foreco
#' @export
Ops.foreco <- function(e1, e2) {
  if (inherits(e1, "foreco")) {
    e1 <- .drop_foreco(e1)
  }
  if (!missing(e2) && inherits(e2, "foreco")) {
    e2 <- .drop_foreco(e2)
  }
  NextMethod()
}

#' @method Math foreco
#' @export
Math.foreco <- function(x, ...) {
  NextMethod(.drop_foreco(x))
}

#' @rdname foreco-class
#' @export
drop_foreco_class <- function(x) {
  .drop_foreco(x)
}

.drop_foreco <- function(x) {
  if (inherits(x, "foreco")) {
    attr(x, "FoReco") <- NULL
    cl <- setdiff(class(x), "foreco")
    if (length(cl) == 0) {
      x <- unclass(x)
    } else {
      class(x) <- cl
    }
  }
  x
}

distr_to_point <- function(x, attr_info) {
  attr(x, "FoReco") <- attr_info
  attr(x, "FoReco")$rtype <- "point"
  class(x) <- c("foreco", class(x))
  return(x)
}


print_foreco <- function(
  x,
  n_row = NULL,
  n_col = NULL,
  .caller = "print",
  .name = NULL,
  ...
) {
  if (is.null(.name)) {
    name_x <- deparse(substitute(x))
  } else {
    name_x <- .name
  }
  class_x <- class(x)
  if (is.matrix(x)) {
    attr(x, "FoReco") <- NULL

    check_null_values <- is.null(n_row) || is.null(n_col)
    nr <- nrow(x)
    nc <- ncol(x)
    n_row <- if (is.null(n_row)) nr else min(as.integer(n_row), nr)
    n_col <- if (is.null(n_col)) nc else min(as.integer(n_col), nc)

    if (n_row < nr || n_col < nc) {
      x <- x[seq_len(n_row), seq_len(n_col), drop = FALSE]
    }

    print(.drop_foreco(x), ...)
    if (n_row < nr || n_col < nc) {
      cat(
        "... (",
        nr - n_row,
        " more row",
        if (nr - n_row != 1L) "s" else "",
        ", ",
        nc - n_col,
        " more column",
        if (nc - n_col != 1L) "s" else "",
        ")\n",
        style_comment(
          paste0(
            "Use `",
            .caller,
            "(",
            name_x,
            ", n_row, n_col)` to see more rows and columns."
          )
        ),
        "\n",
        sep = ""
      )
    } else if (check_null_values) {
      cat(
        style_comment(
          paste0(
            "All rows and columns are shown.\n",
            "Use `print(",
            name_x,
            ", n_row, n_col)` to limit the output."
          )
        ),
        "\n"
      )
    }
  } else if (NCOL(x) == 1) {
    attr(x, "FoReco") <- NULL

    check_null_values <- is.null(n_row)
    nr <- length(x)
    n_row <- if (is.null(n_row)) nr else min(as.integer(n_row), nr)

    if (n_row < nr) {
      x <- x[seq_len(n_row)]
    }

    print(.drop_foreco(x), ...)
    if (n_row < nr) {
      cat(
        "... (",
        nr - n_row,
        " more element",
        if (nr - n_row != 1L) "s" else "",
        ")\n",
        style_comment(
          paste0(
            "Use `print(",
            name_x,
            ", n_row)` to see more elements."
          )
        ),
        "\n",
        sep = ""
      )
    } else if (check_null_values) {
      cat(
        style_comment(
          paste0(
            "All elements are shown. ",
            "Use `print(",
            name_x,
            ", n_row)` to limit the output."
          )
        ),
        "\n"
      )
    }
  } else {
    attr(x, "FoReco") <- NULL
    print(.drop_foreco(x))
  }
  invisible(x)
}
