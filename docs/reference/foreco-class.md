# FoReco Reconciliation Class

The `foreco` class represents reconciled forecasts produced by the
FoReco package. It extends a numeric matrix, vector, or distributional
object with additional attributes that store metadata about the
reconciliation procedure (framework, function used, forecast type, and
other reconciliation-specific information). The class provides dedicated
methods for printing, summarising, extracting components and visualising
reconciled forecasts.

## Usage

``` r
new_foreco_class(reco, framework, rfun, rtype, rinfo = NULL, nninfo = NULL)

# S3 method for class 'foreco'
summary(object, keep_forecasts = TRUE, ...)

# S3 method for class 'summary_foreco'
print(x, n_row = 4L, n_col = 6L, ...)

# S3 method for class 'foreco'
print(x, n_row = NULL, n_col = NULL, ...)

# S3 method for class 'foreco'
plot(x, cs = NULL, te = 1, alpha = 0.95, ...)

# S3 method for class 'foreco'
components(
  object,
  cs = NULL,
  te = NULL,
  keep_names = FALSE,
  temporal_names = NULL,
  simplify = FALSE,
  ...
)

drop_foreco_class(x)
```

## Arguments

- reco:

  A numeric matrix/vector (when `rtype = "point"`) or a distributional
  object (when `rtype = "probabilistic"`).

- framework:

  A character string identifying the reconciliation framework. Must be
  one of `"cross-sectional"`, `"temporal"`, or `"cross-temporal"`.

- rfun:

  A character scalar with the name of the FoReco function that produced
  the reconciled forecasts (e.g. `"csrec"`, `"terec"`, `"ctrec"`).

- rtype:

  A character string indicating the type of reconciled forecasts. Must
  be one of `"point"` or `"probabilistic"`.

- rinfo:

  An optional named list with additional reconciliation information
  (e.g. covariance approximation `comb`, machine-learning approach `ml`,
  non-negativity flag `nn`, cross-sectional size `cs_n`, temporal
  aggregation set `te_set`, forecast horizon `forecast_horizon`,
  level-conditional reconciled forecasts `lcc`).

- nninfo:

  An optional matrix with information about the non-negativity procedure
  applied during reconciliation. Stored as the `info` field of the
  `"FoReco"` attribute.

- keep_forecasts:

  Logical; if `TRUE` (the default), the reconciled forecasts are stored
  in the `rf` element of the returned `summary_foreco` object (and
  therefore printed at the end of the summary). Set to `FALSE` to obtain
  a lighter summary object that omits them.

- ...:

  Additional arguments passed to the underlying methods (e.g.
  [`print()`](https://rdrr.io/r/base/print.html)).

- x, object:

  An object of class `foreco`.

- n_row, n_col:

  Integers giving the maximum number of rows and columns to display when
  printing. If `NULL` (the default) all rows/columns are shown. When the
  matrix is truncated, a summary line reports how many rows and columns
  have been omitted.

- cs:

  Optional integer vector selecting the cross-sectional series to keep.
  If `NULL` (the default) all series are returned.

- te:

  Optional vector (numeric or character) selecting the temporal
  aggregation orders to keep, matched against the elements of `te_set`.
  If `NULL` (the default) all orders are returned.

- alpha:

  Nominal coverage of the prediction interval drawn by `plot.foreco()`
  for probabilistic forecasts. Defaults to `0.95`.

- keep_names:

  Logical. If `TRUE`, the row/column names of the reconciled forecasts
  are preserved in the output of `components.foreco()`. Defaults to
  `FALSE`.

- temporal_names:

  Optional character vector of labels for the temporal aggregation
  orders returned by `components.foreco()`. Its length must match the
  number of returned orders, otherwise a warning is emitted and the
  default `"k-..."` labels are used.

- simplify:

  Logical. If `TRUE` and the result consists of a single component, the
  underlying object (a matrix or vector) is returned directly instead of
  being wrapped in a named list of length one. If `FALSE` (default), the
  output is always a named list.

## Value

A `foreco` object extending the reconciled forecasts/distributions with
reconciliation metadata.

`components.foreco()` returns a named list of reconciled forecasts split
by temporal aggregation order. For the cross-sectional framework the
list has a single element `"k-1"`.

## Details

`new_foreco_class()` is the low-level constructor. It is exported so
that companion packages can produce objects that integrate with FoReco's
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`components()`](https://generics.r-lib.org/reference/components.html)
methods.

`plot.foreco()` draws the reconciled forecasts as line/point plots. For
probabilistic forecasts (`rtype = "probabilistic"`) it also overlays a
shaded central `alpha * 100%` prediction interval, built from the
`(1 - alpha)/2` and `1 - (1 - alpha)/2` quantiles of the distributional
object; the median is shown as a dashed line and the interval limits as
dotted lines.

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y
A   <- t(c(1, 1))
bts <- matrix(rnorm(6, mean = 10), 3, 2)
reco <- csbu(base = bts, agg_mat = A)

# Print and summarise the reconciled forecasts
print(reco)
#>          s-1       s-2      s-3
#> h-1 19.51003  9.439524 10.07051
#> h-2 19.89911  9.769823 10.12929
#> h-3 23.27377 11.558708 11.71506
#> All rows and columns are shown.
#> Use `print(reco, n_row, n_col)` to limit the output. 
print(reco, n_row = 2, n_col = 2)
#>          s-1      s-2
#> h-1 19.51003 9.439524
#> h-2 19.89911 9.769823
#> ... (1 more row, 1 more column)
#> Use `print(reco, n_row, n_col)` to see more rows and columns.
summary(reco)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method 
#> • Function used: `csbu`
#> • Output: (3 x 3) matrix
#> 
#> ── Structure 
#> • Number of cross-sectional series: 3
#> • Forecast horizons (h): 3
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts 
#>          s-1       s-2      s-3
#> h-1 19.51003  9.439524 10.07051
#> h-2 19.89911  9.769823 10.12929
#> h-3 23.27377 11.558708 11.71506
summary(reco, keep_forecasts = FALSE)
#> ✔ Cross-sectional point forecast reconciliation
#> 
#> ── Method 
#> • Function used: `csbu`
#> • Output: (3 x 3) matrix
#> 
#> ── Structure 
#> • Number of cross-sectional series: 3
#> • Forecast horizons (h): 3
#> • Non-negative forecasts (check): `TRUE`

# Extract reconciled forecasts by temporal aggregation order
components(reco)
#> $`k-1`
#>          s-1       s-2      s-3
#> h-1 19.51003  9.439524 10.07051
#> h-2 19.89911  9.769823 10.12929
#> h-3 23.27377 11.558708 11.71506
#> 

# Remove the foreco class
drop_foreco_class(reco)
#>          s-1       s-2      s-3
#> h-1 19.51003  9.439524 10.07051
#> h-2 19.89911  9.769823 10.12929
#> h-3 23.27377 11.558708 11.71506
```
