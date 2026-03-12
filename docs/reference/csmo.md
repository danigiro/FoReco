# Cross-sectional middle-out reconciliation

The middle-out forecast reconciliation (Athanasopoulos et al., 2009)
combines top-down
([cstd](https://danigiro.github.io/FoReco/reference/cstd.md)) and
bottom-up ([csbu](https://danigiro.github.io/FoReco/reference/csbu.md))
for genuine hierarchical/grouped time series. Given the base forecasts
of variables at an intermediate level \\l\\, it performs

- a top-down approach for the levels \\\<l\\;

- a bottom-up approach for the levels \\\>l\\.

## Usage

``` r
csmo(base, agg_mat, weights, id_rows = 1, normalize = TRUE)
```

## Arguments

- base:

  A (\\h \times n_l\\) numeric matrix containing the \\l\\-level base
  forecast; \\n_l\\ is the number of variables at level \\l\\, and \\h\\
  is the forecast horizon.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- weights:

  A (\\h \times n_b\\) numeric matrix containing the proportions for the
  bottom time series; \\h\\ is the forecast horizon, and \\n_b\\ is the
  total number of bottom variables.

- id_rows:

  A numeric vector indicating the \\l\\-level rows of `agg_mat`.

- normalize:

  If `TRUE` (*default*), the `weights` will sum to 1.

## Value

A (\\h \times n\\) numeric matrix of cross-sectional reconciled
forecasts.

## References

Athanasopoulos, G., Ahmed, R. A. and Hyndman, R.J. (2009) Hierarchical
forecasts for Australian domestic tourism. *International Journal of
Forecasting* 25(1), 146–166.
[doi:10.1016/j.ijforecast.2008.07.004](https://doi.org/10.1016/j.ijforecast.2008.07.004)

## See also

Middle-out reconciliation:
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md)

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
# (3 x 2) top base forecasts vector (simulated), forecast horizon = 3
baseL2 <- matrix(rnorm(2*3, 5), 3, 2)
# Same weights for different forecast horizons
fix_weights <- runif(4)
reco <- csmo(base = baseL2, agg_mat = A, id_rows = 2:3, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- matrix(runif(4*3), 3, 4)
recoh <- csmo(base = baseL2, agg_mat = A, id_rows = 2:3, weights = h_weights)
```
