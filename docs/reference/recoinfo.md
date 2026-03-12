# Informations on the reconciliation process

This function extracts reconciliation information from the output of any
reconciled function implemented by FoReco.

## Usage

``` r
recoinfo(x, verbose = TRUE)
```

## Arguments

- x:

  An output from any reconciliation function implemented by FoReco.

- verbose:

  If `TRUE` (*defaults*), reconciliation information are printed.

## Value

A list containing the following reconciliation process information:

- rfun:

  the reconciliation function.

- cs_n:

  the cross-sectional number of variables.

- te_set:

  the set of temporal aggregation orders.

- forecast_horizon:

  the forecast horizon (in temporal and cross-temporal frameworks, for
  the most temporally aggregated series).

- framework:

  the reconciliation framework (cross-sectional, temporal or
  cross-temporal).

- info:

  non-negative reconciled forecast convergence information.

- lcc:

  list of level conditional reconciled forecasts (+ BU) for
  [cslcc](https://danigiro.github.io/FoReco/reference/cslcc.md),
  [telcc](https://danigiro.github.io/FoReco/reference/telcc.md) and
  [ctlcc](https://danigiro.github.io/FoReco/reference/ctlcc.md).

- nn:

  if `TRUE`, all the forecasts are not negative.

- comb:

  the covariance approximation.

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
[`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md),
[`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md),
[`balance_hierarchy()`](https://danigiro.github.io/FoReco/reference/balance_hierarchy.md),
[`commat()`](https://danigiro.github.io/FoReco/reference/commat.md),
[`csprojmat()`](https://danigiro.github.io/FoReco/reference/csprojmat.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md),
[`ctprojmat()`](https://danigiro.github.io/FoReco/reference/ctprojmat.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`df2aggmat()`](https://danigiro.github.io/FoReco/reference/df2aggmat.md),
[`lcmat()`](https://danigiro.github.io/FoReco/reference/lcmat.md),
[`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md),
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)
