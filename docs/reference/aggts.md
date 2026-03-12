# Non-overlapping temporal aggregation of a time series

Non-overlapping temporal aggregation of a time series according to a
specific aggregation order.

## Usage

``` r
aggts(y, agg_order, tew = "sum", align = "end", rm_na = FALSE)
```

## Arguments

- y:

  Univariate or multivariate time series: a vector/matrix or a `ts`
  object.

- agg_order:

  A numeric vector with the aggregation orders to consider.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- align:

  A string or a vector specifying the alignment of `y`. Options include:
  "`end`" (end of the series, *default*), "`start`" (start of the
  series), an integer (or a vector of integers) indicating the starting
  period of the temporally aggregated series.

- rm_na:

  If `TRUE` the missing values are removed.

## Value

A list of vectors or `ts` objects.

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
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
[`recoinfo()`](https://danigiro.github.io/FoReco/reference/recoinfo.md),
[`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md),
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
# Monthly time series (input vector)
y <- ts(rnorm(24), start = 2020, frequency = 12)
# Quarterly time series
x1 <- aggts(y, 3)
# Monthly, quarterly and annual time series
x2 <- aggts(y, c(1, 3, 12))
# All temporally aggregated time series
x3 <- aggts(y)

# Ragged data
y2 <- ts(rnorm(11), start = c(2020, 3), frequency = 4)
# Annual time series: start in 2021
x4 <- aggts(y2, 4, align = 3)
# Semi-annual (start in 2nd semester of 2020) and annual (start in 2021)
# time series
x5 <- aggts(y2, c(2, 4), align = c(1, 3))
```
