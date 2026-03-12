# Temporal top-down reconciliation

Top-down forecast reconciliation for a univariate time series, where the
forecast of the most aggregated temporal level is disaggregated
according to a proportional scheme (weights). Besides fulfilling any
aggregation constraint, the top-down reconciled forecasts should respect
two main properties:

- the top-level value remains unchanged;

- all the bottom time series reconciled forecasts are non-negative.

## Usage

``` r
tetd(base, agg_order, weights, tew = "sum", normalize = TRUE)
```

## Arguments

- base:

  A (\\hm \times 1\\) numeric vector containing the temporal aggregated
  base forecasts of order \\m\\; \\m\\ is the max aggregation order, and
  \\h\\ is the forecast horizon for the lowest frequency time series.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- weights:

  A (\\hm \times 1\\) numeric vector containing the proportions for the
  high-frequency time series; \\m\\ is the max aggregation order, and
  \\h\\ is the forecast horizon for the lowest frequency time series.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- normalize:

  If `TRUE` (*default*), the `weights` will sum to 1.

## Value

A (\\h(k^\ast+m) \times 1\\) numeric vector of temporal reconciled
forecasts.

## See also

Top-down reconciliation:
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
# (2 x 1) top base forecasts vector (simulated), forecast horizon = 2
topf <- rnorm(2, 10)
# Same weights for different forecast horizons
fix_weights <- runif(4)
reco <- tetd(base = topf, agg_order = 4, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- runif(4*2)
recoh <- tetd(base = topf, agg_order = 4, weights = h_weights)
```
