# Temporal middle-out reconciliation

The middle-out forecast reconciliation for temporal hierarchies combines
top-down ([tetd](https://danigiro.github.io/FoReco/reference/tetd.md))
and bottom-up
([tebu](https://danigiro.github.io/FoReco/reference/tebu.md)) methods.
Given the base forecasts of an intermediate temporal aggregation order
\\k\\, it performs

- a top-down approach for the aggregation orders \\\<k\\;

- a bottom-up approach for the aggregation orders \\\>k\\.

## Usage

``` r
temo(base, agg_order, weights, order = max(agg_order), tew = "sum",
     normalize = TRUE)
```

## Arguments

- base:

  A (\\hk \times 1\\) numeric vector containing the temporal aggregated
  base forecasts of order \\k\\; \\k\\ is an aggregation order (a factor
  of \\m\\, and \\1\<k\<m\\), \\m\\ is the max aggregation order, and
  \\h\\ is the forecast horizon for the lowest frequency time series.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- weights:

  A (\\hm \times 1\\) numeric vector containing the proportions for the
  high-frequency time series; \\m\\ is the max aggregation order, and
  \\h\\ is the forecast horizon for the lowest frequency time series.

- order:

  The intermediate fixed aggregation order \\k\\.

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

Middle-out reconciliation:
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
# (6 x 1) base forecasts vector (simulated), forecast horizon = 3
# and intermediate aggregation order k = 2 (max agg order = 4)
basek2 <- rnorm(3*2, 5)
# Same weights for different forecast horizons
fix_weights <- runif(4)
reco <- temo(base = basek2, order = 2, agg_order = 4, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- runif(4*3)
recoh <- temo(base = basek2, order = 2, agg_order = 4, weights = h_weights)
```
