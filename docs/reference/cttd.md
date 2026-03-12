# Cross-temporal top-down reconciliation

Top-down forecast reconciliation for cross-temporal hierarchical/grouped
time series, where the forecast of a \`Total' (top-level series,
expected to be positive) is disaggregated according to a proportional
scheme (weights). Besides fulfilling any aggregation constraint, the
top-down reconciled forecasts should respect two main properties:

- the top-level value remains unchanged;

- all the bottom time series reconciled forecasts are non-negative.

## Usage

``` r
cttd(base, agg_mat, agg_order, weights, tew = "sum", normalize = TRUE)
```

## Arguments

- base:

  A (\\hm \times 1\\) numeric vector containing top- and \\m\\ temporal
  aggregated level base forecasts; \\m\\ is the max aggregation order,
  and \\h\\ is the forecast horizon for the lowest frequency time
  series.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- weights:

  A (\\n_b \times hm\\) numeric matrix containing the proportions for
  each high-frequency bottom time series; \\n_b\\ is the total number of
  high-frequency bottom variables, \\m\\ is the max aggregation order,
  and \\h\\ is the forecast horizon for the lowest frequency time
  series.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- normalize:

  If `TRUE` (*default*), the `weights` will sum to 1.

## Value

A (\\n \times h(k^\ast+m)\\) numeric matrix of cross-temporal reconciled
forecasts.

## See also

Top-down reconciliation:
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

## Examples

``` r
set.seed(123)
# (3 x 1) top base forecasts vector (simulated), forecast horizon = 3
topf <- rnorm(3, 10)
A <- t(c(1,1)) # Aggregation matrix for Z = X + Y

# Same weights for different forecast horizons, agg_order = 4
fix_weights <- matrix(runif(4*2), 2, 4)
reco <- cttd(base = topf, agg_mat = A, agg_order = 4, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- matrix(runif(4*2*3), 2, 3*4)
recoh <- cttd(base = topf, agg_mat = A, agg_order = 4, weights = h_weights)
```
