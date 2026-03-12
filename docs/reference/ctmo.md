# Cross-temporal middle-out reconciliation

The cross-temporal middle-out forecast reconciliation combines top-down
([cttd](https://danigiro.github.io/FoReco/reference/cttd.md)) and
bottom-up ([ctbu](https://danigiro.github.io/FoReco/reference/ctbu.md))
methods in the cross-temporal framework for genuine hierarchical/grouped
time series. Given the base forecasts of an intermediate cross-sectional
level \\l\\ and aggregation order \\k\\, it performs

- a top-down approach for the aggregation orders \\\geq k\\ and
  cross-sectional levels \\\geq l\\;

- a bottom-up approach, otherwise.

## Usage

``` r
ctmo(base, agg_mat, agg_order, weights, id_rows = 1, order = max(agg_order),
     tew = "sum", normalize = TRUE)
```

## Arguments

- base:

  A (\\n_l \times hk\\) numeric matrix containing the \\l\\-level base
  forecasts of temporal aggregation order \\k\\; \\n_l\\ is the number
  of variables at level \\l\\, \\k\\ is an aggregation order (a factor
  of \\m\\, and \\1\<k\<m\\), \\m\\ is the max aggregation order, and
  \\h\\ is the forecast horizon for the lowest frequency time series.

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

- id_rows:

  A numeric vector indicating the \\l\\-level rows of `agg_mat`.

- order:

  The intermediate fixed aggregation order \\k\\.

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

Middle-out reconciliation:
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

## Examples

``` r
set.seed(123)
# Aggregation matrix for Z = X + Y, X = XX + XY and Y = YX + YY
A <- matrix(c(1,1,1,1,1,1,0,0,0,0,1,1), 3, byrow = TRUE)
# (2 x 6) base forecasts matrix (simulated), forecast horizon = 3
# and intermediate aggregation order k = 2 (max agg order = 4)
baseL2k2 <- rbind(rnorm(3*2, 5), rnorm(3*2, 5))

# Same weights for different forecast horizons, agg_order = 4
fix_weights <- matrix(runif(4*4), 4, 4)
reco <- ctmo(base = baseL2k2, id_rows = 2:3, agg_mat = A,
             order = 2, agg_order = 4, weights = fix_weights)

# Different weights for different forecast horizons
h_weights <- matrix(runif(4*4*3), 4, 3*4)
recoh <- ctmo(base = baseL2k2, id_rows = 2:3, agg_mat = A,
             order = 2, agg_order = 4, weights = h_weights)
```
