# Temporal probabilistic reconciliation (sample approach)

This function performs temporal probabilistic forecast reconciliation
using a sample-based approach (Girolimetto et al., 2024) for a single
time series using temporal hierarchies (Athanasopoulos et al., 2017).
Given a \\(L \times h(k^\ast + m))\\ matrix of simulated base forecast
draws, `tesmp()` applies a chosen FoReco temporal reconciliation
independently to each draw, producing a coherent sample distribution of
reconciled forecasts across the temporal hierarchy. Typical choices for
the reconciliation include optimal combination
([terec](https://danigiro.github.io/FoReco/reference/terec.md)) as well
as top-down
([tetd](https://danigiro.github.io/FoReco/reference/tetd.md)),
middle-out
([temo](https://danigiro.github.io/FoReco/reference/temo.md)), bottom-up
([tebu](https://danigiro.github.io/FoReco/reference/tebu.md)), and
level-conditional
([telcc](https://danigiro.github.io/FoReco/reference/telcc.md))
approaches.

## Usage

``` r
tesmp(sample, agg_order, fun = terec, ...)
```

## Arguments

- sample:

  A (\\L \times h(k^\ast + m)\\) numeric matrix containing the base
  forecasts samples to be reconciled; \\m\\ is the max aggregation
  order, \\k^\ast\\ is the sum of (a subset of) (\\p-1\\) factors of
  \\m\\, excluding \\m\\, \\h\\ is the forecast horizon for the lowest
  frequency time series, and \\L\\ is the sample size.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- fun:

  A string specifying the reconciliation function to be used, as
  implemented in FoReco.

- ...:

  Arguments passed on to `fun`

## Value

A
[distributional::dist_sample](https://pkg.mitchelloharawild.com/distributional/reference/dist_sample.html)
object.

## References

Athanasopoulos, G., Hyndman, R.J., Kourentzes, N. and Petropoulos, F.
(2017), Forecasting with Temporal Hierarchies, *European Journal of
Operational Research*, 262, 1, 60-74.
[doi:10.1016/j.ejor.2017.02.046](https://doi.org/10.1016/j.ejor.2017.02.046)

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Probabilistic reconciliation:
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
m <- 4 # from quarterly to annual temporal aggregation

# (100 x 14) base forecasts sample matrix (simulated), m = 4, h = 2
sample <- t(sapply(1:100, function(x) {
  rnorm(14, rep(c(20, 10, 5), 2 * c(1, 2, 4)))
}))
# (70 x 1) in-sample residuals vector (simulated)
res <- rnorm(70)

# Top-down probabilistic reconciliation
reco_dist_td <- tesmp(sample[,c(1:2), drop = FALSE], agg_order = m,
                      fun = tetd, weights = c(0.2, 0.5, 0.3, 0.3))

# Middle-out probabilistic reconciliation
reco_dist_mo <- tesmp(sample[,c(3:6), drop = FALSE], agg_order = m,
                      fun = temo, weights = c(0.2, 0.5, 0.3, 0.3), order = 2)

# Bottom-up probabilistic reconciliation
reco_dist_bu <- tesmp(sample[,-c(1:6)], agg_order = m, fun = tebu)

# Level conditional coherent probabilistic reconciliation
reco_dist_lcc <- tesmp(sample, agg_order = m, fun = telcc)

# Optimal cross-sectional probabilistic reconciliation
reco_dist_opt <- tesmp(sample, agg_order = m, res = res, comb = "wlsv")
```
