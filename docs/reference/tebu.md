# Temporal bottom-up reconciliation

Temporal bottom-up reconciled forecasts at any temporal aggregation
level are computed by appropriate aggregation of the high-frequency base
forecasts, \\\widehat{\mathbf{x}}^{\[1\]}\\: \$\$\widetilde{\mathbf{x}}
= \mathbf{S}\_{te}\widehat{\mathbf{x}}^{\[1\]},\$\$ where
\\\mathbf{S}\_{te}\\ is the temporal structural matrix.

## Usage

``` r
tebu(base, agg_order, tew = "sum", sntz = FALSE, round = FALSE)
```

## Arguments

- base:

  A (\\hm \times 1\\) numeric vector containing the high-frequency base
  forecasts; \\m\\ is the max. temporal aggregation order, and \\h\\ is
  the forecast horizon for the lowest frequency time series.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- sntz:

  Logical. If `TRUE`, the negative base forecasts are set to zero (Di
  Fonzo and Girolimetto, 2023) before applying bottom-up. *Default* is
  `FALSE`.

- round:

  Logical. If `TRUE`, base forecasts are rounded before applying the
  bottom-up reconciliation. *Default* is `FALSE`.

## Value

A (\\h(k^\ast+m) \times 1\\) numeric vector of temporal reconciled
forecasts.

## References

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

## See also

Bottom-up reconciliation:
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md)

Temporal framework:
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
set.seed(123)
# (4 x 1) high frequency base forecasts vector (simulated),
# agg_order = 4 (annual-quarterly)
hfts <- rnorm(4, 5)

reco <- tebu(base = hfts, agg_order = 4)

# Non negative reconciliation
hfts[4] <- -hfts[4] # Making negative one of the quarterly base forecasts
nnreco <- tebu(base = hfts, agg_order = 4, sntz = TRUE)
```
