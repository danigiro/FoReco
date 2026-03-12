# Cross-temporal bottom-up reconciliation

Cross-temporal bottom-up reconciled forecasts for all series at any
temporal aggregation level are computed by appropriate summation of the
high-frequency bottom base forecasts \\\widehat{\mathbf{B}}^{\[1\]}\\:
\$\$\widetilde{\mathbf{X}} =
\mathbf{S}\_{cs}\widehat{\mathbf{B}}^{\[1\]}\mathbf{S}'\_{te},\$\$ where
\\\mathbf{S}\_{cs}\\ and \\\mathbf{S}\_{te}\\ are the cross-sectional
and temporal structural matrices, respectively.

## Usage

``` r
ctbu(base, agg_mat, agg_order, tew = "sum", sntz = FALSE, round = FALSE)
```

## Arguments

- base:

  A (\\n_b \times hm\\) numeric matrix containing high-frequency bottom
  base forecasts; \\n_b\\ is the total number of high-frequency bottom
  variables, \\m\\ is the max aggregation order, and \\h\\ is the
  forecast horizon for the lowest frequency time series.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

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

A (\\n \times h(k^\ast+m)\\) numeric matrix of cross-temporal reconciled
forecasts.

## References

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

## See also

Bottom-up reconciliation:
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md)

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
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
# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
# (2 x 4) high frequency bottom base forecasts matrix (simulated),
# agg_order = 4 (annual-quarterly)
hfbts <- matrix(rnorm(4*2, 2.5), 2, 4)

reco <- ctbu(base = hfbts, agg_mat = A, agg_order = 4)

# Non negative reconciliation
hfbts[1,4] <- -hfbts[1,4] # Making negative one of the quarterly values for X
nnreco <- ctbu(base = hfbts, agg_mat = A, agg_order = 4, sntz = TRUE)
```
