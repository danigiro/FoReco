# Cross-sectional bottom-up reconciliation

This function computes the cross-sectional bottom-up reconciled
forecasts (Dunn et al., 1976) for all series by appropriate summation of
the bottom base forecasts \\\widehat{\mathbf{b}}\\:
\$\$\widetilde{\mathbf{y}} = \mathbf{S}\_{cs}\widehat{\mathbf{b}},\$\$
where \\\mathbf{S}\_{cs}\\ is the cross-sectional structural matrix.

## Usage

``` r
csbu(base, agg_mat, sntz = FALSE, round = FALSE)
```

## Arguments

- base:

  A (\\h \times n_b\\) numeric matrix or multivariate time series (`mts`
  class) containing bottom base forecasts; \\h\\ is the forecast
  horizon, and \\n_b\\ is the total number of bottom variables.

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- sntz:

  Logical. If `TRUE`, the negative base forecasts are set to zero (Di
  Fonzo and Girolimetto, 2023) before applying bottom-up. *Default* is
  `FALSE`.

- round:

  Logical. If `TRUE`, base forecasts are rounded before applying the
  bottom-up reconciliation. *Default* is `FALSE`.

## Value

A (\\h \times n\\) numeric matrix of cross-sectional reconciled
forecasts.

## References

Dunn, D. M., Williams, W. H. and Dechaine, T. L. (1976), Aggregate
versus subaggregate models in local area forecasting, *Journal of the
American Statistical Association* 71(353), 68–71.
[doi:10.1080/01621459.1976.10481478](https://doi.org/10.1080/01621459.1976.10481478)

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

## See also

Bottom-up reconciliation:
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md)

Cross-sectional framework:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)

## Examples

``` r
set.seed(123)
# (3 x 2) bottom base forecasts matrix (simulated), Z = X + Y
bts <- matrix(rnorm(6, mean = c(10, 10)), 3, byrow = TRUE)

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
reco <- csbu(base = bts, agg_mat = A)

# Non negative reconciliation
bts[2,2] <- -bts[2,2] # Making negative one of the base forecasts for Y
nnreco <- csbu(base = bts, agg_mat = A, sntz = TRUE)
```
