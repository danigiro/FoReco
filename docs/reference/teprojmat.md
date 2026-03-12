# Projection matrix for optimal combination temporal reconciliation

This function computes the projection or the mapping matrix
\\\mathbf{M}\\ and \\\mathbf{G}\\, respectively, such that
\\\widetilde{\mathbf{y}} = \mathbf{M}\widehat{\mathbf{y}} =
\mathbf{S}\_{te}\mathbf{G}\widehat{\mathbf{y}}\\, where
\\\widetilde{\mathbf{y}}\\ is the vector of the reconciled forecasts,
\\\widehat{\mathbf{y}}\\ is the vector of the base forecasts,
\\\mathbf{S}\_{te}\\ is the temporal structural matrix, and \\\mathbf{M}
= \mathbf{S}\_{te}\mathbf{G}\\. For further information regarding on the
structure of these matrices, refer to Girolimetto et al. (2023).

## Usage

``` r
teprojmat(agg_order, comb = "ols", res = NULL, mat = "M", tew = "sum", ...)
```

## Arguments

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- comb:

  A string specifying the reconciliation method. For a complete list,
  see [tecov](https://danigiro.github.io/FoReco/reference/tecov.md).

- res:

  A (\\N(k^\ast+m) \times 1\\) optional numeric vector containing the
  in-sample residuals or validation errors ordered from the lowest
  frequency to the highest frequency. This vector is used to compute
  come covariance matrices.

- mat:

  A string specifying which matrix to return: "`M`" (*default*) for
  \\\mathbf{M}\\ and "`G`" for \\\mathbf{G}\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- ...:

  Arguments passed on to
  [`tecov`](https://danigiro.github.io/FoReco/reference/tecov.md)

  `mse`

  :   If `TRUE` (*default*) the errors used to compute the covariance
      matrix are not mean-corrected.

  `shrink_fun`

  :   Shrinkage function of the covariance matrix,
      [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
      (*default*)

## Value

The projection matrix \\\mathbf{M}\\ (`mat = "M"`) or the mapping matrix
\\\mathbf{G}\\ (`mat = "G"`).

## References

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

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
[`recoinfo()`](https://danigiro.github.io/FoReco/reference/recoinfo.md),
[`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md),
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
# Temporal framework (annual-quarterly)
Mte <- teprojmat(agg_order = 4, comb = "ols")
Gte <- teprojmat(agg_order = 4, comb = "ols", mat = "G")
```
