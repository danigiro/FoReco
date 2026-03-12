# Projection matrix for optimal combination cross-sectional reconciliation

This function computes the projection or the mapping matrix
\\\mathbf{M}\\ and \\\mathbf{G}\\, respectively, such that
\\\widetilde{\mathbf{y}} = \mathbf{M}\widehat{\mathbf{y}} =
\mathbf{S}\_{cs}\mathbf{G}\widehat{\mathbf{y}}\\, where
\\\widetilde{\mathbf{y}}\\ is the vector of the reconciled forecasts,
\\\widehat{\mathbf{y}}\\ is the vector of the base forecasts,
\\\mathbf{S}\_{cs}\\ is the cross-sectional structural matrix, and
\\\mathbf{M} = \mathbf{S}\_{cs}\mathbf{G}\\. For further information
regarding on the structure of these matrices, refer to Girolimetto et
al. (2023).

## Usage

``` r
csprojmat(agg_mat, cons_mat, comb = "ols", res = NULL, mat = "M", ...)
```

## Arguments

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- cons_mat:

  A (\\n_a \times n\\) numeric matrix representing the cross-sectional
  zero constraints: each row represents a constraint equation, and each
  column represents a variable. The matrix can be of full rank, meaning
  the rows are linearly independent, but this is not a strict
  requirement, as the function allows for redundancy in the constraints.

- comb:

  A string specifying the reconciliation method. For a complete list,
  see [cscov](https://danigiro.github.io/FoReco/reference/cscov.md).

- res:

  An (\\N \times n\\) optional numeric matrix containing the in-sample
  residuals or validation errors. This matrix is used to compute some
  covariance matrices.

- mat:

  A string specifying which matrix to return: "`M`" (*default*) for
  \\\mathbf{M}\\ and "`G`" for \\\mathbf{G}\\.

- ...:

  Arguments passed on to
  [`cscov`](https://danigiro.github.io/FoReco/reference/cscov.md)

  `mse`

  :   If `TRUE` (*default*) the errors used to compute the covariance
      matrix are not mean-corrected.

  `shrink_fun`

  :   Shrinkage function of the covariance matrix,
      [shrink_estim](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
      (*default*).

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
# Cross-sectional framework
A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
Mcs <- csprojmat(agg_mat = A, comb = "ols")
Gcs <- csprojmat(agg_mat = A, comb = "ols", mat = "G")
```
