# Linear combination (aggregation) matrix for a general linearly constrained multiple time series

This function transforms a general (possibly redundant) zero constraints
matrix into a linear combination (aggregation) matrix
\\\mathbf{A}\_{cs}\\. When working with a general linearly constrained
multiple (\\n\\-variate) time series, getting a linear combination
matrix \\\mathbf{A}\_{cs}\\ is a critical step to obtain a
structural-like representation such that \$\$\mathbf{C}\_{cs} =
\[\mathbf{I} \quad -\mathbf{A}\],\$\$ where \\\mathbf{C}\_{cs}\\ is the
full rank zero constraints matrix (Girolimetto and Di Fonzo, 2023).

## Usage

``` r
lcmat(cons_mat, method = "rref", tol = sqrt(.Machine$double.eps),
       verbose = FALSE, sparse = TRUE)
```

## Arguments

- cons_mat:

  A (\\r \times n\\) numeric matrix representing the cross-sectional
  zero constraints.

- method:

  Method to use: "`rref`" for the Reduced Row Echelon Form through
  Gauss-Jordan elimination (*default*), or "`qr`" for the (pivoting) QR
  decomposition (Strang, 2019).

- tol:

  Tolerance for the "`rref`" or "`qr`" method.

- verbose:

  If `TRUE`, intermediate steps are printed (*default* is `FALSE`).

- sparse:

  Option to return a sparse matrix (*default* is `TRUE`).

## Value

A list with two elements: (i) the linear combination (aggregation)
matrix (`agg_mat`) and (ii) the vector of the column permutations
(`pivot`).

## References

Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic
forecast reconciliation for general linearly constrained multiple time
series, *Statistical Methods & Applications*, 33, 581-607.
[doi:10.1007/s10260-023-00738-6](https://doi.org/10.1007/s10260-023-00738-6)
.

Strang, G. (2019), *Linear algebra and learning from data*, Wellesley,
Cambridge Press.

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
## Two hierarchy sharing the same top-level variable, but not sharing the bottom variables
#        X            X
#    |-------|    |-------|
#    A       B    C       D
#  |---|
# A1   A2
# 1) X = C + D,
# 2) X = A + B,
# 3) A = A1 + A2.
cons_mat <- matrix(c(1,-1,-1,0,0,0,0,
               1,0,0,-1,-1,0,0,
               0,0,0,1,0,-1,-1), nrow = 3, byrow = TRUE)
obj <- lcmat(cons_mat = cons_mat, verbose = TRUE)
#> ! A pivot is performed. Remember to apply the pivot also to the base forecast.
#> ℹ E.g. `base[, pivot]` in cross-sectional or `base[pivot, ]` in cross-temporal.
agg_mat <- obj$agg_mat # linear combination matrix
pivot <- obj$pivot # Pivot vector
```
