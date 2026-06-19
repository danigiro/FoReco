# Shrinkage of the Covariance Matrix

Shrinkage of the covariance matrix according to Schäfer and Strimmer
(2005).

## Usage

``` r
shrink_estim(x, mse = TRUE)
```

## Arguments

- x:

  A numeric matrix containing the in-sample residuals or validation
  errors.

- mse:

  If `TRUE` (*default*), the residuals used to compute the covariance
  matrix are not mean-corrected.

## Value

A shrunk covariance matrix.

## References

Schäfer, J.L. and Strimmer, K. (2005), A shrinkage approach to
large-scale covariance matrix estimation and implications for functional
genomics, *Statistical Applications in Genetics and Molecular Biology*,
4, 1.
[doi:10.2202/1544-6115.1175](https://doi.org/10.2202/1544-6115.1175)

## See also

Utilities:
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
[`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md),
[`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md),
[`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md),
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
set.seed(123)
# Simulated in-sample residuals: 50 observations of a 4-variate process
res <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
# Schafer-Strimmer shrunk covariance matrix
shr <- shrink_estim(res)
# Shrinkage intensity selected by the procedure
attr(shr, "lambda")
#> [1] 1
```
