# Shrinkage of the covariance matrix using the Oracle approximation

Shrinkage of the covariance matrix according to the Oracle Approximating
Shrinkage (OAS) of Chen et al. (2009) and Ando and Xiao (2023).

## Usage

``` r
shrink_oasd(x, mse = TRUE)
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

Ando, S., and Xiao, M. (2023), High-dimensional covariance matrix
estimation: shrinkage toward a diagonal target. *IMF Working Papers*,
2023(257), A001.

Chen, Y., Wiesel, A., and Hero, A. O. (2009), Shrinkage estimation of
high dimensional covariance matrices, *2009 IEEE international
conference on acoustics, speech and signal processing*, 2937–2940. IEEE.

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
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)
