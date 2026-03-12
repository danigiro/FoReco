# Commutation matrix

This function returns the (\\r c \times r c\\) commutation matrix
\\\mathbf{P}\\ such that \\\mathbf{P} \mbox{vec}(\mathbf{Y}) =
\mbox{vec}(\mathbf{Y}'),\\ where \\\mathbf{Y}\\ is a (\\r \times c\\)
matrix (Magnus and Neudecker, 2019).

## Usage

``` r
# Commutation matrix
commat(r, c)

# Vector of indexes for the rows of commutation matrix
commat_index(r, c)
```

## Arguments

- r:

  Number of rows of \\\mathbf{Y}\\.

- c:

  Number of columns of \\\mathbf{Y}\\.

## Value

A sparse (\\r c \times r c\\) matrix \\\mathbf{P}\\ (commat), or the
vector of indexes for the rows of commutation matrix \\\mathbf{P}\\
(commat_index)

## References

Magnus, J.R. and Neudecker, H. (2019), Matrix Differential Calculus with
Applications in Statistics and Econometrics, third edition, New York,
Wiley, pp. 54-55.

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
[`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md),
[`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md),
[`balance_hierarchy()`](https://danigiro.github.io/FoReco/reference/balance_hierarchy.md),
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
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
Y <- matrix(rnorm(30), 5, 6)
P <- commat(5, 6)
P %*% as.vector(Y) == as.vector(t(Y)) # check
#> 30 x 1 Matrix of class "lgeMatrix"
#>       [,1]
#>  [1,] TRUE
#>  [2,] TRUE
#>  [3,] TRUE
#>  [4,] TRUE
#>  [5,] TRUE
#>  [6,] TRUE
#>  [7,] TRUE
#>  [8,] TRUE
#>  [9,] TRUE
#> [10,] TRUE
#> [11,] TRUE
#> [12,] TRUE
#> [13,] TRUE
#> [14,] TRUE
#> [15,] TRUE
#> [16,] TRUE
#> [17,] TRUE
#> [18,] TRUE
#> [19,] TRUE
#> [20,] TRUE
#> [21,] TRUE
#> [22,] TRUE
#> [23,] TRUE
#> [24,] TRUE
#> [25,] TRUE
#> [26,] TRUE
#> [27,] TRUE
#> [28,] TRUE
#> [29,] TRUE
#> [30,] TRUE
```
