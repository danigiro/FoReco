# Aggregation matrix of a (possibly) unbalanced hierarchy in balanced form

A hierarchy with \\L\\ upper levels is said to be balanced if each
variable at level \\l\\ has at least one child at level \\l+1\\. When
this doesn't hold, the hierarchy is unbalanced. This function transforms
an aggregation matrix of an unbalanced hierarchy into an aggregation
matrix of a balanced one. This function is used to reconcile forecasts
with [cslcc](https://danigiro.github.io/FoReco/reference/cslcc.md),
which operates exclusively with balanced hierarchies.

## Usage

``` r
balance_hierarchy(agg_mat, nodes = "auto", sparse = TRUE)
```

## Arguments

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- nodes:

  A (\\L \times 1\\) numeric vector indicating the number of variables
  in each of the upper \\L\\ levels of the hierarchy. The *default*
  value is the string "`auto`" which calculates the number of variables
  in each level.

- sparse:

  Option to return sparse matrices (*default* is `TRUE`).

## Value

A list containing four elements:

- bam:

  The balanced aggregation matrix.

- agg_mat:

  The input matrix.

- nodes:

  A (\\L \times 1\\) numeric vector indicating the number of variables
  in each of the \\L\\ upper levels of the balanced hierarchy.

- id:

  The identification number of each variable in the balanced hierarchy.
  It may contains duplicated values.

## See also

Utilities:
[`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md),
[`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md),
[`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md),
[`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md),
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
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md),
[`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)

## Examples

``` r
#    Unbalanced     ->      Balanced
#        T                     T
#    |-------|             |-------|
#    A       |             A       B
#  |---|     |           |---|     |
# AA   AB    B          AA   AB    BA
A <- matrix(c(1, 1, 1,
              1, 1, 0), 2, byrow = TRUE)
obj <- balance_hierarchy(agg_mat = A, nodes = c(1, 1))
obj$bam
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>           
#> [1,] 1 1 1
#> [2,] 1 1 .
#> [3,] . . 1
```
