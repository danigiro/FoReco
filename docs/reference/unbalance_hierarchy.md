# Aggregation matrix of a balanced hierarchy in (possibly) unbalanced form

A hierarchy with \\L\\ upper levels is said to be balanced if each
variable at level \\l\\ has at least one child at level \\l+1\\. When
this doesn't hold, the hierarchy is unbalanced. This function transforms
an aggregation matrix of a balanced hierarchy into an aggregation matrix
of an unbalanced one, by removing possible duplicated series.

## Usage

``` r
unbalance_hierarchy(agg_mat, more_info = FALSE, sparse = TRUE)
```

## Arguments

- agg_mat:

  A (\\n_a \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix. It maps the \\n_b\\ bottom-level (free) variables
  into the \\n_a\\ upper (constrained) variables.

- more_info:

  If `TRUE`, it returns only the aggregation matrix of the unbalanced
  hierarchy. *Default* is `FALSE`.

- sparse:

  Option to return sparse matrices (*default* is `TRUE`).

## Value

A list containing four elements (`more_info = TRUE`):

- ubm:

  The aggregation matrix of the unbalanced hierarchy.

- agg_mat:

  The input matrix.

- idrm:

  The identification number of the duplicated variables (row numbers of
  the aggregation matrix `agg_mat`).

- id:

  The identification number of each variable in the balanced hierarchy.
  It may contains duplicated values.

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
[`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)

## Examples

``` r
#     Balanced     ->     Unbalanced
#        T                    T
#    |-------|            |-------|
#    A       B            A       |
#  |---|     |          |---|     |
# AA   AB    BA        AA   AB    BA
A <- matrix(c(1, 1, 1,
              1, 1, 0,
              0, 0, 1), 3, byrow = TRUE)
obj <- unbalance_hierarchy(agg_mat = A)
obj
#> 2 x 3 sparse Matrix of class "dgCMatrix"
#>           
#> [1,] 1 1 1
#> [2,] 1 1 .
```
