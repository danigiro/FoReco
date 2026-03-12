# Cross-temporal reconciliation tools

Some useful tools for the cross-temporal forecast reconciliation of a
linearly constrained (e.g., hierarchical/grouped) multiple time series.

## Usage

``` r
cttools(agg_mat, cons_mat, agg_order, tew = "sum", fh = 1, sparse = TRUE)
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

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- tew:

  A string specifying the type of temporal aggregation. Options include:
  "`sum`" (simple summation, *default*), "`avg`" (average), "`first`"
  (first value of the period), and "`last`" (last value of the period).

- fh:

  Forecast horizon for the lowest frequency (most temporally aggregated)
  time series (*default* is `1`).

- sparse:

  Option to return sparse matrices (*default* is `TRUE`).

## Value

A list with four elements:

- dim:

  A vector containing information about the number of series for the
  complete system (`n`), for upper levels (`na`) and bottom level
  (`nb`), the maximum aggregation order (`m`), the number of factor
  (`p`), the partial (`ks`) and total sum (`kt`) of factors.

- set:

  The vector of the temporal aggregation orders (in decreasing order).

- agg_mat:

  The cross-temporal aggregation matrix.

- strc_mat:

  The cross-temporal structural matrix.

- cons_mat:

  The cross-temporal zero constraints matrix.

## See also

Cross-temporal framework:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

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
# Cross-temporal framework
A <- t(c(1,1)) # Aggregation matrix for Z = X + Y
m <- 4 # from quarterly to annual temporal aggregation
cttools(agg_mat = A, agg_order = m)
#> $dim
#>  n na nb  m  p ks kt 
#>  3  1  2  4  3  3  7 
#> 
#> $set
#> [1] 4 2 1
#> 
#> $agg_mat
#> 13 x 8 sparse Matrix of class "dgCMatrix"
#>                      
#>  [1,] 1 1 1 1 1 1 1 1
#>  [2,] 1 1 . . 1 1 . .
#>  [3,] . . 1 1 . . 1 1
#>  [4,] 1 . . . 1 . . .
#>  [5,] . 1 . . . 1 . .
#>  [6,] . . 1 . . . 1 .
#>  [7,] . . . 1 . . . 1
#>  [8,] 1 1 1 1 . . . .
#>  [9,] 1 1 . . . . . .
#> [10,] . . 1 1 . . . .
#> [11,] . . . . 1 1 1 1
#> [12,] . . . . 1 1 . .
#> [13,] . . . . . . 1 1
#> 
#> $strc_mat
#> 21 x 8 sparse Matrix of class "dgCMatrix"
#>                      
#>  [1,] 1 1 1 1 1 1 1 1
#>  [2,] 1 1 . . 1 1 . .
#>  [3,] . . 1 1 . . 1 1
#>  [4,] 1 . . . 1 . . .
#>  [5,] . 1 . . . 1 . .
#>  [6,] . . 1 . . . 1 .
#>  [7,] . . . 1 . . . 1
#>  [8,] 1 1 1 1 . . . .
#>  [9,] 1 1 . . . . . .
#> [10,] . . 1 1 . . . .
#> [11,] 1 . . . . . . .
#> [12,] . 1 . . . . . .
#> [13,] . . 1 . . . . .
#> [14,] . . . 1 . . . .
#> [15,] . . . . 1 1 1 1
#> [16,] . . . . 1 1 . .
#> [17,] . . . . . . 1 1
#> [18,] . . . . 1 . . .
#> [19,] . . . . . 1 . .
#> [20,] . . . . . . 1 .
#> [21,] . . . . . . . 1
#> 
#> $cons_mat
#> 13 x 21 sparse Matrix of class "dgCMatrix"
#>                                                            
#>  [1,] . . .  1  .  .  . . . . -1  .  .  . . . . -1  .  .  .
#>  [2,] . . .  .  1  .  . . . .  . -1  .  . . . .  . -1  .  .
#>  [3,] . . .  .  .  1  . . . .  .  . -1  . . . .  .  . -1  .
#>  [4,] . . .  .  .  .  1 . . .  .  .  . -1 . . .  .  .  . -1
#>  [5,] 1 . . -1 -1 -1 -1 . . .  .  .  .  . . . .  .  .  .  .
#>  [6,] . 1 . -1 -1  .  . . . .  .  .  .  . . . .  .  .  .  .
#>  [7,] . . 1  .  . -1 -1 . . .  .  .  .  . . . .  .  .  .  .
#>  [8,] . . .  .  .  .  . 1 . . -1 -1 -1 -1 . . .  .  .  .  .
#>  [9,] . . .  .  .  .  . . 1 . -1 -1  .  . . . .  .  .  .  .
#> [10,] . . .  .  .  .  . . . 1  .  . -1 -1 . . .  .  .  .  .
#> [11,] . . .  .  .  .  . . . .  .  .  .  . 1 . . -1 -1 -1 -1
#> [12,] . . .  .  .  .  . . . .  .  .  .  . . 1 . -1 -1  .  .
#> [13,] . . .  .  .  .  . . . .  .  .  .  . . . 1  .  . -1 -1
#> 
```
