# Cross-sectional aggregation matrix of a dataframe

This function allows the user to easily build the (\\n_a \times n_b\\)
cross-sectional aggregation matrix starting from a data frame.

## Usage

``` r
df2aggmat(formula, data, sep = "_", sparse = TRUE, top_label = "Total",
          verbose = TRUE)
```

## Arguments

- formula:

  Specification of the hierarchical structure: grouped hierarchies are
  specified using `~ g1 * g2` and nested hierarchies are specified using
  `~ parent / child`. Mixtures of the two formulations are also
  possible, like `~ g1 * (grandparent / parent / child)`.

- data:

  A dataset in which each column contains the values of the variables in
  the formula and each row identifies a bottom level time series.

- sep:

  Character to separate the names of the aggregated series, (*default*
  is "`_`").

- sparse:

  Option to return sparse matrices (*default* is `TRUE`).

- top_label:

  Label of the top level variable (*default* is "`Total`").

- verbose:

  If `TRUE` (*default*), hierarchy informations are printed.

## Value

A (`na x nb`) matrix.

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
## Balanced hierarchy
#         T
#    |--------|
#    A        B
#  |---|   |--|--|
# AA   AB  BA BB BC
# Names of the bottom level variables
data_bts <- data.frame(X1 = c("A", "A", "B", "B", "B"),
                       X2 = c("A", "B", "A", "B", "C"),
                       stringsAsFactors = FALSE)
# Cross-sectional aggregation matrix
agg_mat <- df2aggmat(~ X1 / X2, data_bts, sep = "", verbose = FALSE)

## Unbalanced hierarchy
#                 T
#       |---------|---------|
#       A         B         C
#     |---|     |---|     |---|
#    AA   AB   BA   BB   CA   CB
#  |----|         |----|
# AAA  AAB       BBA  BBB
# Names of the bottom level variables
data_bts <- data.frame(X1 = c("A", "A", "A", "B", "B", "B", "C", "C"),
                       X2 = c("A", "A", "B", "A", "B", "B", "A", "B"),
                       X3 = c("A", "B", NA, NA, "A", "B", NA, NA),
                       stringsAsFactors = FALSE)
# Cross-sectional aggregation matrix
agg_mat <- df2aggmat(~ X1 / X2 / X3, data_bts, sep = "", verbose = FALSE)

## Group of two hierarchies
#     T          T         T | A  | B  | C
#  |--|--|  X  |---|  ->  ---+----+----+----
#  A  B  C     M   F       M | AM | BM | CM
#                          F | AF | BF | CF
# Names of the bottom level variables
data_bts <- data.frame(X1 = c("A", "A", "B", "B", "C", "C"),
                       Y1 = c("M", "F", "M", "F", "M", "F"),
                       stringsAsFactors = FALSE)
# Cross-sectional aggregation matrix
agg_mat <- df2aggmat(~ Y1 * X1, data_bts, sep = "", verbose = FALSE)
```
