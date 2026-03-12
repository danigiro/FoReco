# Temporal joint block bootstrap

Joint block bootstrap for generating probabilistic base forecasts that
take into account the correlation between different temporal aggregation
orders (Girolimetto et al. 2023).

## Usage

``` r
teboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL, ...)
```

## Arguments

- model_list:

  A list of all the \\(k^\ast+m)\\ base forecasts models ordered from
  the lowest frequency (most temporally aggregated) to the highest
  frequency. A [`simulate()`](https://rdrr.io/r/stats/simulate.html)
  function for each model has to be available and implemented according
  to the package
  [forecast](https://CRAN.R-project.org/package=forecast), with the
  following mandatory parameters: *object*, *innov*, *future*, and
  *nsim*.

- boot_size:

  The number of bootstrap replicates.

- agg_order:

  Highest available sampling frequency per seasonal cycle (max. order of
  temporal aggregation, \\m\\), or a vector representing a subset of
  \\p\\ factors of \\m\\.

- block_size:

  Block size of the bootstrap, which is typically equivalent to the
  forecast horizon for the most temporally aggregated series.

- seed:

  An integer seed.

- ...:

  Additional arguments for the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function.

## Value

A list with two elements: the seed used to sample the errors and a
(\\\text{boot\\size}\times (k^\ast+m)\text{block\\size}\\) matrix.

## References

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2023), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40(3), 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Bootstrap samples:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md)

Temporal framework:
[`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md),
[`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md),
[`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md),
[`temo()`](https://danigiro.github.io/FoReco/reference/temo.md),
[`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md),
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md),
[`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md),
[`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md),
[`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)
