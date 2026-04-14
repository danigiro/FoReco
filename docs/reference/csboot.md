# Cross-sectional joint block bootstrap

Joint block bootstrap for generating probabilistic base forecasts that
take into account the correlation between different time series
(Panagiotelis et al. 2023).

## Usage

``` r
csboot(model_list, boot_size, block_size, seed = NULL, xreg = NULL, ...)
```

## Arguments

- model_list:

  A list of all the \\n\\ base forecasts models. A
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function for
  each model has to be available and implemented according to the
  package [forecast](https://CRAN.R-project.org/package=forecast), with
  the following mandatory parameters: *object*, *innov*, *future*, and
  *nsim*.

- boot_size:

  The number of bootstrap replicates.

- block_size:

  Block size of the bootstrap, which is typically equivalent to the
  forecast horizon.

- seed:

  An integer seed.

- xreg:

  An optional 3-d numeric array of dimensions (\\\text{block\\size}
  \times n \times N\_{xreg}\\) containing the new values of `xreg` to be
  used for forecasting. It can contain `NA`s.

- ...:

  Additional arguments for the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function.

## Value

A list with two elements: the seed used to sample the errors and a 3-d
array (\\\text{block\\size} \times n \times \text{boot\\size}\\).

## References

Panagiotelis, A., Gamakumara, P., Athanasopoulos, G. and Hyndman, R.J.
(2023), Probabilistic forecast reconciliation: Properties, evaluation
and score optimisation, *European Journal of Operational Research*,
306(2), 693–706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

## See also

Bootstrap samples:
[`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md),
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md)

Cross-sectional framework:
[`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
[`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
[`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md),
[`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
[`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
[`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
[`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
[`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
[`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)
