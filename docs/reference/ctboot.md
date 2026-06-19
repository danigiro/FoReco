# Cross-temporal Joint Block Bootstrap

Joint block bootstrap for generating probabilistic base forecasts that
take into account the correlation between variables at different
temporal aggregation orders (Girolimetto et al. 2023).

## Usage

``` r
ctboot(model_list, boot_size, agg_order, block_size = 1, seed = NULL,
       xreg = NULL, ...)
```

## Arguments

- model_list:

  A list of \\p\\ elements, one for each temporal level ordered from the
  lowest frequency (most temporally aggregated) to the highest
  frequency. Each elements is a list with the \\n\\ base forecasts
  models for each cross-sectional series. A
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function for
  each model has to be available and implemented according to the
  package. [forecast](https://CRAN.R-project.org/package=forecast), with
  the following mandatory parameters: *object*, *innov*, *future*, and
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

- xreg:

  An optional 3-d numeric array of dimensions (\\n \times
  \text{boot\\size}(k^\ast+m) \times N\_{xreg}\\) containing the new
  values of `xreg` to be used for forecasting ordered from the lowest
  frequency to the highest frequency (columns) for each variable (rows).
  It can contain `NA`s.

- ...:

  Additional arguments for the
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) function.

## Value

A list with \\\text{boot\\size}\\ matrix of size
(\\n\times(k^\ast+m)\text{block\\size}\\).

## References

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T. and Hyndman, R.J.
(2024), Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, 40, 3, 1134-1151.
[doi:10.1016/j.ijforecast.2023.10.003](https://doi.org/10.1016/j.ijforecast.2023.10.003)

## See also

Bootstrap samples:
[`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
[`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md)

Cross-temporal framework:
[`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md),
[`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md),
[`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md),
[`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md),
[`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md),
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md),
[`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md),
[`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md),
[`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md),
[`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
[`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)

## Examples

``` r
set.seed(123)

# Minimal example functions: each "model" stores Gaussian residuals;
# simulate() draws new innovations (innov=NULL) or uses the supplied
# ones (innov given).
simple_model <- function(res) {
  structure(list(residuals = res, sigma = sd(res)), class = "simple_model")
}
simulate.simple_model <- function(object, nsim = 1, innov = NULL,
                                  future = TRUE, ...) {
  if (is.null(innov)) {
    rnorm(nsim, mean = 0, sd = object$sigma)
  } else {
    as.numeric(innov)[seq_len(nsim)]
  }
}

# Cross-temporal hierarchy:
# - cross-sectional: Z = X + Y  => n = 3
# - temporal: annual-quarterly  => m = 4, kset = c(4, 2, 1), p = 3 levels
n <- 3
m <- 4
kset <- c(4, 2, 1)
n_obs_per_k <- 40 / kset

# Nested list: outer level = p temporal aggregation orders (low -> high
# frequency), inner level = n cross-sectional series.
model_list <- lapply(seq_along(kset), function(i) {
  lapply(seq_len(n), function(j) {
    simple_model(rnorm(n_obs_per_k[i]))
  })
})

# Joint block bootstrap: 50 replicates, block_size = 1
boot <- ctboot(model_list = model_list, boot_size = 50, agg_order  = m,
               block_size = 1, seed = 1)
#> Error in UseMethod("simulate"): no applicable method for 'simulate' applied to an object of class "simple_model"
```
