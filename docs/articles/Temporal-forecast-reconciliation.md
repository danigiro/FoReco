# Temporal forecast reconciliation

## Introduction

In this vignette, we explore the process of temporal reconciliation
using the `FoReco` package, with a focus on the Capital Country monthly
time series (`ADB`) from the `vndata` dataset ([Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)). Temporal reconciliation is a method
used to ensure that forecasts are coherent across different time
periods, providing a consistent view from high-frequency (e.g., monthly)
to aggregated levels (e.g., yearly). We will first obtain base forecasts
using an ETS model with log transformation and then proceed with several
reconciliation methods: bottom-up, top-down, and optimal combination.
Each method has its own strengths and practical challenges, such as
ensuring non-negativity and incorporating a priori constraints.

### Packages

First, we load the necessary packages.

``` r

library(FoReco)   # -> To perform reconciliation
library(forecast)  # -> To obtain base forecasts
```

## Capital Country monthly time series

In the temporal framework, we are working with a single series. In
particular for this vignette we consider the Capital Country monthly
time series (`ADB`) from the `vndata` dataset ([Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)). See the [dataset
vignette](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.md)
for more details.

``` r

data(vndata)      # dataset
y <- vndata[, "GBDOth"]
te_set <- tetools(12)$set
m <- max(te_set)
```

### Base Forecast

To obtained the base forecasts, we fit an ETS model with log
transformation. We obtain twelve-, six-, four-, three-, two-, and
one-step-ahead base forecasts from the monthly data and the aggregation
over 2, 3, 4, 6, and 12 months.

``` r

data_k <- aggts(y, te_set)
model <- setNames(vector(mode='list', length=length(te_set)), paste0("k-", te_set))
fc_obj <- setNames(vector(mode='list', length=length(te_set)), paste0("k-", te_set))

# ETS model with log transformation
ets_log <- function(x, ...){
  x[x==0] <- min(x[x!=0])/2
  ets(x, lambda = 0, ...)
}

for(k in te_set){
  model[[paste0("k-", k)]] <- ets_log(data_k[[paste0("k-", k)]])
  fc_obj[[paste0("k-", k)]] <- forecast(model[[paste0("k-", k)]], h = m/k)
}
```

We extract the point forecasts and residuals from the fitted models.

``` r

# Point forecasts
base <- lapply(fc_obj, function(x) x$mean)
str(base, give.attr = FALSE)
#> List of 6
#>  $ k-12: Time-Series [1:1] from 2017 to 2017: 1.31
#>  $ k-6 : Time-Series [1:2] from 2017 to 2018: 0.557 0.557
#>  $ k-4 : Time-Series [1:3] from 2017 to 2018: 0.356 0.356 0.356
#>  $ k-3 : Time-Series [1:4] from 2017 to 2018: 0.243 0.243 0.243 0.243
#>  $ k-2 : Time-Series [1:6] from 2017 to 2018: 0.225 0.225 0.225 0.225 0.225 ...
#>  $ k-1 : Time-Series [1:12] from 2017 to 2018: 0.178 0.178 0.178 0.178 0.178 ...
# Residuals
res <- lapply(fc_obj, residuals, type = "response")
str(res, give.attr = FALSE)
#> List of 6
#>  $ k-12: Time-Series [1:19] from 1998 to 2016: 18.15 -1.18 1.13 31.64 13.37 ...
#>  $ k-6 : Time-Series [1:38] from 1998 to 2016: -0.422 18.911 -0.422 -0.422 -0.422 ...
#>  $ k-4 : Time-Series [1:57] from 1998 to 2017: -0.221 14.333 4.423 -0.221 -0.221 ...
#>  $ k-3 : Time-Series [1:76] from 1998 to 2017: -0.243 -0.228 14.338 4.377 -0.306 ...
#>  $ k-2 : Time-Series [1:114] from 1998 to 2017: -0.0904 -0.0904 -0.0904 14.4633 -0.0905 ...
#>  $ k-1 : Time-Series [1:228] from 1998 to 2017: -0.0433 -0.0433 -0.0433 -0.0433 -0.0433 ...
```

### Point Reconciliation

Bottom-up reconciliation ([Dunn et al., 1976](#ref-Dunn1976-kv))
aggregates the high-frequency forecasts from the monthly level to higher
temporal levels ([Girolimetto et al., 2024](#ref-Girolimetto2023-jm)).

``` r

fc_bts <- base$`k-1`
rf_bu <- tebu(fc_bts, agg_order = m)
summary(rf_bu)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `tebu`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#> k-12 h-1  k-6 h-1  k-6 h-2  k-4 h-1 
#> 2.137377 1.068688 1.068688 0.712459 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

To obtain a list of forecasts at different orders of aggregation, we can
use the `components` function.

``` r

str(components(rf_bu), give.attr=FALSE)
#> List of 6
#>  $ k-12: num 2.14
#>  $ k-6 : num [1:2] 1.07 1.07
#>  $ k-4 : num [1:3] 0.712 0.712 0.712
#>  $ k-3 : num [1:4] 0.534 0.534 0.534 0.534
#>  $ k-2 : num [1:6] 0.356 0.356 0.356 0.356 0.356 ...
#>  $ k-1 : num [1:12] 0.178 0.178 0.178 0.178 0.178 ...
```

In top-down reconciliation for hierarchical time series, the forecast
for the annual series (\\k = 12\\) is distributed proportionally to
ensure the top-level value stays the same and all bottom-level forecasts
are non-negative ([Gross & Sohl, 1990](#ref-Gross1990-uf)).

``` r

y_mat <- matrix(data_k$`k-1`, ncol = m, byrow = TRUE)
x12 <- data_k$`k-12`
fc_x12 <- base$`k-12`

# Average historical proportions - Gross-Sohl method A
p_gsa <- colMeans(apply(y_mat, 2, function(x) x/x12), na.rm = TRUE)
rf_td_gsa <- tetd(fc_x12, agg_order = m, weights = p_gsa)
summary(rf_td_gsa)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `tetd`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.3142157 0.4965682 0.8176475 0.2156414 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.

# Proportions of the historical averages - Gross-Sohl method F
p_gsf <- colMeans(y_mat)/mean(x12)
rf_td_gsf <- tetd(fc_x12, agg_order = m, weights = p_gsf)
summary(rf_td_gsf)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `tetd`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.3142157 0.5007542 0.8134615 0.0790347 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

To perform temporal reconciliation with `FoReco` using the base
forecasts at any temporal aggregation order, it is necessary to arrange
base forecasts (and residuals) in vector form.

``` r

base_vec <- unlist(base, use.names = FALSE)
res_vec <- unlist(res, use.names = FALSE)
```

The level conditional coherent reconciliation (LCC) is a generalization
of the original proposal by Hollyman et al.
([2021](#ref-Hollyman2021-zq)) and Di Fonzo & Girolimetto
([2024](#ref-DiFonzo2024-ijf)) to include the temporal framework.

``` r

rf_lcc <- telcc(base = base_vec, agg_order = m,
                res = res_vec, comb = "wlsv")
summary(rf_lcc)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `telcc`
#> • Covariance approximation approach: `wlsv`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.3258648 0.6629324 0.6629324 0.4419549 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

Finally we can obtained the optimal (in least squares sense) combination
temporal reconciled forecast ([Athanasopoulos et al.,
2017](#ref-Athanasopoulos2017-zh)).

``` r

rf_opt <- terec(base = base_vec, agg_order = m,
                res = res_vec, comb = "sar1")
summary(rf_opt)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `sar1`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.3291002 0.6645501 0.6645501 0.4432082 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

The following table shows some options for the optimal combination
temporal reconciliation function
[`terec()`](https://danigiro.github.io/FoReco/reference/terec.md).

[TABLE]

### Practical challenges

#### Non negativity issues

It is not always true that reconciled forecasts will remain non-negative
even if the base forecasts are non-negative. For example, consider the
case of an shrinkage covariance matrix (`"shr"`).

``` r

rf_opt_shr <- terec(base = base_vec, agg_order = m,
                    res = res_vec, comb = "shr")
summary(rf_opt_shr)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `shr`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 0.9952037 0.5299384 0.4652653 0.4565095 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

To address this issue, we can use two approaches:

- State-of-the-art numerical optimization procedure, **osqp** ([Stellato
  et al., 2020](#ref-Stellato2020-jd)).

``` r

rf_osqp <- terec(base = base_vec, agg_order = m,
                res = res_vec, comb = "shr", nn = "osqp")
summary(rf_osqp)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `shr`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>     obj_val    run_time iter     prim_res status status_polish
#> 1 -16.58821 0.000254749   50 5.551115e-17      1             1
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.0625336 0.5934914 0.4690422 0.4557664 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

- Simple heuristic strategy: set-negative-to-zero, **sntz** ([Di Fonzo &
  Girolimetto, 2023](#ref-Di_Fonzo2023-ae)).

``` r

rf_sntz <- terec(base = base_vec, agg_order = m,
                res = res_vec, comb = "shr", nn = "sntz")
summary(rf_sntz)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `shr`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 2.861023e-06 1.490116e-08    0
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.0632623 0.5954494 0.4678130 0.4565095 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

#### A priori constrained (immutable) forecasts

Sometimes we may wish to incorporate a priori knowledge during the
reconciliation process ([Zhang et al., 2023](#ref-Zhang2023-yj)) in
order to improve the accuracy of the reconciled forecasts. For example,
suppose we want to fix the forecasts of the annual level at the base
forecasts values.

``` r

rf_imm <- terec(base = base_vec, agg_order = m,
                res = res_vec, comb = "sar1", immutable = rbind(c(12,1)))
summary(rf_imm)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `sar1`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1 
#> 1.3142157 0.6571078 0.6571078 0.4382244 
#> ... (24 more elements)
#> Use `print(x, n_row)` to see more elements.
```

``` r

rf_imm[1] - base_vec[1]
#> k-12 h-1 
#>        0
```

#### Exploring a subset of temporal aggregation orders

Our approaches so far have involved considering all factors of \\m\\ as
potential aggregation orders. Nevertheless, it is worth noting that we
could also focus only on a given subset of these factors. For example,
we could be interested only on monthly, quarterly and annual forecasts.

``` r

te_subset <- c(12, 3, 1)
base_vec2 <- unlist(base[paste0("k-", te_subset)], use.names = FALSE)
res_vec2 <- unlist(res[paste0("k-", te_subset)], use.names = FALSE)
rf_sub <- terec(base = base_vec2, agg_order = te_subset,
                res = res_vec2, comb = "sar1")
summary(rf_sub)
#> ✔ Temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `terec`
#> • Covariance approximation approach: `sar1`
#> • Output: unknown object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 3, and 1
#> • Forecast horizons (h) per k: 1, 4, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>  k-12 h-1   k-3 h-1   k-3 h-2   k-3 h-3 
#> 1.4803010 0.3710480 0.3691025 0.3691025 
#> ... (13 more elements)
#> Use `print(x, n_row)` to see more elements.
```

### Probabilistic Reconciliation

Following Panagiotelis et al. ([2023](#ref-Panagiotelis2023-se)) and
Girolimetto et al. ([2024](#ref-Girolimetto2023-jm)), the reconciliation
of probabilistic forecasts is a two-step process: first, we sample from
an incoherent distributio, and then we reconcile each sample.

We can use a non-parametric method, the joint block bootstrap to
simulate \\B\\ samples and then reconciled them.

``` r

# Base forecasts' sample:
# we simulate from the base models by sampling errors 
# while keeping the cross-sectional dimension fixed.
B <- 100
base_tejb <- teboot(model, B, m)
dim(base_tejb)
#> [1] 100  28
reco_tejb <- tesmp(base_tejb, agg_order = m, res = res_vec, nn = "sntz", comb = "sar1")
summary(reco_tejb)
#> ✔ Temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `tesmp(terec)`
#> • Covariance approximation approach: `sar1`
#> • Output: distributional object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>        run_time          tol iter
#> 11 6.914139e-06 1.490116e-08    0
#> 12 6.914139e-06 1.490116e-08    0
#> 13 6.914139e-06 1.490116e-08    0
#> 21 6.914139e-06 1.490116e-08    0
#> 22 6.914139e-06 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>       tao-1 
#> sample[100]

# Extracts mean:
str(as_tevector(mean(reco_tejb), agg_order = m), give.attr = FALSE)
#>  Named num [1:28] 5.968 2.046 3.922 0.548 4.145 ...
```

A parametric method assumes a normal distribution (Gaussian) to generate
the incoherent sample set of forecasts.

``` r

# Gaussian reconciled distribution
reco_teg_sar1 <- temvn(base = base_vec, agg_order = m, comb = "sar1", res = res_vec)
summary(reco_teg_sar1)
#> ✔ Temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `temvn`
#> • Covariance approximation approach: `sar1`
#> • Output: distributional object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 2, 4, 6, 8, 12, and 24
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>   tao-1 
#> MVN[28]

# Gaussian reconciled distribution with different base covariance matrix
# Multi-step residuals
hres <- lapply(model, function(fit)
  lapply(1:frequency(fit$x), function(h) residuals(fit, type='response', h = h)))
hres <- Reduce("c", lapply(hres, arrange_hres))
# Re-arrenge multi-step residuals in a matrix form
mres <- as_hstack_telayout(hres, agg_order = m)
cov_base <- shrink_estim(mres)

reco_teg <- temvn(base = base_vec, agg_order = m, comb = "sar1", 
                    res = res_vec, comb_base = cov_base)
summary(reco_teg) 
#> ✔ Temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `temvn`
#> • Covariance approximation approach: `sar1`
#> • Output: distributional object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 2, 4, 6, 8, 12, and 24
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>   tao-1 
#> MVN[28]

# Reconciled sample distribution starting from gaussian base forecasts
base_tes <- MASS::mvrnorm(n = B, mu = base_vec, Sigma = shrink_estim(mres))
reco_tes <- tesmp(base_tes, agg_order = m, res = res_vec, nn = "sntz", comb = "sar1")
summary(reco_tes)
#> ✔ Temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `tesmp(terec)`
#> • Covariance approximation approach: `sar1`
#> • Output: distributional object
#> 
#> ── Structure
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 1.001358e-05 1.490116e-08    0
#> 2 1.001358e-05 1.490116e-08    0
#> 3 1.001358e-05 1.490116e-08    0
#> 4 1.001358e-05 1.490116e-08    0
#> 5 1.001358e-05 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>       tao-1 
#> sample[100]
```

## References

Athanasopoulos, G., Hyndman, R. J., Kourentzes, N., & Petropoulos, F.
(2017). Forecasting with temporal hierarchies. *European Journal of
Operational Research*, *262*(1), 60–74.
<https://doi.org/10.1016/j.ejor.2017.02.046>

Di Fonzo, T., & Girolimetto, D. (2023). Spatio-temporal reconciliation
of solar forecasts. *Solar Energy*, *251*, 13–29.
<https://doi.org/10.1016/j.solener.2023.01.003>

Di Fonzo, T., & Girolimetto, D. (2024). Forecast combination-based
forecast reconciliation: Insights and extensions. *International Journal
of Forecasting*, *40*(2), 490–514.
<https://doi.org/10.1016/j.ijforecast.2022.07.001>

Dunn, D. M., Williams, W. H., & Dechaine, T. L. (1976). Aggregate versus
subaggregate models in local area forecasting. *Journal of the American
Statistical Association*, *71*(353), 68–71.
<https://doi.org/10.1080/01621459.1976.10481478>

Girolimetto, D., Athanasopoulos, G., Di Fonzo, T., & Hyndman, R. J.
(2024). Cross-temporal probabilistic forecast reconciliation:
Methodological and practical issues. *International Journal of
Forecasting*, *40*(3), 1134–1151.
<https://doi.org/10.1016/j.ijforecast.2023.10.003>

Gross, C. W., & Sohl, J. E. (1990). Disaggregation methods to expedite
product line forecasting. *Journal of Forecasting*, *9*(3), 233–254.
<https://doi.org/10.1002/for.3980090304>

Hollyman, R., Petropoulos, F., & Tipping, M. E. (2021). Understanding
forecast reconciliation. *European Journal of Operational Research*,
*294*(1), 149–160. <https://doi.org/10.1016/j.ejor.2021.01.017>

Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., & Hyndman, R. J.
(2023). Probabilistic forecast reconciliation: Properties, evaluation
and score optimisation. *European Journal of Operational Research*,
*306*(2), 693–706. <https://doi.org/10.1016/j.ejor.2022.07.040>

Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020).
OSQP: An operator splitting solver for quadratic programs. *Mathematical
Programming Computation*, *12*(4), 637–672.
<https://doi.org/10.1007/s12532-020-00179-2>

Wickramasuriya, S. L., Athanasopoulos, G., & Hyndman, R. J. (2018).
Optimal Forecast Reconciliation for Hierarchical and Grouped Time Series
Through Trace Minimization. *Journal of the American Statistical
Association*, *114*(526), 804–819.
<https://doi.org/10.1080/01621459.2018.1448825>

Zhang, B., Kang, Y., Panagiotelis, A., & Li, F. (2023). Optimal
reconciliation with immutable forecasts. *European Journal of
Operational Research*, *308*(2), 650–660.
<https://doi.org/10.1016/j.ejor.2022.11.035>
