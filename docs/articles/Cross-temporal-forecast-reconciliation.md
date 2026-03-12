# Cross-temporal forecast reconciliation

## Introduction

This vignette demonstrates the process of cross-temporal forecast
reconciliation using the `FoReco` package. The vignette covers the
following steps:

1.  Preparing and loading the necessary packages and data.
2.  Generating base forecasts for grouped time series.
3.  Reconciling point forecasts.
4.  Addressing practical challenges such as non-negativity issues.
5.  Exploring probabilistic forecast reconciliation.

### Packages

First, we load the necessary packages.

``` r

library(FoReco)   # -> To perform reconciliation
library(forecast)  # -> To obtain base forecasts
library(GMCM)      # -> Sample from a multivariate normal distibution
```

## `vndata`: Groupped monthly time series

We will use the `vndata` dataset ([Wickramasuriya et al.,
2018](#ref-Wickramasuriya2019)), which contains grouped monthly time
series data, and `vnaggmat`, which is the corresponding aggregation
matrix. See the [dataset
vignette](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.md)
for more details.

``` r

data(vndata)      # dataset
data(vnaggmat)    # Agg mat matrix
```

### Base forecast

To obtained the base forecasts, we fit an ETS model with log
transformation to each series. We handle zeros by replacing them with
half the minimum non-zero value in the series ([Wickramasuriya et al.,
2020](#ref-Wickramasuriya2020-zk)), then fit the ETS model and generate
forecasts. We obtain twelve-, six-, four-, three-, two-, and
one-step-ahead base forecasts from the monthly data and the aggregation
over 2, 3, 4, 6, and 12 months.

``` r

te_set <- tetools(12)$set
m <- max(te_set)
data_k <- aggts(vndata, te_set)
model <- setNames(rep(list(setNames(vector(mode='list', length=NCOL(vndata)), colnames(vndata))), 
                      length(te_set)), paste0("k-", te_set))
fc_obj <- setNames(rep(list(setNames(vector(mode='list', length=NCOL(vndata)), colnames(vndata))), 
                      length(te_set)), paste0("k-", te_set))

# ETS model with log transformation
ets_log <- function(x, ...){
  x[x==0] <- min(x[x!=0])/2
  ets(x, lambda = 0, ...)
}

for(k in te_set){
  idk <- paste0("k-", k)
  for(i in 1:NCOL(vndata)){
    ids <- colnames(vndata)[i]
    model[[idk]][[ids]] <- ets_log(data_k[[idk]][,i])
    fc_obj[[idk]][[ids]] <- forecast(model[[idk]][[ids]], h = m/k)
  }
  cat(k, " ")
}
#> 12  6  4  3  2  1
```

We extract the point forecasts and residuals from the fitted models.

``` r

# Point forecasts
base <- lapply(fc_obj, function(x) rbind(sapply(x, function(y) y$mean)))
str(base, give.attr=FALSE)
#> List of 6
#>  $ k-12: num [1, 1:525] 327178 99476 63854 79348 20809 ...
#>  $ k-6 : num [1:2, 1:525] 171018 157987 52963 46291 35704 ...
#>  $ k-4 : num [1:3, 1:525] 126347 98437 105568 41869 29839 ...
#>  $ k-3 : num [1:4, 1:525] 96845 75196 79409 79498 30443 ...
#>  $ k-2 : num [1:6, 1:525] 72593 54242 45402 53122 55872 ...
#>  $ k-1 : num [1:12, 1:525] 50651 21336 24567 29800 22846 ...

# Residuals
res <- lapply(fc_obj, function(x) rbind(sapply(x, residuals, type = "response")))
str(res, give.attr=FALSE)
#> List of 6
#>  $ k-12: num [1:19, 1:525] 3241 3444 -1935 -4632 8946 ...
#>  $ k-6 : num [1:38, 1:525] -140.6 -408 4588.6 -4072.2 -71.5 ...
#>  $ k-4 : num [1:57, 1:525] -20.7 19.8 -1219.6 4592.3 -706.1 ...
#>  $ k-3 : num [1:76, 1:525] 192 -542 674 -1212 3679 ...
#>  $ k-2 : num [1:114, 1:525] -329 -720 -453 538 -693 ...
#>  $ k-1 : num [1:228, 1:525] 2143 -970 -115 133 951 ...
```

### Point forecast reconciliation

Within `FoReco`, a range of reconciliation strategies are available,
including bottom-up, top-down, level conditional coherent forecast
reconciliation, and cross-temporal heuristics.

Bottom-up reconciliation ([Dunn et al., 1976](#ref-Dunn1976-kv))
aggregates the high-frequency forecasts from the lowest cross-sectional
level to higher cross-temporal levels ([Girolimetto et al.,
2024](#ref-Girolimetto2023-jm)).

``` r

fc_bts <- t(base$`k-1`[, colnames(vnaggmat)])
rf_bu <- ctbu(fc_bts, agg_order = m, agg_mat = vnaggmat)
str(rf_bu, give.attr=FALSE)
#>  num [1:525, 1:28] 280206 89744 55472 69366 16847 ...
```

To obtain a list of forecasts at different orders of aggregation, we can
use the `FoReco2matrix` function.

``` r

str(FoReco2matrix(rf_bu), give.attr=FALSE)
#> List of 6
#>  $ k-12: num [1, 1:525] 280206 89744 55472 69366 16847 ...
#>  $ k-6 : num [1:2, 1:525] 147023 133182 48055 41689 31372 ...
#>  $ k-4 : num [1:3, 1:525] 108461 82962 88782 36173 25354 ...
#>  $ k-3 : num [1:4, 1:525] 82898 64126 66292 66891 27564 ...
#>  $ k-2 : num [1:6, 1:525] 62313 46149 38562 44400 46889 ...
#>  $ k-1 : num [1:12, 1:525] 44094 18218 20585 25563 19385 ...
```

In top-down reconciliation for hierarchical time series, the forecast
for the top-level series (Total) is distributed proportionally to ensure
the top-level value stays the same and all bottom-level forecasts are
non-negative ([Gross & Sohl, 1990](#ref-Gross1990-uf)).

``` r

bts_mat <- data_k$`k-1`[, colnames(vnaggmat)]
tot_12 <- data_k$`k-12`[,1]
fc_tot_12 <- base$`k-12`[,1]

# Average historical proportions - Gross-Sohl method A
p_gsa <- apply(bts_mat, 2, function(y){
  colMeans(apply(matrix(y, ncol = m, byrow = TRUE), 2, function(x) x/tot_12))
})
rf_td_gsa <- cttd(fc_tot_12, agg_order = m, weights = t(p_gsa), agg_mat = vnaggmat)
str(rf_td_gsa, give.attr=FALSE)
#>  num [1:525, 1:28] 327178 105757 63173 85397 22176 ...

# Proportions of the historical averages - Gross-Sohl method F
p_gsf <- apply(bts_mat, 2, function(y){
  colMeans(matrix(y, ncol = m, byrow = TRUE))/mean(tot_12)
})
rf_td_gsf <- cttd(fc_tot_12, agg_order = m, weights = t(p_gsf), agg_mat = vnaggmat)
str(rf_td_gsf, give.attr=FALSE)
#>  num [1:525, 1:28] 327178 105656 63188 85290 22160 ...
```

To perform cross-temporal reconciliation with FoReco using the complete
set of base forecasts (any cross-sectional and temporal level), it is
necessary to arrange base forecasts (and residuals) in matrix form. The
rows of the matrix represent the cross-sectional variables, while the
columns the temporal dimension.

``` r

base_mat <- t(Reduce(rbind, base))
res_mat <- t(Reduce(rbind, res))
```

The level conditional coherent reconciliation (LCC) is a generalization
of the original proposal by Hollyman et al.
([2021](#ref-Hollyman2021-zq)) and Di Fonzo & Girolimetto
([2024](#ref-DiFonzo2024-ijf)) to include the cross-temporal framework

``` r

rf_lcc <- ctlcc(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                res = res_mat, comb = "wlsv")
str(rf_lcc, give.attr=FALSE)
#>  num [1:525, 1:28] 316724 98441 62491 79421 19591 ...
```

The iterative procedure described in Di Fonzo & Girolimetto
([2023a](#ref-Di_Fonzo2023-dg)) produces cross-temporally reconciled
forecasts by alternating forecast reconciliation along one single
dimension (either cross-sectional or temporal) at each iteration step.

``` r

rf_ite <- iterec(base = base_mat, res = res_mat,
                 cslist = list(agg_mat = vnaggmat, comb = "shr"), 
                 telist = list(agg_order = m, comb = "wlsv"))
#> ── Iterative heuristic cross-temporal forecast reconciliation ──────────────────
#> Legend: i = iteration; s = step. Norm = "inf".
#> 
#>   i.s |        Temporal | Cross-sectional |
#>     0 |         4869.77 |        23213.13 |
#>   1.1 |            0.00 |        34261.97 |
#>   1.2 |         6994.04 |            0.00 |
#>   2.1 |            0.00 |         1112.09 |
#>   2.2 |          342.64 |            0.00 |
#>   3.1 |            0.00 |           47.06 |
#>   3.2 |           16.48 |            0.00 |
#>   4.1 |            0.00 |            2.15 |
#>   4.2 |        7.90e-01 |            0.00 |
#>   5.1 |            0.00 |        1.01e-01 |
#>   5.2 |        3.79e-02 |            0.00 |
#>   6.1 |            0.00 |        4.81e-03 |
#>   6.2 |        1.81e-03 |            0.00 |
#>   7.1 |            0.00 |        2.29e-04 |
#>   7.2 |        8.68e-05 |            0.00 |
#>   8.1 |            0.00 |        1.10e-05 |
#>   8.2 |        4.15e-06 |            0.00 |
#> 
#> ✔ Convergence achieved at iteration 8.
#> ────────────────────────────────────────────────────────────────────────────────
str(rf_ite, give.attr=FALSE)
#>  num [1:525, 1:28] 324066 98945 64400 80221 20703 ...
```

The cross-temporal method by Kourentzes & Athanasopoulos
([2019](#ref-Kourentzes2019-dj)), involves three steps: first,
reconciling forecasts for each time series at different temporal
aggregation levels using temporal hierarchies; second, performing
cross-sectional reconciliation at each temporal aggregation order; and
third, averaging the projection matrices from the second step and using
them to cross-sectionally reconcile the forecasts from the first step.
In contrast, we can reverses these steps starting with cross-sectional
reconciliation followed by temporal reconciliation ([Di Fonzo &
Girolimetto, 2023a](#ref-Di_Fonzo2023-dg)).

``` r

rf_tcs <- tcsrec(base = base_mat, res = res_mat,
                 cslist = list(agg_mat = vnaggmat, comb = "shr"), 
                 telist = list(agg_order = m, comb = "wlsv"))
str(rf_tcs, give.attr=FALSE)
#>  num [1:525, 1:28] 323554 98860 64315 80091 20660 ...
rf_cst <- cstrec(base = base_mat, res = res_mat,
                 cslist = list(agg_mat = vnaggmat, comb = "shr"), 
                 telist = list(agg_order = m, comb = "wlsv"))
str(rf_cst, give.attr=FALSE)
#>  num [1:525, 1:28] 324708 99123 64462 80416 20802 ...
```

Finally we can obtained the optimal (in least squares sense) combination
cross-temporal reconciled forecast ([Di Fonzo & Girolimetto,
2023a](#ref-Di_Fonzo2023-dg); [Girolimetto et al.,
2024](#ref-Girolimetto2023-jm)).

``` r

rf_opt <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                res = res_mat, comb = "wlsv", approach = "strc")
str(rf_opt, give.attr=FALSE)
#>  num [1:525, 1:28] 314298 97383 62630 77794 19695 ...
```

The following table shows some options for the optimal combination
cross-temporal reconciliation function
[`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md).

[TABLE]

### Practical challenges

#### Non negativity issues

It is not always true that reconciled forecasts will remain non-negative
even if the base forecasts are non-negative. For example, consider the
case of an identity covariance matrix (`"ols"`).

``` r

rf_ols <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                comb = "ols", approach = "strc")
recoinfo(rf_ols)
#> ✔ Optimal Cross-temporal Forecast Reconciliation
#> ℹ FoReco function: `ctrec`
#> ℹ Covariance approximation: `ols`
#> ℹ Non-negative forecasts: `FALSE`
```

To address this issue, we can use two approaches:

- State-of-the-art numerical optimization procedure, **osqp** ([Stellato
  et al., 2020](#ref-Stellato2020-jd)).

``` r

rf_osqp <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat, 
                 approach = "strc", comb = "ols", nn = "osqp")
tmp <- recoinfo(rf_osqp)
#> ✔ Optimal Cross-temporal Forecast Reconciliation
#> ℹ FoReco function: `ctrec`
#> ℹ Covariance approximation: `ols`
#> ℹ Non-negative forecasts: `TRUE`
tmp$info # OSQP information matrix 
#>         obj_val run_time iter     prim_res status status_polish
#> 1 -220421550178 11.75266  100 7.152467e-22      1             1
```

- Simple heuristic strategy: set-negative-to-zero, **sntz** ([Di Fonzo &
  Girolimetto, 2023b](#ref-Di_Fonzo2023-ae)).

``` r

rf_sntz <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                 approach = "strc", comb = "ols", nn = "sntz")
recoinfo(rf_sntz)
#> ✔ Optimal Cross-temporal Forecast Reconciliation
#> ℹ FoReco function: `ctrec`
#> ℹ Covariance approximation: `ols`
#> ℹ Non-negative forecasts: `TRUE`
```

#### A priori constrained (immutable) forecasts

Sometimes we may wish to incorporate a priori knowledge during the
reconciliation process ([Zhang et al., 2023](#ref-Zhang2023-yj)) in
order to improve the accuracy of the reconciled forecasts. For example,
suppose we want to fix the forecasts of the states level for the annual
series at the base forecasts values.

``` r

rf_imm <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat, approach = "strc",
                res = res_mat, comb = "wlsv", immutable = cbind(2:8, 12, 1))
str(rf_imm, give.attr=FALSE)
#>  num [1:525, 1:28] 321957 99476 63854 79348 20809 ...
```

``` r

round(rf_imm[2:8,1] - base_mat[2:8,1], 6)
#> A B C D E F G 
#> 0 0 0 0 0 0 0
```

#### Exploring a subset of temporal aggregation orders

Our approaches so far have involved considering all factors of \\m\\ as
potential aggregation orders. Nevertheless, it is worth noting that we
could also focus only on a given subset of these factors. For example,
we could be interested only on monthly, quarterly and annual forecasts.

``` r

te_subset <- c(12, 3, 1)
base_mat2 <- t(Reduce(rbind, base[paste0("k-", te_subset)]))
res_mat2 <- t(Reduce(rbind, res[paste0("k-", te_subset)]))
rf_sub <- ctrec(base = base_mat2, agg_order = te_subset, agg_mat = vnaggmat,
                res = res_mat2, comb = "wlsv", approach = "strc")
str(rf_sub, give.attr=FALSE)
#>  num [1:525, 1:17] 313401 97327 62517 77557 19511 ...
```

### Probabilistic forecast reconciliation

Girolimetto et al. ([2024](#ref-Girolimetto2023-jm)) extends to the
cross-temporal framework the probabilistic results presented by
Panagiotelis et al. ([2023](#ref-Panagiotelis2023-se)) for the
cross-sectional reconciliation. The reconciliation of probabilistic
forecasts is a two-step process: first, we sample from an incoherent
distribution, and then we reconcile each sample.

We can use a non-parametric method, the joint block bootstrap to
simulate \\B\\ samples and then reconciled them.

``` r

# Base forecasts' sample
B <- 100
# Base forecasts' sample
base_ctjb <- ctboot(model, B, agg_order = m)$sample 
str(base_ctjb[1:3], give.attr=FALSE)
#> List of 3
#>  $ : num [1:525, 1:28] 322620 94800 62004 83040 21827 ...
#>  $ : num [1:525, 1:28] 335643 97756 64587 83180 21143 ...
#>  $ : num [1:525, 1:28] 356882 102430 72814 86461 21771 ...

reco_ctjb <- ctsmp(simplify2array(base_ctjb), agg_order = m, agg_mat = vnaggmat,
                      res = res_mat, comb = "wlsv", approach = "strc", nn = "sntz")
class(reco_ctjb) # distribution object
#> [1] "distribution" "vctrs_vctr"   "list"

# Extracts mean:
mean_ctjb <- mean(reco_ctjb)
mean_ctjb <- as_ctmatrix(mean_ctjb, agg_order = m, 
                         n = sum(dim(vnaggmat)), 
                         row_names = unlist(dimnames(vnaggmat)))
str(mean_ctjb, give.attr = FALSE)
#>  num [1:525, 1:28] 325593 98890 65031 80722 20362 ...
```

A parametric method assumes a normal distribution (Gaussian), to
generate the incoherent sample set of forecasts. Since we have to
simulate from a multivariate normal distribution with a size of 14700,
we will use a diagonal covariance matrix in this vignette. However, it’s
important to note that this choice will result in a significantly narrow
variance for the reconciled forecasts.

``` r

# Gaussian reconciled distribution
reco_ctg_wlsv <- ctmvn(base = base_mat, agg_mat = vnaggmat, 
                         agg_order = m, comb = "wlsv", res = res_mat)
#> Warning in asMethod(object): sparse->dense coercion: allocating vector of size
#> 1.6 GiB
reco_ctg_wlsv # distribution object
#> <distribution[1]>
#>      tao-1 
#> MVN[14700]

# Gaussian reconciled distribution with different base covariance matrix
# Multi-step residuals
hres <- lapply(model, function(fit)
  lapply(1:frequency(fit[[1]]$x), function(h) 
    sapply(fit, residuals, type='response', h = h)))
hres_ct <- t(Reduce("rbind", lapply(hres, arrange_hres)))
# Re-arrenge multi-step residuals in a matrix form
mres <- as_hstack_ctlayout(hres_ct, agg_order = m)
# cov_shr <- shrink_estim(na.omit(mres)) # Time and computational intensive to use, but the better one
cov_wls <- diag(x = diag(cov(na.omit(mres))))

# Computational intensive:
# reco_ctg_wls <- ctmvn(base = base_mat, agg_mat = vnaggmat, 
#                       agg_order = m, comb = "wlsv", res = res_mat, 
#                       comb_base = cov_wls)
# reco_ctg_wls # distribution object

# Reconciled sample distribution starting from Gaussian base forecasts - faster approach
base_cts <- GMCM::rmvnormal(B, mu = as_hstack_ctlayout(base_mat, agg_order = m), 
                            sigma = cov_wls)
base_cts <- apply(base_cts, 1, as_ctmatrix, agg_order = m, 
                  n = sum(dim(vnaggmat)),
                  row_names = unlist(dimnames(vnaggmat)),
                  simplify = FALSE)
base_cts <- simplify2array(base_cts)
reco_cts <- ctsmp(base_cts, agg_order = m, agg_mat = vnaggmat,
                     res = res_mat, comb = "wlsv", approach = "strc", nn = "sntz")
reco_cts # distribution object
#> <distribution[1]>
#>       tao-1 
#> sample[100]
```

## `itagdp`: general linearly constrained multiple quarterly time series

In this section, we work with the `itagdp` dataset and the corresponding
zero-constrained matrix `gdpconsmat`. This dataset illustrates
reconciliation under more complex linear constraints where an unique
aggregation is not available. See the [dataset
vignette](https://danigiro.github.io/FoReco/articles/Dataset-vndata-and-itagdp.md)
for more details.

``` r

data(itagdp)      # dataset
data(gdpconsmat)    # Agg mat matrix
```

### Base forecast

We fit ARIMA models to each series and generate base forecasts. The data
is a quarterly multivariate time series so we obtain four-, two-, and
one-step-ahead base forecasts from the quarterly data and the
aggregation over 2, and 4 quarters.

``` r

te_set <- tetools(4)$set
m <- max(te_set)
data_k <- aggts(itagdp, te_set)
model <- setNames(rep(list(setNames(vector(mode='list', length=NCOL(itagdp)), colnames(itagdp))), 
                      length(te_set)), paste0("k-", te_set))
fc_obj <- setNames(rep(list(setNames(vector(mode='list', length=NCOL(itagdp)), colnames(itagdp))), 
                      length(te_set)), paste0("k-", te_set))

for(k in te_set){
  idk <- paste0("k-", k)
  for(i in 1:NCOL(itagdp)){
    ids <- colnames(itagdp)[i]
    model[[idk]][[ids]] <- ets(data_k[[idk]][,i])
    fc_obj[[idk]][[ids]] <- forecast(model[[idk]][[ids]], h = m/k)
  }
  cat(k, " ")
}
#> 4  2  1
```

``` r

# Point forecasts
base <- lapply(fc_obj, function(x) rbind(sapply(x, function(y) y$mean)))
str(base, give.attr = FALSE)
#> List of 3
#>  $ k-4: num [1, 1:21] 1811793 736509 1748222 1421782 326061 ...
#>  $ k-2: num [1:2, 1:21] 884015 928287 351866 379716 855453 ...
#>  $ k-1: num [1:4, 1:21] 431993 451829 449546 480841 167984 ...

# Residuals
res <- lapply(fc_obj, function(x) rbind(sapply(x, residuals, type = "response")))
str(res, give.attr = FALSE)
#> List of 3
#>  $ k-4: num [1:20, 1:21] -45068 22792 8057 8051 22851 ...
#>  $ k-2: num [1:40, 1:21] -14138 421 9758 -4257 2225 ...
#>  $ k-1: num [1:80, 1:21] -2778 -1745 -448 4068 3088 ...
```

### Point forecast reconciliation

We apply the optimal reconciliation method to the base forecasts,
considering the linear constraints defined by `gdpconsmat`.

``` r

base_mat <- t(Reduce(rbind, base))
res_mat <- t(Reduce(rbind, res))
rf_opt <- ctrec(base = base_mat, agg_order = m, cons_mat = gdpconsmat,
                res = res_mat, comb = "bdshr")
str(rf_opt, give.attr = FALSE)
#>  num [1:21, 1:7] 1807755 730746 1739243 1417381 321862 ...
```

#### Practical challenge: immutable forecast

In this case, we want to fix the forecasts of the top level series
(\\GDP\\) for the annual temporally aggreated series (\\k = 4\\) at the
base forecasts values.

``` r

rf_imm <- ctrec(base = base_mat, agg_order = m, cons_mat = gdpconsmat,
                res = res_mat, comb = "wlsv", immutable = rbind(c(1,4,1)))
str(rf_imm, give.attr = FALSE)
#>  num [1:21, 1:7] 1811793 731220 1743712 1418685 325027 ...
```

``` r

rf_imm[1,1] - base_mat[1,1]
#> GDP 
#>   0
```

### Probabilistic forecast reconciliation

We can use a non-parametric method, the joint block bootstrap to
simulate \\B\\ samples and then reconciled them.

``` r

# Base forecasts' sample
B <- 100
# Base forecasts' sample
base_ctjb <- ctboot(model, B, agg_order = m)$sample 
str(base_ctjb[1:3], give.attr=FALSE)
#> List of 3
#>  $ : num [1:21, 1:7] 1821647 743804 1770804 1442052 333871 ...
#>  $ : num [1:21, 1:7] 1822669 743494 1766008 1422906 349383 ...
#>  $ : num [1:21, 1:7] 1715859 708966 1655955 1382431 281299 ...

reco_ctjb <- ctsmp(simplify2array(base_ctjb), agg_order = m, cons_mat = gdpconsmat,
                      res = res_mat, comb = "wlsv")
reco_ctjb # distribution object
#> <distribution[1]>
#>       tao-1 
#> sample[100]

# Extracts mean:
mean_ctjb <- mean(reco_ctjb)
mean_ctjb <- as_ctmatrix(mean_ctjb, agg_order = m, 
                         n = NCOL(gdpconsmat), 
                         row_names = colnames(gdpconsmat))
str(mean_ctjb, give.attr = FALSE)
#>  num [1:21, 1:7] 1817045 734150 1748060 1419745 328316 ...
```

Alternatively, we can use a parametric method.

``` r

# Gaussian reconciled distribution
reco_ctg_wlsv <- ctmvn(base = base_mat, cons_mat = gdpconsmat,
                         agg_order = m, comb = "wlsv", res = res_mat)
reco_ctg_wlsv # distribution object
#> <distribution[1]>
#>    tao-1 
#> MVN[147]

# Gaussian reconciled distribution with different base covariance matrix
# Multi-step residuals
hres <- lapply(model, function(fit)
  lapply(1:frequency(fit[[1]]$x), function(h) 
    sapply(fit, residuals, type='response', h = h)))
hres_ct <- t(Reduce("rbind", lapply(hres, arrange_hres)))
# Re-arrenge multi-step residuals in a matrix form
mres <- as_hstack_ctlayout(hres_ct, agg_order = m)
# cov_shr <- shrink_estim(na.omit(mres)) # Time and computational intensive to use, but the better one
cov_wls <- diag(x = diag(cov(na.omit(mres))))

reco_ctg_wls <- ctmvn(base = base_mat, cons_mat = gdpconsmat,
                        agg_order = m, comb = "wlsv", res = res_mat, 
                        comb_base = cov_wls)
reco_ctg_wls # distribution object
#> <distribution[1]>
#>    tao-1 
#> MVN[147]

# Reconciled sample distribution starting from Gaussian base forecasts - faster approach
base_cts <- GMCM::rmvnormal(B, mu = as_hstack_ctlayout(base_mat, agg_order = m), 
                            sigma = cov_wls)
base_cts <- apply(base_cts, 1, as_ctmatrix, agg_order = m, 
                  n = NCOL(gdpconsmat), 
                  row_names = colnames(gdpconsmat),
                  simplify = FALSE)
base_cts <- simplify2array(base_cts)
reco_cts <- ctsmp(base_cts, agg_order = m, cons_mat = gdpconsmat,
                     res = res_mat, comb = "wlsv")
reco_cts # distribution object
#> <distribution[1]>
#>       tao-1 
#> sample[100]
```

## References

Di Fonzo, T., & Girolimetto, D. (2023a). Cross-temporal forecast
reconciliation: Optimal combination method and heuristic alternatives.
*International Journal of Forecasting*, *39*(1), 39–57.
<https://doi.org/10.1016/j.ijforecast.2021.08.004>

Di Fonzo, T., & Girolimetto, D. (2023b). Spatio-temporal reconciliation
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

Kourentzes, N., & Athanasopoulos, G. (2019). Cross-temporal coherent
forecasts for australian tourism. *Annals Of Tourism Research*, *75*,
393–409. <https://doi.org/10.1016/j.annals.2019.02.001>

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

Wickramasuriya, S. L., Turlach, B. A., & Hyndman, R. J. (2020). Optimal
non-negative forecast reconciliation. *Statistics and Computing*,
*30*(5), 1167–1182. <https://doi.org/10.1007/s11222-020-09930-0>

Zhang, B., Kang, Y., Panagiotelis, A., & Li, F. (2023). Optimal
reconciliation with immutable forecasts. *European Journal of
Operational Research*, *308*(2), 650–660.
<https://doi.org/10.1016/j.ejor.2022.11.035>
