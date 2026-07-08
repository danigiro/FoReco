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
summary(rf_bu)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctbu`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2  k-4 h-3
#> Total 280205.66 147023.22 133182.45 108461.14 82962.06 88782.47
#> A      89744.37  48054.93  41689.44  36172.94 25354.01 28217.43
#> B      55472.04  31371.67  24100.37  24039.84 14914.38 16517.81
#> C      69366.03  32688.48  36677.55  22773.06 22923.39 23669.58
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
```

To obtain a list of forecasts at different orders of aggregation, we can
use the `components` function.

``` r

str(components(rf_bu), give.attr=FALSE)
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
summary(rf_td_gsa)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cttd`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 327178.11 169754.43 157423.68 124816.10 97513.30 104848.71
#> A     105756.54  56443.82  49312.73  42417.28 29878.99  33460.27
#> B      63172.87  35416.52  27756.35  26961.07 17218.82  18992.98
#> C      85397.48  39956.82  45440.66  27996.67 28007.38  29393.44
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.

# Proportions of the historical averages - Gross-Sohl method F
p_gsf <- apply(bts_mat, 2, function(y){
  colMeans(matrix(y, ncol = m, byrow = TRUE))/mean(tot_12)
})
rf_td_gsf <- cttd(fc_tot_12, agg_order = m, weights = t(p_gsf), agg_mat = vnaggmat)
summary(rf_td_gsf)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cttd`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 327178.11 169749.01 157429.10 124783.34 97511.04 104883.72
#> A     105655.72  56376.47  49279.25  42338.06 29876.08  33441.59
#> B      63188.00  35423.12  27764.88  26968.92 17211.76  19007.32
#> C      85289.86  39911.00  45378.86  27963.87 27969.69  29356.30
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_lcc)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctlcc`
#> • Covariance approximation approach: `wlsv`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 316724.44 164732.15 151992.29 120461.49 95042.00 101220.94
#> A      98441.32  52330.36  46110.96  39107.44 28131.83  31202.05
#> B      62490.88  34962.18  27528.70  26536.98 17133.05  18820.85
#> C      79421.25  37508.95  41912.30  26046.70 26248.17  27126.39
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_ite)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `iterec`
#> • Covariance approximation approach: `te-wlsv` and `cs-shr`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 324065.94 168163.59 155902.34 122974.31 97266.86 103824.77
#> A      98945.29  52542.01  46403.28  39243.88 28248.40  31453.01
#> B      64400.49  35944.50  28456.00  27276.52 17656.91  19467.06
#> C      80221.08  37794.23  42426.85  26292.03 26508.79  27420.25
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_tcs)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `tcsrec`
#> • Covariance approximation approach: `te-wlsv` and `cs-shr`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 323553.55 167866.29 155687.26 122761.77 97155.82 103635.96
#> A      98860.45  52503.98  46356.46  39233.59 28212.93  31413.94
#> B      64314.67  35894.65  28420.03  27239.58 17635.83  19439.26
#> C      80090.80  37731.19  42359.60  26253.76 26455.23  27381.81
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
rf_cst <- cstrec(base = base_mat, res = res_mat,
                 cslist = list(agg_mat = vnaggmat, comb = "shr"), 
                 telist = list(agg_order = m, comb = "wlsv"))
summary(rf_cst)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `cstrec`
#> • Covariance approximation approach: `cs-shr` and `te-wlsv`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 324708.28 168486.31 156221.98 123231.54 97397.10 104079.64
#> A      99123.23  52608.99  46514.24  39280.79 28312.19  31530.25
#> B      64462.48  35965.55  28496.94  27296.45 17662.60  19503.43
#> C      80416.33  37890.65  42525.68  26369.61 26560.39  27486.33
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
```

Finally we can obtained the optimal (in least squares sense) combination
cross-temporal reconciled forecast ([Di Fonzo & Girolimetto,
2023a](#ref-Di_Fonzo2023-dg); [Girolimetto et al.,
2024](#ref-Girolimetto2023-jm)).

``` r

rf_opt <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                res = res_mat, comb = "wlsv", approach = "strc")
summary(rf_opt)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `wlsv`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 314297.92 163231.20 151066.72 119528.49 94308.64 100460.79
#> A      97383.15  51787.97  45595.18  38774.36 27716.53  30892.26
#> B      62630.42  35071.27  27559.16  26676.44 17104.91  18849.07
#> C      77793.96  36597.82  41196.14  25458.94 25729.03  26606.00
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_ols)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `ols`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 327273.55 169942.50 157331.05 124860.62 97691.28 104721.64
#> A     100541.58  53506.35  47035.22  40159.79 28415.26  31966.53
#> B      64915.19  36333.74  28581.45  27707.87 17688.32  19519.00
#> C      79771.73  37507.25  42264.49  26215.22 26231.87  27324.65
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
```

To address this issue, we can use two approaches:

- State-of-the-art numerical optimization procedure, **osqp** ([Stellato
  et al., 2020](#ref-Stellato2020-jd)).

``` r

rf_osqp <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat, 
                 approach = "strc", comb = "ols", nn = "osqp")
summary(rf_osqp) # OSQP information matrix 
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `ols`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>         obj_val run_time iter     prim_res status status_polish
#> 1 -220421550178  12.5692  100 7.152467e-22      1             1
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 327274.34 169943.31 157331.02 124861.19 97691.27 104721.87
#> A     100540.92  53505.68  47035.24  40159.31 28415.27  31966.34
#> B      64914.55  36333.08  28581.47  27707.40 17688.33  19518.82
#> C      79772.05  37507.45  42264.60  26215.47 26231.76  27324.83
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
```

- Simple heuristic strategy: set-negative-to-zero, **sntz** ([Di Fonzo &
  Girolimetto, 2023b](#ref-Di_Fonzo2023-ae)).

``` r

rf_sntz <- ctrec(base = base_mat, agg_order = m, agg_mat = vnaggmat,
                 approach = "strc", comb = "ols", nn = "sntz")
summary(rf_sntz)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `ols`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 1.096725e-05 1.490116e-08    0
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 327467.74 170060.55 157407.19 124949.42 97724.57 104793.75
#> A     100541.58  53506.35  47035.22  40159.79 28415.26  31966.53
#> B      64915.19  36333.74  28581.45  27707.87 17688.32  19519.00
#> C      79855.06  37564.75  42290.31  26261.46 26243.13  27350.47
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_imm)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `wlsv`
#> • Output: (525 x 28) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1   k-6 h-1   k-6 h-2   k-4 h-1  k-4 h-2   k-4 h-3
#> Total 321957.07 167060.78 154896.29 122081.54 96861.69 103013.84
#> A      99475.83  52834.31  46641.52  39471.92 28414.09  31589.83
#> B      63853.54  35682.83  28170.72  27084.15 17512.62  19256.78
#> C      79347.67  37374.67  41973.00  25976.84 26246.93  27123.90
#> ... (521 more rows, 22 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
summary(rf_sub)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `wlsv`
#> • Output: (525 x 17) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 3, and 1
#> • Forecast horizons (h) per k: 1, 4, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Reconciled forecasts
#>        k-12 h-1  k-3 h-1  k-3 h-2  k-3 h-3  k-3 h-4   k-1 h-1
#> Total 313401.11 90989.24 71707.27 75420.23 75284.37 47568.514
#> A      97326.81 29573.17 22191.15 21835.63 23726.84 15990.347
#> B      62516.98 20929.41 14110.75 12826.44 14650.38 11532.116
#> C      77556.52 19098.27 17339.99 22163.56 18954.70  9779.126
#> ... (521 more rows, 11 more columns)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
base_ctjb <- ctboot(model, B, agg_order = m) 
str(base_ctjb[1:3], give.attr=FALSE)
#> List of 3
#>  $ : num [1:525, 1:28] 305455 93053 59753 76427 18089 ...
#>  $ : num [1:525, 1:28] 331164 98131 60354 87924 20222 ...
#>  $ : num [1:525, 1:28] 340214 101263 67203 82541 20915 ...

reco_ctjb <- ctsmp(simplify2array(base_ctjb), agg_order = m, agg_mat = vnaggmat,
                      res = res_mat, comb = "wlsv", approach = "strc", nn = "sntz")
summary(reco_ctjb)
#> ✔ Cross-temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctsmp(ctrec)`
#> • Covariance approximation approach: `wlsv`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 0.0007648468 1.490116e-08    0
#> 2 0.0007648468 1.490116e-08    0
#> 3 0.0007648468 1.490116e-08    0
#> 4 0.0007648468 1.490116e-08    0
#> 5 0.0007648468 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>       tao-1 
#> sample[100]

# Extracts mean:
mean_ctjb <- mean(reco_ctjb)
mean_ctjb <- as_ctmatrix(mean_ctjb, agg_order = m, 
                         n = sum(dim(vnaggmat)), 
                         row_names = unlist(dimnames(vnaggmat)))
str(mean_ctjb, give.attr = FALSE)
#>  num [1:525, 1:28] 323733 98517 64665 80830 20098 ...
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
#> All elements are shown. Use `print(x, n_row)` to limit the output.

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
summary(reco_cts)
#> ✔ Cross-temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctsmp(ctrec)`
#> • Covariance approximation approach: `wlsv`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 525
#> • Temporal orders (k): 12, 6, 4, 3, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, 3, 4, 6, and 12
#> • Non-negative forecasts (check): `TRUE`
#> 
#> ── Non-negative reconciliation diagnostics
#>       run_time          tol iter
#> 1 0.0008330345 1.490116e-08    0
#> 2 0.0008330345 1.490116e-08    0
#> 3 0.0008330345 1.490116e-08    0
#> 4 0.0008330345 1.490116e-08    0
#> 5 0.0008330345 1.490116e-08    0
#> ℹ Showing the first 5 rows of the non-negativity diagnostics info matrix.
#> 
#> ── Reconciled forecasts
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
summary(rf_opt)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `bdshr`
#> • Output: (21 x 7) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Temporal orders (k): 4, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, and 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>         k-4 h-1  k-2 h-1  k-2 h-2  k-1 h-1  k-1 h-2  k-1 h-3
#> GDP   1807755.5 882541.6 925213.9 431619.2 450922.3 447344.1
#> D1     730746.5 351702.2 379044.3 168133.2 183569.0 170574.4
#> P3_P5 1739243.0 851902.5 887340.5 420175.3 431727.2 426642.1
#> P3    1417380.6 696115.0 721265.6 345387.2 350727.8 353386.0
#> ... (17 more rows, 1 more column)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
```

#### Practical challenge: immutable forecast

In this case, we want to fix the forecasts of the top level series
(\\GDP\\) for the annual temporally aggreated series (\\k = 4\\) at the
base forecasts values.

``` r

rf_imm <- ctrec(base = base_mat, agg_order = m, cons_mat = gdpconsmat,
                res = res_mat, comb = "wlsv", immutable = rbind(c(1,4,1)))
summary(rf_imm)
#> ✔ Cross-temporal point forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctrec`
#> • Covariance approximation approach: `wlsv`
#> • Output: (21 x 7) matrix
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Temporal orders (k): 4, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, and 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#>         k-4 h-1  k-2 h-1  k-2 h-2  k-1 h-1  k-1 h-2  k-1 h-3
#> GDP   1811793.0 884122.9 927670.2 432344.7 451778.2 448468.2
#> D1     731219.8 351870.6 379349.2 168155.7 183715.0 170652.7
#> P3_P5 1743711.7 853543.5 890168.2 420927.2 432616.3 427971.6
#> P3    1418685.0 696685.6 721999.4 345652.1 351033.4 353726.5
#> ... (17 more rows, 1 more column)
#> Use `print(x, n_row, n_col)` to see more rows and columns.
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
base_ctjb <- ctboot(model, B, agg_order = m) 
str(base_ctjb[1:3], give.attr=FALSE)
#> List of 3
#>  $ : num [1:21, 1:7] 1825237 738841 1757021 1426951 331837 ...
#>  $ : num [1:21, 1:7] 1805001 737769 1751547 1431916 323729 ...
#>  $ : num [1:21, 1:7] 1803527 742817 1734769 1406326 330324 ...

reco_ctjb <- ctsmp(simplify2array(base_ctjb), agg_order = m, cons_mat = gdpconsmat,
                      res = res_mat, comb = "wlsv")
summary(reco_ctjb)
#> ✔ Cross-temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctsmp(ctrec)`
#> • Covariance approximation approach: `wlsv`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Temporal orders (k): 4, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, and 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
#> <distribution[1]>
#>       tao-1 
#> sample[100]

# Extracts mean:
mean_ctjb <- mean(reco_ctjb)
mean_ctjb <- as_ctmatrix(mean_ctjb, agg_order = m, 
                         n = NCOL(gdpconsmat), 
                         row_names = colnames(gdpconsmat))
str(mean_ctjb, give.attr = FALSE)
#>  num [1:21, 1:7] 1813092 732263 1742744 1417217 325527 ...
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
#> All elements are shown. Use `print(x, n_row)` to limit the output.

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
summary(reco_ctg_wls)
#> ✔ Cross-temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctmvn`
#> • Covariance approximation approach: `wlsv`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Temporal orders (k): 4, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, and 4
#> 
#> ── Reconciled forecasts
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
summary(reco_cts)
#> ✔ Cross-temporal probabilistic forecast reconciliation
#> 
#> ── Method
#> • Function used: `ctsmp(ctrec)`
#> • Covariance approximation approach: `wlsv`
#> • Output: distributional object
#> 
#> ── Structure
#> • Number of cross-sectional series: 21
#> • Temporal orders (k): 4, 2, and 1
#> • Forecast horizons (h) per k: 1, 2, and 4
#> • Non-negative forecasts (check): `FALSE`
#> 
#> ── Reconciled forecasts
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
