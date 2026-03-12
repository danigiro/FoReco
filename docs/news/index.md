# Changelog

## FoReco 1.2.0

- (CRAN Package Check Results) Fixed compatibility issues with
  [`osqp`](https://CRAN.R-project.org/package=osqp) 1.0. The package is
  now fully compatible with any version.
- Added
  [`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md),
  [`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md), and
  [`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md) for
  Gaussian probabilistic forecast reconciliation in the cross-sectional,
  temporal, and cross-temporal frameworks using the
  [`distributional`](https://CRAN.R-project.org/package=distributional)
  package;
- Added
  [`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md),
  [`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md), and
  [`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md) for
  sample-based probabilistic forecast reconciliation in the
  cross-sectional, temporal, and cross-temporal frameworks using the
  [`distributional`](https://CRAN.R-project.org/package=distributional)
  package;
- Added
  [`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md)
  and `as_horizon_stacked_ctmatrix()` functions to convert between
  horizon-stacked (cross-temporal version) and cross-temporal layouts;
- Added
  [`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md)
  and `as_horizon_stacked_tematrix()` functions to convert between
  horizon-stacked (temporal version) and temporal layouts;
- Added non-negative forecast reconciliation algorithms *bpv* (block
  principal pivoting algorithm), *nfca* (negative forecasts correction
  algorithm), *nnic* (iterative non-negative reconciliation with
  immutable constraints), and *sntz* (set-negative-to-zero with
  bottom-up and top-down alternatives) based on:
  - Girolimetto, D. (2025), Non-negative forecast reconciliation:
    Optimal methods and operational solutions. *arXiv*;
  - Kourentzes, N. and Athanasopoulos, G. (2021) Elucidate structure in
    intermittent demand series. *European Journal of Operational
    Research*, 288, 141-152.
    [doi:10.1016/j.ejor.2020.05.046](https://doi.org/10.1016/j.ejor.2020.05.046);
  - Wickramasuriya, S. L., Turlach, B. A., and Hyndman, R. J. (2020),
    “Optimal non-negative forecast reconciliation”, *Statistics and
    Computing*, 30(5), 1167–1182.
    [doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0);
- Added `...` for [`simulate()`](https://rdrr.io/r/stats/simulate.html)
  additional arguments in
  [`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md),
  [`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md)
  and
  [`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md);
- Fixed bugs and improved stability.

## FoReco 1.1.0

CRAN release: 2025-06-07

#### New Features

- Added new `bpv` non-negative forecast reconciliation algorithm
  (experimental) based on the cross-sectional framework presented in
  Wickramasuriya et al. (2020) and, now, extended for the temporal and
  cross-temporal framework:
  - Wickramasuriya, S. L., Turlach, B. A., and Hyndman, R. J. (2020),
    “Optimal non-negative forecast reconciliation”, *Statistics and
    Computing*, 30(5), 1167–1182.
    [doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0);
- New `oasd` cross-sectional covariance matrix (experimental),
  implementing an oracle shrunk covariance estimation (Ando and Xiao,
  2023):
  - Ando, S., and Xiao, M. (2023), “High-dimensional covariance matrix
    estimation: shrinkage toward a diagonal target”, *IMF Working
    Papers*, 2023(257), A001;
- Redesigned `bounds` parameter to enable bounded forecast
  reconciliation for
  [`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
  [`terec()`](https://danigiro.github.io/FoReco/reference/terec.md), and
  [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md)
  functions;
- Introduced new
  [`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md)
  function to define custom bounds for reconciliation.

#### Bug Fixes

- Fixed bug when only a subset of `agg_order` factors was selected in
  [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md),
  [`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md),
  and
  [`cstrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)
  functions.

## FoReco 1.0.0

CRAN release: 2024-08-20

> **Note** – The latest release of FoReco introduces significant changes
> to its function notation and adds several new features. This major
> update, FoReco 1.0, is not compatible with previous versions due to
> the substantial changes made to the package’s core structure. The
> previous version is available on
> [Github](https://github.com/danigiro/FoReco026)
> ([docs](https://danigiro.github.io/FoReco026/)).
>
> Due to the significant changes in FoReco 1.0, users are advised to
> carefully review the updated documentation and examples before using
> the new version. The latest documentation and release notes are
> available on
> [danigiro.github.io/FoReco/](https://danigiro.github.io/FoReco/)

- **Updated Function Notation:** All functions related to
  cross-sectional, temporal, and cross-temporal frameworks now use the
  prefixes **cs**, **te**, and **ct**, respectively. For example, the
  optimal combination reconciliation functions are now
  [`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
  [`terec()`](https://danigiro.github.io/FoReco/reference/terec.md), and
  [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md).

- **Simplified Function Outputs:** Reconciliation functions now return
  only matrices. Additional information can be accessed using
  `attr(., "FoReco")` or the
  [`recoinfo()`](https://danigiro.github.io/FoReco/reference/recoinfo.md)
  function.

- **New Datasets:** Two new datasets, `itagdp` (Italian Quarterly
  National Accounts) and `vndata` (Australian Tourism Demand), are
  included along with their respective aggregation or constraint
  matrices.

- **Classic Approach:** The middle-out approach
  ([`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md),
  [`temo()`](https://danigiro.github.io/FoReco/reference/temo.md), and
  [`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md)) has
  been implemented alongside the classic bottom-up
  ([`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md),
  [`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md), and
  [`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md)) and
  top-down
  ([`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md),
  [`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md), and
  [`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md))
  methods.

- **Level Conditional Coherent Reconciliation:** Level conditional
  coherent reconciliation is now available for all constraints:
  [`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md)
  (cross-sectional),
  [`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md)
  (temporal), and
  [`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md)
  (cross-temporal).

- **Immutable reconciliation:** The `immutable()` parameter has been
  added to the reconciliation functions
  ([`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md),
  [`terec()`](https://danigiro.github.io/FoReco/reference/terec.md), and
  [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md)) to
  prevent the base forecasts from being modified with both the
  structural (`approach='strc'`) and the zero-constrained
  (`approach='proj'`) approach.

- **Balanced and unbalanced hierarchy:** added
  [`balance_hierarchy()`](https://danigiro.github.io/FoReco/reference/balance_hierarchy.md)
  and
  [`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)
  for dealing with balanced and unbalanced hierarchies.

- **Projection Matrix Functions:** Functions
  [`csprojmat()`](https://danigiro.github.io/FoReco/reference/csprojmat.md),
  [`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md),
  and
  [`ctprojmat()`](https://danigiro.github.io/FoReco/reference/ctprojmat.md)
  have been added to obtain projection matrices.

- **Covariance Matrix Functions:** Functions
  [`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md),
  [`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md), and
  [`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md) have
  been added to obtain covariance matrices.

- **Function Renaming:** Several functions have been renamed to improve
  consistency and clarity

  - `Cmatrix()` -\>
    [`df2aggmat()`](https://danigiro.github.io/FoReco/reference/df2aggmat.md)
  - `hts_tools()` -\>
    [`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)
  - `thf_tools()` -\>
    [`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)
  - `ctf_tools()` -\>
    [`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md)
  - `agg_ts()` -\>
    [`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md)
  - `residuals_matrix()` -\>
    [`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md)
  - `boot_cs()` -\>
    [`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md)
  - `boot_te()` -\>
    [`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md)
  - `boot_ct()` -\>
    [`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md)
  - `htsrec()` -\>
    [`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md)
  - `thfrec()` -\>
    [`terec()`](https://danigiro.github.io/FoReco/reference/terec.md)
  - `octrec()` -\>
    [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md)
  - `lccrec()` -\>
    [`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md)

## FoReco 0.2.6

CRAN release: 2023-05-16

###### Major changes: probabilistic forecast reconciliation

- Added `boot_cs()`, `boot_te()` and `boot_ct()` to draw samples from,
  respectively, cross-sectional, temporal and cross-temporal joint
  (block) bootstrap.

###### Minor changes

- Fixed deprecation warnings in
  [Matrix](https://CRAN.R-project.org/package=Matrix) (v. 1.5-0);
- Improved docs and bug fixes;
- Fixed [`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md)
  inputs;
- Added
  [`FoReco2matrix()`](https://danigiro.github.io/FoReco/reference/FoReco2matrix.md)
  to transform FoReco forecasts input and output in a list of
  matrix/vector class;
- Added `agg_ts()`: non-overlapping temporal aggregation of a time
  series according to a specific aggregation order;
- Added
  [`arrange_hres()`](https://danigiro.github.io/FoReco/reference/residuals.md)
  and `residuals_matrix()` functions to arrange residuals for the
  covariance matrix under Gaussianity.

## FoReco 0.2.5

CRAN release: 2022-07-04

###### Minor changes

- Fixed the not negative reconciliation “**sntz**” for `octrec()`;
- Fixed documentation.

## FoReco 0.2.4

CRAN release: 2022-06-16

###### Major changes

- Added
  [`lcmat()`](https://danigiro.github.io/FoReco/reference/lcmat.md)
  function.

###### Minor changes

- Fixed BU approach when the number of columns of basef is equal to the
  number of bottom time series `htsrec()`;
- Fixed `score_index()`;
- Fixed the `bounds` param when `type = "S"` in `htsrec()`, `thfrec()`
  and `octrec()`;
- Add the possibility to fix base forecasts through the `v` param in
  `htsrec()`, `thfrec()` and `octrec()` - experimental;
- Add two new type of optimal cross-temporal reconciliation
  (**cs_struc** and **t_struc**);
- Improved docs and bug fixes.

## FoReco 0.2.2

CRAN release: 2022-02-17

###### Minor changes

- Fixed documentation;
- Removed `ut2c()` and `srref()`.

## FoReco 0.2.1

CRAN release: 2021-07-23

###### Minor changes

- Fixed a bug in the output of `lccrec()` (now the function returns the
  Level Conditional Coherent or the Combined Conditional Coherent
  forecasts);
- Fixed the not negative reconciliation “**KAnn**” when `keep = "list"`.

## FoReco 0.2.0

CRAN release: 2021-05-21

###### Major changes

- It’s possible to use a subset of factors of *m* (max. order of
  temporal aggregation);
- Added the possibility for `htsrec()`, `thfrec()` and `octrec()` to
  introduce a list of *h* covariance matrices in the parameters `W` and
  `Omega`, where *h* stands for the forecast horizon (note that for
  `thfrec()` and `octrec()` this is the forecast horizon of the entire
  cycle);
- Param `Sstruc` is no more avaible in `octrec()` and `ctf_tools()`.
  FoReco uses a fast algorithm to compute **Scheck**, so no external
  input is needed;
- Modified output of `ctf_tools()` (added `Ccheck`, `Htcheck`, `Scheck`,
  removed `Cstruc`, `Sstruc`), `hts_tools()` (added `C`) and
  `thf_tools()` (added `m`);
- Added two new not negative reconciliation techniques (“**KAnn**” and
  “**sntz**”) with a new parameter (`nn_type`) in `htsrec()`, `thfrec()`
  and `octrec()`;
- Added the top-down reconciliation function `tdrec()`;
- Added the level conditional forecast reconciliation (with and without
  not-negative constraints) for genuine hierarchical/grouped time series
  `levrec()` (cross-sectional, temporal and cross-temporal).

###### Minor changes

- Now in `octrec()` it is also possible to introduce the **Ω**
  covariance matrix variant through the `Omega` parameter and not only
  the **W** variant with the `W` parameter;
- Updated
  [`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md),
  [`cstrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)
  and
  [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md).
  In the
  [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md)
  function the `maxit` parameter has been replaced by `itmax`, however
  for the moment `maxit` is still supported;
- Now FoReco removes null rows from the cross-sectional aggregation
  matrix **C** and it warns the user if the balanced version of an
  unbalanced hierarchy is considering duplicated variables;
- Redesigned the console output and added a new convergence norm as
  *default* for
  [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md)
  (`norm` parameter).

###### Experimental

- Add the possibility to introduce constraints through the `bounds`
  param in `htsrec()`, `thfrec()` and `octrec()`;
- Add a function `oct_bounds()` to organize the bounds on a specific
  dimension (i.e. only cross-sectional or only temporal) in a
  cross-temporal framework;
- Added `ut2c()` and `srref()` to develop a cross-sectional structural
  representation starting from a zero constraints kernel matrix;
- Added in `score_index()` the calculation of multiple forecast horizons
  index (like 1:6) and multiple cross-sectional levels for a forecasting
  experiment.

## FoReco 0.1.1

CRAN release: 2020-10-17

Minore release, fixing some bugs and the documentation

- Fixed a bug in
  [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md)
  when calculating the incoherence
- Fixed documentation
- Changed the contact mail (now it’s <daniele.girolimetto@phd.unipd.it>)
- Corrected the second section of the vignette
  “`Average relative accuracy indices`”

## FoReco 0.1.0

CRAN release: 2020-10-01

- Release on github
