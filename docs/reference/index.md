# Package index

## Cross-sectional reconciliation

All the reconciliation functions for the cross-sectional framework

- [`csbu()`](https://danigiro.github.io/FoReco/reference/csbu.md) :
  Cross-sectional Bottom-up Reconciliation
- [`csmo()`](https://danigiro.github.io/FoReco/reference/csmo.md) :
  Cross-sectional Middle-out Reconciliation
- [`cstd()`](https://danigiro.github.io/FoReco/reference/cstd.md) :
  Cross-sectional Top-down Reconciliation
- [`cslcc()`](https://danigiro.github.io/FoReco/reference/cslcc.md) :
  Level Conditional Coherent Reconciliation for Genuine
  Hierarchical/Grouped Time Series
- [`csrec()`](https://danigiro.github.io/FoReco/reference/csrec.md) :
  Optimal Combination Cross-sectional Reconciliation
- [`cssmp()`](https://danigiro.github.io/FoReco/reference/cssmp.md) :
  Cross-sectional Probabilistic Reconciliation (Sample Approach)
- [`csmvn()`](https://danigiro.github.io/FoReco/reference/csmvn.md) :
  Cross-sectional Gaussian Probabilistic Reconciliation

## Temporal reconciliation

All the reconciliation functions for the temporal framework

- [`tebu()`](https://danigiro.github.io/FoReco/reference/tebu.md) :
  Temporal Bottom-up Reconciliation
- [`temo()`](https://danigiro.github.io/FoReco/reference/temo.md) :
  Temporal Middle-out Reconciliation
- [`tetd()`](https://danigiro.github.io/FoReco/reference/tetd.md) :
  Temporal Top-down Reconciliation
- [`telcc()`](https://danigiro.github.io/FoReco/reference/telcc.md) :
  Level Conditional Coherent Reconciliation for Temporal Hierarchies
- [`terec()`](https://danigiro.github.io/FoReco/reference/terec.md) :
  Optimal Combination Temporal Reconciliation
- [`tesmp()`](https://danigiro.github.io/FoReco/reference/tesmp.md) :
  Temporal Probabilistic Reconciliation (Sample Approach)
- [`temvn()`](https://danigiro.github.io/FoReco/reference/temvn.md) :
  Temporal Gaussian Probabilistic Reconciliation

## Cross-temporal reconciliation

All the reconciliation functions for the cross-temporal framework

- [`ctbu()`](https://danigiro.github.io/FoReco/reference/ctbu.md) :
  Cross-temporal Bottom-up Reconciliation
- [`ctmo()`](https://danigiro.github.io/FoReco/reference/ctmo.md) :
  Cross-temporal Middle-out Reconciliation
- [`cttd()`](https://danigiro.github.io/FoReco/reference/cttd.md) :
  Cross-temporal Top-down Reconciliation
- [`ctlcc()`](https://danigiro.github.io/FoReco/reference/ctlcc.md) :
  Level Conditional Coherent Reconciliation for Cross-temporal
  Hierarchies
- [`ctrec()`](https://danigiro.github.io/FoReco/reference/ctrec.md) :
  Optimal Combination Cross-temporal Reconciliation
- [`tcsrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)
  [`cstrec()`](https://danigiro.github.io/FoReco/reference/heuristic-reco.md)
  : Heuristic Cross-temporal Reconciliation
- [`iterec()`](https://danigiro.github.io/FoReco/reference/iterec.md) :
  Iterative Cross-temporal Reconciliation
- [`ctsmp()`](https://danigiro.github.io/FoReco/reference/ctsmp.md) :
  Cross-temporal Probabilistic Reconciliation (Sample Approach)
- [`ctmvn()`](https://danigiro.github.io/FoReco/reference/ctmvn.md) :
  Cross-temporal Gaussian Probabilistic Reconciliation

## Datasets

- [`itagdp`](https://danigiro.github.io/FoReco/reference/itagdp.md)
  [`outside`](https://danigiro.github.io/FoReco/reference/itagdp.md)
  [`expside`](https://danigiro.github.io/FoReco/reference/itagdp.md)
  [`incside`](https://danigiro.github.io/FoReco/reference/itagdp.md)
  [`gdpconsmat`](https://danigiro.github.io/FoReco/reference/itagdp.md)
  : Italian Quarterly National Accounts Dataset
- [`vndata`](https://danigiro.github.io/FoReco/reference/vndata.md)
  [`vnaggmat`](https://danigiro.github.io/FoReco/reference/vndata.md) :
  Australian Tourism Demand Dataset

## Reconciliation tools

Functions to compute the reconciliation matrices

- [`cstools()`](https://danigiro.github.io/FoReco/reference/cstools.md)
  : Cross-sectional Reconciliation Tools
- [`csprojmat()`](https://danigiro.github.io/FoReco/reference/csprojmat.md)
  : Projection Matrix for Optimal Combination Cross-sectional
  Reconciliation
- [`tetools()`](https://danigiro.github.io/FoReco/reference/tetools.md)
  : Temporal Reconciliation Tools
- [`teprojmat()`](https://danigiro.github.io/FoReco/reference/teprojmat.md)
  : Projection Matrix for Optimal Combination Temporal Reconciliation
- [`cttools()`](https://danigiro.github.io/FoReco/reference/cttools.md)
  : Cross-temporal Reconciliation Tools
- [`ctprojmat()`](https://danigiro.github.io/FoReco/reference/ctprojmat.md)
  : Projection Matrix for Optimal Combination Cross-temporal
  Reconciliation

## Covariance estimation

Functions to compute the covariance matrix estimators

- [`cscov()`](https://danigiro.github.io/FoReco/reference/cscov.md) :
  Cross-sectional Covariance Matrix Approximation
- [`tecov()`](https://danigiro.github.io/FoReco/reference/tecov.md) :
  Temporal Covariance Matrix Approximation
- [`ctcov()`](https://danigiro.github.io/FoReco/reference/ctcov.md) :
  Cross-temporal Covariance Matrix Approximation
- [`shrink_estim()`](https://danigiro.github.io/FoReco/reference/shrink_estim.md)
  : Shrinkage of the Covariance Matrix
- [`shrink_oasd()`](https://danigiro.github.io/FoReco/reference/shrink_oasd.md)
  : Shrinkage of the Covariance Matrix Using the Oracle Approximation

## Utilities

Functions to manipulate the data and the forecasts

- [`aggts()`](https://danigiro.github.io/FoReco/reference/aggts.md) :
  Non-Overlapping Temporal Aggregation of a Time Series
- [`as_ctmatrix()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md)
  [`as_hstack_ctlayout()`](https://danigiro.github.io/FoReco/reference/ctmatrix_layouts.md)
  : Convert Between Horizon-stacked and Cross-temporal Layouts
- [`as_tevector()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md)
  [`as_hstack_telayout()`](https://danigiro.github.io/FoReco/reference/tematrix_layouts.md)
  : Convert Between Horizon-stacked and Temporal Layouts
- [`balance_hierarchy()`](https://danigiro.github.io/FoReco/reference/balance_hierarchy.md)
  : Aggregation Matrix of a (Possibly) Unbalanced Hierarchy in Balanced
  Form
- [`commat()`](https://danigiro.github.io/FoReco/reference/commat.md)
  [`commat_index()`](https://danigiro.github.io/FoReco/reference/commat.md)
  : Commutation Matrix
- [`csboot()`](https://danigiro.github.io/FoReco/reference/csboot.md) :
  Cross-sectional Joint Block Bootstrap
- [`csprojmat()`](https://danigiro.github.io/FoReco/reference/csprojmat.md)
  : Projection Matrix for Optimal Combination Cross-sectional
  Reconciliation
- [`ctboot()`](https://danigiro.github.io/FoReco/reference/ctboot.md) :
  Cross-temporal Joint Block Bootstrap
- [`df2aggmat()`](https://danigiro.github.io/FoReco/reference/df2aggmat.md)
  : Cross-sectional Aggregation Matrix of a Dataframe
- [`lcmat()`](https://danigiro.github.io/FoReco/reference/lcmat.md) :
  Linear Combination (Aggregation) Matrix for a General Linearly
  Constrained Multiple Time Series
- [`res2matrix()`](https://danigiro.github.io/FoReco/reference/residuals.md)
  [`arrange_hres()`](https://danigiro.github.io/FoReco/reference/residuals.md)
  **\[deprecated\]** : One-Step and Multi-Step Residuals
- [`set_bounds()`](https://danigiro.github.io/FoReco/reference/set_bounds.md)
  : Set Bounds for Bounded Forecast Reconciliation
- [`teboot()`](https://danigiro.github.io/FoReco/reference/teboot.md) :
  Temporal Joint Block Bootstrap
- [`unbalance_hierarchy()`](https://danigiro.github.io/FoReco/reference/unbalance_hierarchy.md)
  : Aggregation Matrix of a Balanced Hierarchy in (Possibly) Unbalanced
  Form
- [`new_foreco_class()`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`summary(`*`<foreco>`*`)`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`print(`*`<summary_foreco>`*`)`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`print(`*`<foreco>`*`)`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`plot(`*`<foreco>`*`)`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`components(`*`<foreco>`*`)`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  [`drop_foreco_class()`](https://danigiro.github.io/FoReco/reference/foreco-class.md)
  : FoReco Reconciliation Class

## Package

- [`FoReco`](https://danigiro.github.io/FoReco/reference/FoReco-package.md)
  [`FoReco-package`](https://danigiro.github.io/FoReco/reference/FoReco-package.md)
  : FoReco: Forecast Reconciliation
