# rd2d

**Boundary Regression Discontinuity Designs in Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![Version: 1.1.0](https://img.shields.io/badge/Version-1.1.0-green.svg)](CHANGELOG.md)

## Overview

`rd2d` implements local polynomial estimation and inference for boundary
discontinuity designs in Stata. Standard RD designs use a scalar running
variable with a single cutoff. Boundary RD extends this to settings where a
bivariate score (X1, X2) determines treatment assignment through a
one-dimensional boundary curve, as in geographic RD or multi-score threshold
designs. The causal estimand is the boundary average treatment effect curve:

> τ(b) = E[Y(1) − Y(0) | X = b],  b ∈ B

The package supports two approaches to estimation, each motivated by different
data availability and geometric constraints:

- **Location-based** (`rd2d`, `rdbw2d`): uses bivariate local polynomial
  regression directly on (X1, X2) relative to boundary points. This approach
  provides MSE-optimal bandwidth selection and remains valid regardless of
  boundary geometry, including at kinks.
- **Distance-based** (`rd2d_dist`, `rdbw2d_dist`): transforms the bivariate
  score into a signed distance to the boundary. This simplifies estimation to
  univariate local polynomial regression with cutoff at zero, but requires
  caution near boundary kinks where the misspecification bias is irreducible at
  rate O(h) regardless of polynomial order.

Both approaches offer pointwise and uniform robust bias-corrected inference.
Post-estimation commands handle visualization, aggregation, bandwidth
sensitivity, table construction, and diagnostics.

## Installation

```stata
net install rd2d, from("https://raw.githubusercontent.com/gorgeousfish/rd2d/main") replace
net get rd2d, from("https://raw.githubusercontent.com/gorgeousfish/rd2d/main") replace
```

`net install` places ado and help files on the adopath. `net get` retrieves
example do-files into the current working directory.

## Quick Start

### Location-Based Estimation

```stata
. rd2d y x1 x2 t, at(0 0)

Location RD estimation
-------------------------------------------------------------------------------
  Evaluation points:         1    Observations:     20000
  VCE: hc1    RBC: on    Side: two
  Bandwidth source: automatic    Kernel: triangular    Std. Vars: On (default)
-------------------------------------------------------------------------------
     Point      Est.q      SE.q    CI.low   CI.high    CB.low   CB.high
-------------------------------------------------------------------------------
       at1      .6635    .07854     .5096     .8175     .5039     .8231
-------------------------------------------------------------------------------
  Uniform bands: on    Critical value:     2.032    Repetitions:      1000
```

Multiple boundary points with uniform confidence bands:

```stata
. rd2d y x1 x2 t, at(0 0 0.5 0.5 -0.5 -0.5)

Location RD estimation
-------------------------------------------------------------------------------
  Evaluation points:         3    Observations:     20000
  VCE: hc1    RBC: on    Side: two
  Bandwidth source: automatic    Kernel: triangular    Std. Vars: On (default)
-------------------------------------------------------------------------------
     Point      Est.q      SE.q    CI.low   CI.high    CB.low   CB.high
-------------------------------------------------------------------------------
       at1      .6635    .07854     .5096     .8175      .497     .8301
       at2      .6764    .06205     .5548      .798     .5448      .808
       at3      .6444     .1034     .4419      .847     .4252     .8636
-------------------------------------------------------------------------------
  Uniform bands: on    Critical value:     2.121    Repetitions:      1000
```

### Distance-Based Estimation

```stata
. rd2d_dist y D

Distance RD estimation
-------------------------------------------------------------------------------
  Evaluation points:         1    Observations:     20000
  VCE: hc1    RBC: on    Side: two
  Bandwidth source: automatic    Kernel: triangular    Kink: off
-------------------------------------------------------------------------------
     Point      Est.q      SE.q    CI.low   CI.high    CB.low   CB.high
-------------------------------------------------------------------------------
         D      .6954    .09922     .5009     .8898     .4926     .8982
-------------------------------------------------------------------------------
  Uniform bands: on    Critical value:     2.044    Repetitions:      1000
```

### Post-Estimation

```stata
* Visualize treatment effects along the boundary
rd2d_plot

* Structured summary with diagnostics
rd2d_summary

* Aggregate effects (weighted boundary average treatment effect)
rd2d_aggregate, method(wbate)

* Bandwidth sensitivity
rd2d_bwsens

* Publication table via esttab
rd2d_table
```

## Commands

### Estimation

| Command         | Purpose                                              | Minimal Syntax                 |
| :-------------- | :--------------------------------------------------- | :----------------------------- |
| `rdbw2d`      | MSE-optimal bandwidth selection (location)           | `yvar x1 x2 tvar, at(b1 b2)` |
| `rd2d`        | Treatment effect estimation and inference (location) | `yvar x1 x2 tvar, at(b1 b2)` |
| `rdbw2d_dist` | Bandwidth selection (distance)                       | `yvar distvar`               |
| `rd2d_dist`   | Treatment effect estimation and inference (distance) | `yvar distvar`               |

### Post-Estimation

| Command              | Purpose                                              |
| :------------------- | :--------------------------------------------------- |
| `rd2d_plot`        | Effect plots and boundary heat maps                  |
| `rd2d_summary`     | Structured result display with bandwidth diagnostics |
| `rd2d_aggregate`   | Aggregation: WBATE, AATE, LBATE                      |
| `rd2d_table`       | Publication table construction (esttab-compatible)   |
| `rd2d_diagnostics` | Fit conditioning, support, and covariance checks     |
| `rd2d_bwsens`      | Bandwidth sensitivity analysis                       |

## Key Options

| Option           | Default        | Description                                                             |
| :--------------- | :------------- | :---------------------------------------------------------------------- |
| `p(#)`         | 1              | Local polynomial order for point estimation                             |
| `q(#)`         | p+1            | Polynomial order for bias correction                                    |
| `kernel()`     | `triangular` | Kernel:`triangular`, `uniform`, `epanechnikov`, `gaussian`      |
| `bwselect()`   | `mserd`      | Bandwidth selector:`mserd`, `imserd`, `msetwo`, `imsetwo`       |
| `vce()`        | `hc1`        | Variance estimator:`hc0`–`hc3`; cluster-robust for `hc0`/`hc1` |
| `cluster()`    | —             | Cluster variable for cluster-robust inference                           |
| `cbands`       | on             | Gaussian-simulation uniform confidence bands                            |
| `kink(on)`     | off            | Nonsmooth-boundary convention (distance commands)                       |
| `stdvars`      | off            | Standardize coordinates for bandwidth selection                         |
| `masspoints()` | `check`      | Mass-point handling:`check`, `adjust`, `off`                      |
| `tangvec()`    | —             | Directional derivative target (location commands)                       |
| `scaleregul()` | 3              | Regularization scale for bandwidth selection                            |

## Stored Results

Estimation commands post results in `e()`:

| Object             | Content                                                                |
| :----------------- | :--------------------------------------------------------------------- |
| `e(results)`     | Point estimates, SEs, pointwise CIs, and band endpoints per target     |
| `e(bws)`         | Final bandwidths (h01, h02, h11, h12 or h0, h1) and local sample sizes |
| `e(diagnostics)` | Fit conditioning and generalized-inverse flags                         |
| `e(cov_q)`       | Q-order covariance matrix for uniform inference                        |
| `e(cb_crit)`     | Simulated critical value for confidence bands                          |
| `e(masspoints)`  | Mass-point support diagnostics                                         |
| `e(b)`, `e(V)` | Stata estimation conventions for postestimation                        |

Bandwidth selectors post results in `r()`:

| Object            | Content                                      |
| :---------------- | :------------------------------------------- |
| `r(bws)`        | Selected side-specific bandwidths            |
| `r(mseconsts)`  | MSE expansion constants used by the selector |
| `r(masspoints)` | Support diagnostics                          |

## Scope and Limitations

The package requires:

- A continuous bivariate running variable or a user-constructed signed-distance
  score.
- Known boundary location (specified via `at()` for location commands).
- Stata 16 or newer. No R runtime dependency.

Diagnostics are reporting aids, not identification proofs. They help a table
script decide what to disclose or rerun before publishing a row, but they do not
replace the assumptions, data provenance, or design-specific sensitivity checks
required for an empirical application.

## References

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025). Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods. *arXiv preprint arXiv:2505.05670*.

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025). rd2d: Causal inference in boundary discontinuity designs. *arXiv preprint arXiv:2505.07989*.

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2026). Estimation and inference in boundary discontinuity designs: Distance-based methods. *Journal of Econometrics*, 256, 106266.

### BibTeX

```bibtex
@unpublished{cattaneo2025boundaryrd,
  title={Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods},
  author={Cattaneo, Matias D. and Titiunik, Rocio and Yu, Ruiqi},
  year={2025},
  note={arXiv:2505.05670}
}

@unpublished{cattaneo2025rd2d,
  title={{rd2d}: Causal Inference in Boundary Discontinuity Designs},
  author={Cattaneo, Matias D. and Titiunik, Rocio and Yu, Ruiqi},
  year={2025},
  note={arXiv:2505.07989}
}

@article{cattaneo2026distance,
  title={Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods},
  author={Cattaneo, Matias D. and Titiunik, Rocio and Yu, Ruiqi},
  journal={Journal of Econometrics},
  volume={256},
  pages={106266},
  year={2026}
}
```

## Authors

**Methodology:**

- Matias D. Cattaneo, Princeton University
- Rocio Titiunik, Princeton University
- Ruiqi Yu, Princeton University

**Implementation:**

- Xuanyu Cai, City University of Macau
- Wenli Xu, City University of Macau
