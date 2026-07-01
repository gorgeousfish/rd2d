*! version 1.1.0 01jul2026
*! rd2d_tutorial_workflow.do
*! Complete research workflow tutorial for the rd2d package
*! Demonstrates: data → bandwidth → estimation → diagnostics → sensitivity → tables → figures
*!
*! Reference: Cattaneo, Titiunik, and Yu (2025). "rd2d: Boundary Regression
*!   Discontinuity Designs in Two-Dimensional Running Variable Spaces."
*!   Working Paper, arXiv:2502.19788.
*!
*! NOTE: This tutorial requires the bundled data under rd2d-stata/data/.
*!   Run from the examples/ or rd2d-stata/ directory, or from the repository root.

version 16.0
clear all
set more off

* --- Setup: adopath resolution ---
capture which rd2d
if _rc {
    foreach pkgpath in "rd2d-stata/ado" "rd2d-stata/pkg/ado" "ado" "pkg/ado" "../ado" "../pkg/ado" {
        capture confirm file "`pkgpath'/rd2d.ado"
        if !_rc {
            adopath ++ "`pkgpath'"
        }
    }
    foreach helppath in "rd2d-stata/help" "rd2d-stata/pkg/help" "help" "pkg/help" "../help" "../pkg/help" {
        capture confirm file "`helppath'/rd2d.sthlp"
        if !_rc {
            adopath ++ "`helppath'"
        }
    }
}

* ============================================================================
* 1. DATA PREPARATION
* ============================================================================
* The boundary RD design requires: an outcome y, two running variables
* (x1, x2) defining location relative to a treatment boundary, and a
* binary treatment indicator. The bundled dataset contains 400 observations
* from a simulated geographic boundary design.

* --- Locate the data file across possible working directories ---
local datadir ""
foreach dpath in "rd2d-stata/data" "data" "../data" "../rd2d-stata/data" {
    capture confirm file "`dpath'/data_rd2d.csv"
    if !_rc {
        local datadir "`dpath'"
        continue, break
    }
}
if "`datadir'" == "" {
    display as err "Cannot locate data_rd2d.csv; run from repository root or examples/ directory"
    exit 601
}

import delimited "`datadir'/data_rd2d.csv", varnames(1) clear asdouble

* The treatment variable is imported as string ("TRUE"/"FALSE"); recode to numeric
gen byte d = (t == "TRUE")
drop t

* Rename variables to standard Stata conventions
rename x_1 x1
rename x_2 x2

* Inspect the data structure
describe
summarize y x1 x2 d

* Treatment frequency: approximately balanced by design
tabulate d

* ============================================================================
* 2. BANDWIDTH SELECTION
* ============================================================================
* The MSE-optimal bandwidth minimizes integrated mean squared error of the
* local polynomial estimator at the boundary (Cattaneo, Titiunik, Yu 2025,
* Section 3). Two methods are available:
*   DPI  - direct plug-in estimator (preferred for inference)
*   ROT  - rule-of-thumb (faster, used as initial guess)
*
* Bandwidth is selected independently for each dimension and each side of
* the boundary, yielding four components: h01, h02 (control), h11, h12 (treated).

* --- Method 1: DPI (direct plug-in) ---
di as txt _newline "=== DPI Bandwidth Selection ==="
rdbw2d y x1 x2 d, at(0 0) method(dpi) kernel(triangular) ktype(prod) bwcheck(20)
matrix BW_dpi = r(bws)
display "  DPI bandwidths: h01=" %6.3f BW_dpi[1,3] "  h02=" %6.3f BW_dpi[1,4]

* --- Method 2: ROT (rule-of-thumb) for comparison ---
di as txt _newline "=== ROT Bandwidth Selection ==="
rdbw2d y x1 x2 d, at(0 0) method(rot) kernel(triangular) ktype(prod) bwcheck(20)
matrix BW_rot = r(bws)
display "  ROT bandwidths: h01=" %6.3f BW_rot[1,3] "  h02=" %6.3f BW_rot[1,4]

* ============================================================================
* 3. MAIN ESTIMATION: SINGLE BOUNDARY POINT
* ============================================================================
* Estimate the treatment effect at a single boundary point (0, 0).
* The local polynomial estimator uses:
*   p=1 (local linear for point estimation)
*   q=2 (local quadratic for bias correction, enabling robust inference)
*   Triangular product kernel (adapts to boundary geometry)
*
* The cbands option computes uniform confidence bands across evaluation
* points and stores the cross-point covariance matrix, which is required
* for aggregation (Section 4 of CTY 2025).

di as txt _newline "=== Single-Point Estimation ==="
rd2d y x1 x2 d, at(0 0) p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(20) ///
    masspoints(check) scaleregul(3) level(95)

* Structured display of results
rd2d_summary

* Bandwidth diagnostics: effective sample sizes within each bandwidth window
rd2d_summary, output(bw)

* ============================================================================
* 4. MULTI-POINT ESTIMATION WITH UNIFORM INFERENCE
* ============================================================================
* In boundary RD, treatment effects may vary along the boundary (spatial
* heterogeneity). We evaluate at multiple points spanning the boundary to
* characterize this heterogeneity. The cbands option provides simultaneous
* coverage across all evaluation points (Theorem 5, CTY 2025).

di as txt _newline "=== Multi-Point Estimation with Confidence Bands ==="
rd2d y x1 x2 d, at(0 0 0 15 0 27.5 0 40 10 0 22.5 0 35 0) ///
    p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(20) ///
    masspoints(check) scaleregul(3) level(95) cbands repp(500)

* Display results with uniform confidence bands
rd2d_summary, cbuniform

* ============================================================================
* 5. DIAGNOSTICS
* ============================================================================
* Numerical diagnostics assess local polynomial fit quality: matrix
* condition numbers, rank deficiency, pseudo-inverse fallback usage,
* and mass-point concentration. These checks guard against numerical
* instability from sparse boundary neighborhoods.

di as txt _newline "=== Estimation Diagnostics ==="

* Compact summary: one row per evaluation point with status flags
rd2d_diagnostics, output(summary)

* Detailed view: expand all four fit combinations per point
rd2d_diagnostics, output(full)

* Show only problematic points with remedial suggestions
rd2d_diagnostics, output(warnings)

* ============================================================================
* 6. BANDWIDTH SENSITIVITY ANALYSIS
* ============================================================================
* Bandwidth sensitivity assesses whether conclusions are robust to the
* specific bandwidth chosen. We scale the baseline bandwidth by multipliers
* {0.5, ..., 2.0} and re-estimate. Stable estimates across the grid
* indicate that results are not artifacts of a particular bandwidth.

di as txt _newline "=== Bandwidth Sensitivity ==="

* Default grid (0.5 0.75 0.9 1.0 1.1 1.25 1.5 2.0)
rd2d_bwsens, subset(1)

* Finer grid near baseline for the first evaluation point
rd2d_bwsens, grid(0.7 0.8 0.9 0.95 1.0 1.05 1.1 1.2 1.3) subset(1)

* Generate a sensitivity plot
rd2d_bwsens, subset(1) plot

* ============================================================================
* 7. VISUALIZATION
* ============================================================================
* Two complementary displays:
*   effect: coefficient plot showing estimates with uncertainty intervals
*           ordered by evaluation point index (for spotting heterogeneity)
*   heat:   bubble scatter in (b1, b2) space with marker area proportional
*           to |tau| (for visualizing spatial treatment effect patterns)

di as txt _newline "=== Visualization ==="

* Effect plot: point estimates with both CI and uniform confidence bands
rd2d_plot, type(effect) interval(both)

* Heat map: spatial distribution of treatment effects along boundary
rd2d_plot, type(heat)

* ============================================================================
* 8. PUBLICATION TABLE GENERATION
* ============================================================================
* rd2d_table repacks multi-point results into standard eclass form,
* enabling seamless integration with esttab/estout. This produces
* tables conforming to Stata Journal formatting requirements.

di as txt _newline "=== Table Generation ==="

* Display formatted results via ereturn display
rd2d_table

* Export to LaTeX (requires estout package)
capture rd2d_table using "rd2d_tutorial_table.tex", replace tex ///
    title("Boundary RD Treatment Effect Estimates")

* ============================================================================
* 9. AGGREGATION: BOUNDARY-WIDE TREATMENT EFFECTS
* ============================================================================
* Point-by-point estimates characterize heterogeneity; aggregated measures
* summarize the overall treatment effect along the boundary.
*
*   AATE:  Arithmetic Average Treatment Effect (equal weights; Theorem 6)
*   WBATE: Weighted Boundary ATE (user-supplied weights)
*   LBATE: Largest Boundary ATE (max estimate with uniform band CI)

di as txt _newline "=== Aggregated Treatment Effects ==="

* AATE: unweighted average across all evaluation points
rd2d_aggregate, method(aate) level(95)
display "  AATE = " %9.4f r(estimate) "  SE = " %9.4f r(se) ///
    "  95% CI = [" %9.4f r(ci_lower) ", " %9.4f r(ci_upper) "]"

* LBATE: largest boundary treatment effect with uniform-band CI
rd2d_aggregate, method(lbate) level(95)
display "  LBATE = " %9.4f r(estimate) ///
    "  95% CI = [" %9.4f r(ci_lower) ", " %9.4f r(ci_upper) "]"

* ============================================================================
* 10. DISTANCE-BASED ESTIMATION (rd2d_dist)
* ============================================================================
* An alternative approach uses signed distances to the boundary as running
* variables. Negative distances denote control; nonnegative denote treated.
* This is appropriate when the boundary geometry is complex or when the
* researcher has a natural distance metric (Cattaneo, Titiunik, Yu 2025,
* Section 5).

di as txt _newline "=== Distance-Based Estimation ==="

* Construct signed distances (positive = treated side)
gen double dist = cond(d == 1, sqrt(x1^2 + x2^2), -sqrt(x1^2 + x2^2))

* Estimate via signed-distance local polynomial
rd2d_dist y dist, p(1) q(2) kernel(triangular) vce(hc1) level(95)

* Display results
rd2d_summary

* Diagnostics for distance-based estimation
rd2d_diagnostics, output(summary)

* ============================================================================
* END OF TUTORIAL
* ============================================================================
di _n as txt "{hline 64}"
di as txt " Tutorial complete."
di as txt " Commands demonstrated: rdbw2d, rd2d, rd2d_summary, rd2d_diagnostics,"
di as txt "   rd2d_bwsens, rd2d_plot, rd2d_table, rd2d_aggregate, rd2d_dist"
di as txt " See {help rd2d} for full documentation."
di as txt "{hline 64}"
