version 16.0
clear all
set more off

* ==============================================================================
* rd2d_sensitivity_example.do
* Demonstrates bandwidth sensitivity analysis using rd2d_bwsens
* ==============================================================================

* --- Setup: adopath ---
capture which rd2d_bwsens
if _rc {
    foreach pkgpath in "rd2d-stata/ado" "rd2d-stata/pkg/ado" "ado" "pkg/ado" "../ado" "../pkg/ado" {
        capture confirm file "`pkgpath'/rd2d_bwsens.ado"
        if !_rc {
            adopath ++ "`pkgpath'"
        }
    }
    foreach helppath in "rd2d-stata/help" "rd2d-stata/pkg/help" "help" "pkg/help" "../help" "../pkg/help" {
        capture confirm file "`helppath'/rd2d_bwsens.sthlp"
        if !_rc {
            adopath ++ "`helppath'"
        }
    }
}

* ==============================================================================
* Part 1: Location-based RD (rd2d) sensitivity analysis
* ==============================================================================

* Generate simulated 2D RD data
set obs 500
generate double row = mod(_n - 1, 25) - 12
generate double col = floor((_n - 1) / 25) - 9.5
generate double x1 = row / 5
generate double x2 = col / 5
generate byte treat = x1 >= 0
generate double y = 1.0 + 0.5*x1 - 0.3*x2 + 1.5*treat ///
    + 0.1*x1*x2 + 0.05*mod(_n, 7) + invnormal(uniform())*0.3

* Step 1: Run the base estimation
di as txt _newline "=== Step 1: Base rd2d estimation ==="
rd2d y x1 x2 treat, at(0 0) h(1.5) p(1) q(2) kernel(triangular) ///
    ktype(prod) bwcheck(20)

* Step 2: Bandwidth sensitivity with default grid
di as txt _newline "=== Step 2: Bandwidth sensitivity (default grid) ==="
rd2d_bwsens

* Display the returned matrix
di as txt _newline "Returned sensitivity matrix:"
matrix list r(sens_results), format(%9.4f)

* Step 3: Custom finer grid around baseline
di as txt _newline "=== Step 3: Fine-resolution sensitivity ==="
rd2d_bwsens, grid(0.6 0.7 0.8 0.9 0.95 1.0 1.05 1.1 1.2 1.3)

* ==============================================================================
* Part 2: Distance-based RD (rd2d_dist) sensitivity analysis
* ==============================================================================

* Generate distance data
drop _all
set obs 400
generate double dist1 = invnormal(uniform()) * 2
generate double dist2 = invnormal(uniform()) * 2
generate byte treat = (dist1 > 0)
generate double y = 0.5 + 0.3*dist1 - 0.2*dist2 + 1.2*treat ///
    + invnormal(uniform())*0.4

* Step 4: Run rd2d_dist estimation
di as txt _newline "=== Step 4: Base rd2d_dist estimation ==="
rd2d_dist y dist1 dist2, h(1.0) p(1) q(2) kernel(epanechnikov)

* Step 5: Sensitivity analysis for distance-based RD
di as txt _newline "=== Step 5: Bandwidth sensitivity for rd2d_dist ==="
rd2d_bwsens

* ==============================================================================
* Part 3: Interpreting results
* ==============================================================================

di as txt _newline "=== Interpretation Guide ==="
di as txt "1. Look at the Est.q column across multipliers."
di as txt "   Stable values indicate robustness to bandwidth choice."
di as txt ""
di as txt "2. Check whether confidence intervals overlap across the grid."
di as txt "   Non-overlapping CIs at adjacent multipliers suggest fragility."
di as txt ""
di as txt "3. The baseline row (marked with *) corresponds to your original"
di as txt "   estimation. Deviations from this are informative about sensitivity."
di as txt ""
di as txt "4. For publication: report the MSE-optimal bandwidth estimate."
di as txt "   The sensitivity table belongs in an appendix or supplementary materials."

di as txt _newline "=== Example complete ==="
