*! version 1.1.0 24jun2026
* rd2d_table_example_esttab.do
* Demonstrates rd2d_table for esttab/estout integration workflow
* Requires: rd2d package, estout package (ssc install estout)
version 16.0
clear all
set more off

// =====================================================================
// rd2d_table esttab Integration Example
// Shows how to convert rd2d multi-point results into esttab-compatible
// format using the rd2d_table post-estimation command.
// =====================================================================

// ─── Setup: locate package ado files ─────────────────────────────────────────
capture which rdbw2d
if _rc {
    foreach pkgpath in "rd2d-stata/ado" "rd2d-stata/pkg/ado" "ado" "pkg/ado" "../ado" "../pkg/ado" {
        capture confirm file "`pkgpath'/rdbw2d.ado"
        if !_rc {
            adopath ++ "`pkgpath'"
        }
    }
    foreach helppath in "rd2d-stata/help" "rd2d-stata/pkg/help" "help" "pkg/help" "../help" "../pkg/help" {
        capture confirm file "`helppath'/rdbw2d.sthlp"
        if !_rc {
            adopath ++ "`helppath'"
        }
    }
}

// ─── 1. Load data and run estimation ─────────────────────────────────────────

// Locate data files
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
gen byte d = t == "TRUE"
drop t

// Define evaluation points (5 points along boundary)
local atlist "0 50 0 27.5 0 0 22.5 0 47.5 0"

// Run rd2d with confidence bands (provides full covariance for rd2d_table)
rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
    scaleregul(3) stdvars level(95) cbands repp(500)

// ─── 2. Basic usage: display results as standard Stata table ─────────────────

display _n as txt "=== Example 1: Basic rd2d_table display ==="
rd2d_table

// ─── 3. Store and compare specifications ─────────────────────────────────────

display _n as txt "=== Example 2: Multi-specification comparison ==="

// Specification 1: bias-corrected (default)
rd2d_table
estimates store spec_rbc

// Re-run estimation (rd2d_table overwrites e(), so re-estimate)
rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
    scaleregul(3) stdvars level(95) cbands repp(500)

// Specification 2: conventional estimates
rd2d_table, estimate(p)
estimates store spec_conv

// Display comparison (requires esttab)
capture which esttab
if !_rc {
    display _n as txt "--- Comparison: RBC vs Conventional estimates ---"
    esttab spec_rbc spec_conv, se nostar ///
        mtitles("Bias-corrected" "Conventional") ///
        title("RD Treatment Effects: Specification Comparison")
}
else {
    display as txt "  esttab not installed; using estimates table instead"
    estimates table spec_rbc spec_conv, se
}

// ─── 4. Subset selection ─────────────────────────────────────────────────────

display _n as txt "=== Example 3: Subset of evaluation points ==="

// Re-run estimation
rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
    scaleregul(3) stdvars level(95) cbands repp(500)

// Post only points 1, 3, 5
rd2d_table, subset(1 3 5)

// ─── 5. Export to LaTeX and CSV (requires esttab) ────────────────────────────

display _n as txt "=== Example 4: File export ==="

// Re-estimate
rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
    scaleregul(3) stdvars level(95) cbands repp(500)

capture which esttab
if !_rc {
    // LaTeX export
    rd2d_table using "rd2d_esttab_output.tex", replace tex ///
        title("Boundary RD Treatment Effect Estimates")

    // CSV export
    rd2d_table

    // Re-estimate for CSV (rd2d_table overwrites e())
    rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
        bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
        scaleregul(3) stdvars level(95) cbands repp(500)

    rd2d_table using "rd2d_esttab_output.csv", replace csv
}
else {
    display as txt "  esttab not installed; skipping file export."
    display as txt "  Install via: ssc install estout"
}

// ─── 6. Advanced: lincom for cross-point inference ───────────────────────────

display _n as txt "=== Example 5: Cross-point inference with lincom ==="

// Re-estimate with cbands for full covariance
rd2d y x1 x2 d, at(`atlist') p(1) q(2) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(dpi) vce(hc1) bwcheck(52) masspoints(check) ///
    scaleregul(3) stdvars level(95) cbands repp(500)

rd2d_table

// Test difference between first and last evaluation point
// (only valid with full covariance from cbands)
capture noisily lincom _b[tau_0_50] - _b[tau_47d5_0]

// ─── Report ──────────────────────────────────────────────────────────────────
display _n as txt "=== rd2d_table esttab example complete ==="
estimates clear

di as txt "EXIT:0"
