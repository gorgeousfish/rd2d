*! version 1.1.0 24jun2026
*! rd2d_plot: Post-estimation plotting for rd2d and rd2d_dist
*! Authors: Cattaneo, Titiunik, Yu

program define rd2d_plot
    version 16.0
    syntax [, Type(string) ESTimate(string) INTerval(string) ///
             SUBset(numlist integer >0) LABEL ///
             SAVing(string) TITLE(string) SUBTitle(string) ///
             SCHEME(string) NAME(string) REPLACE ///
             NODRaw]

    // ===================================================================
    // 1. CHECK e() RESULTS EXIST
    // ===================================================================
    local ecmd "`e(cmd)'"
    if !inlist("`ecmd'", "rd2d", "rd2d_dist") {
        di as err "rd2d_plot requires estimation results from {bf:rd2d} or {bf:rd2d_dist}"
        exit 198
    }

    // ===================================================================
    // 2. PARAMETER DEFAULTS AND VALIDATION
    // ===================================================================
    // type
    local type = lower("`type'")
    if ("`type'" == "") local type "effect"
    if !inlist("`type'", "effect", "heat", "heterogeneity", "hetero") {
        di as err "type() must be {bf:effect}, {bf:heat}, or {bf:heterogeneity}; got '`type''"
        exit 198
    }
    if "`type'" == "hetero" local type "heterogeneity"

    // estimate
    local estimate = lower("`estimate'")
    if ("`estimate'" == "") local estimate "q"
    if !inlist("`estimate'", "q", "p") {
        di as err "estimate() must be {bf:q} (bias-corrected) or {bf:p} (conventional); got '`estimate''"
        exit 198
    }

    // interval
    local interval = lower("`interval'")
    if ("`interval'" == "") local interval "ci"
    if !inlist("`interval'", "ci", "cb", "both", "none") {
        di as err "interval() must be {bf:ci}, {bf:cb}, {bf:both}, or {bf:none}; got '`interval''"
        exit 198
    }

    // CB availability check
    local cbands_on "`e(cbands)'"
    if inlist("`interval'", "cb", "both") & "`cbands_on'" != "on" {
        di as err "confidence bands requested but {bf:cbands} was not specified in estimation"
        di as err "rerun {bf:`ecmd'} with the {bf:cbands} option, or use interval(ci) or interval(none)"
        exit 198
    }

    // CB not available for estimate(p)
    if "`estimate'" == "p" & inlist("`interval'", "cb", "both") {
        di as txt "(note: confidence bands are defined only for bias-corrected estimates;"
        di as txt " interval changed to {bf:ci})"
        if "`interval'" == "cb" local interval "ci"
        if "`interval'" == "both" local interval "ci"
    }

    // heat map requires boundary coordinates
    if "`type'" == "heat" & "`ecmd'" == "rd2d_dist" {
        // rd2d_dist stores b1/b2 as columns 1-2 of e(results), but for
        // distance-based estimation these may be degenerate (all zeros or
        // missing) if the user supplied only distance variables
        tempname _chkmat
        matrix `_chkmat' = e(results)
        local _nchk = rowsof(`_chkmat')
        local _has_bdy 0
        forvalues _i = 1/`_nchk' {
            if (`_chkmat'[`_i', 1] < . & `_chkmat'[`_i', 2] < .) ///
                & (`_chkmat'[`_i', 1] != 0 | `_chkmat'[`_i', 2] != 0) {
                local _has_bdy 1
                continue, break
            }
        }
        if !`_has_bdy' {
            di as err "type(heat) requires boundary coordinates (b1, b2),"
            di as err "but {bf:rd2d_dist} results contain no location information"
            exit 198
        }
    }

    // neval and subset validation
    local neval = e(neval)
    if "`subset'" != "" {
        foreach _s of numlist `subset' {
            if `_s' > `neval' {
                di as err "subset() value `_s' exceeds number of evaluation points (`neval')"
                exit 198
            }
        }
    }

    // ===================================================================
    // 3. EXTRACT DATA FROM e() INTO TEMPORARY VARIABLES
    // ===================================================================
    tempname resmat
    matrix `resmat' = e(results)
    local level = e(level)
    local side "`e(side)'"

    // Determine which rows to plot
    if "`subset'" != "" {
        local plotrows "`subset'"
        local nplot : word count `subset'
    }
    else {
        local plotrows ""
        forvalues _i = 1/`neval' {
            local plotrows "`plotrows' `_i'"
        }
        local nplot `neval'
    }

    // ===================================================================
    // 4. COMPUTE z-CRITICAL VALUE FOR estimate(p) CI RECALCULATION
    // ===================================================================
    tempname zcrit
    if "`side'" == "two" {
        scalar `zcrit' = invnormal((100 + `level') / 200)
    }
    else {
        scalar `zcrit' = invnormal(`level' / 100)
    }

    // ===================================================================
    // 5. BUILD PLOT DATA IN TEMPORARY FRAME
    // ===================================================================
    preserve
    quietly {
        clear
        set obs `nplot'

        // Position variable (sequential index for X axis)
        tempvar plot_pos orig_idx
        gen int `plot_pos' = _n
        gen int `orig_idx' = .

        // Data variables
        tempvar estimate_v ci_lo ci_hi cb_lo cb_hi b1v b2v
        gen double `estimate_v' = .
        gen double `ci_lo' = .
        gen double `ci_hi' = .
        gen double `cb_lo' = .
        gen double `cb_hi' = .
        gen double `b1v' = .
        gen double `b2v' = .

        // Fill in data from results matrix
        local _obs 0
        foreach _row of numlist `plotrows' {
            local _obs = `_obs' + 1
            replace `orig_idx' = `_row' in `_obs'
            replace `b1v' = `resmat'[`_row', 1] in `_obs'
            replace `b2v' = `resmat'[`_row', 2] in `_obs'

            if "`estimate'" == "q" {
                // Bias-corrected
                replace `estimate_v' = `resmat'[`_row', 5] in `_obs'
                replace `ci_lo' = `resmat'[`_row', 9] in `_obs'
                replace `ci_hi' = `resmat'[`_row', 10] in `_obs'
                replace `cb_lo' = `resmat'[`_row', 11] in `_obs'
                replace `cb_hi' = `resmat'[`_row', 12] in `_obs'
            }
            else {
                // Conventional: recalculate CI from Est_p and Se_p
                local _est = `resmat'[`_row', 3]
                local _se  = `resmat'[`_row', 4]
                replace `estimate_v' = `_est' in `_obs'
                if "`side'" == "two" {
                    replace `ci_lo' = `_est' - `zcrit' * `_se' in `_obs'
                    replace `ci_hi' = `_est' + `zcrit' * `_se' in `_obs'
                }
                else if "`side'" == "left" {
                    replace `ci_lo' = . in `_obs'
                    replace `ci_hi' = `_est' + `zcrit' * `_se' in `_obs'
                }
                else {
                    // right
                    replace `ci_lo' = `_est' - `zcrit' * `_se' in `_obs'
                    replace `ci_hi' = . in `_obs'
                }
                // CB not applicable for estimate(p)
                replace `cb_lo' = . in `_obs'
                replace `cb_hi' = . in `_obs'
            }
        }

        // Handle one-sided intervals: replace infinite bounds with missing
        // (Stata stores ±maxdouble for unbounded; we use . for open ends)
        if "`side'" == "left" & "`estimate'" == "q" {
            replace `ci_lo' = . if `ci_lo' <= -c(maxdouble)/2
            replace `cb_lo' = . if `cb_lo' <= -c(maxdouble)/2
        }
        if "`side'" == "right" & "`estimate'" == "q" {
            replace `ci_hi' = . if `ci_hi' >= c(maxdouble)/2
            replace `cb_hi' = . if `cb_hi' >= c(maxdouble)/2
        }

        // Generate value labels for X axis
        if "`label'" != "" | "`type'" == "effect" {
            // Label X axis with original evaluation-point indices
            tempname xlbl
            label define `xlbl' 0 " ", add
            local _obs 0
            foreach _row of numlist `plotrows' {
                local _obs = `_obs' + 1
                label define `xlbl' `_obs' "`_row'", add
            }
            label values `plot_pos' `xlbl'
        }
    }

    // ===================================================================
    // 6. CONSTRUCT GRAPH OPTIONS
    // ===================================================================
    // Saving option
    local saving_opt ""
    if `"`saving'"' != "" {
        if "`replace'" != "" {
            local saving_opt `"saving(`saving', replace)"'
        }
        else {
            local saving_opt `"saving(`saving')"'
        }
    }

    // Name option
    local name_opt ""
    if "`name'" != "" {
        if "`replace'" != "" {
            local name_opt `"name(`name', replace)"'
        }
        else {
            local name_opt `"name(`name')"'
        }
    }

    // Scheme
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt `"scheme(`scheme')"'
    }

    // Nodraw
    local nodraw_opt ""
    if "`nodraw'" != "" {
        local nodraw_opt "nodraw"
    }

    // Title defaults
    if `"`title'"' == "" {
        if "`type'" == "effect" {
            if "`estimate'" == "q" {
                local title "Bias-corrected treatment effect estimates"
            }
            else {
                local title "Conventional treatment effect estimates"
            }
        }
        else if "`type'" == "heterogeneity" {
            local title "Treatment Effect Heterogeneity along Boundary"
        }
        else {
            local title "Treatment effects at boundary"
        }
    }

    // Subtitle default
    if `"`subtitle'"' == "" {
        local _krnl "`e(kernel)'"
        local _vce  "`e(vce)'"
        local subtitle "Kernel: `_krnl', VCE: `_vce', Level: `level'%"
    }

    // Y-axis title
    local ytitle "Treatment effect"
    if "`estimate'" == "q" {
        local ytitle "Treatment effect (bias-corrected)"
    }
    else {
        local ytitle "Treatment effect (conventional)"
    }

    // ===================================================================
    // 7. DRAW PLOT
    // ===================================================================
    if "`type'" == "effect" {
        // Determine x-axis range
        local xmin 0.5
        local xmax = `nplot' + 0.5

        // Build xlabel list
        local xlabels ""
        forvalues _i = 1/`nplot' {
            local xlabels "`xlabels' `_i'"
        }

        // Construct twoway layers
        local layers ""

        // Layer 0: zero reference line
        local layers `"`layers' (function y=0, range(`xmin' `xmax') lcolor(gs10) lpattern(dash) lwidth(thin))"'

        // Layer 1: Confidence bands (rarea) if requested
        if inlist("`interval'", "cb", "both") {
            local layers `"`layers' (rarea `cb_lo' `cb_hi' `plot_pos', color("0 114 178%30") lwidth(none) fintensity(100))"'
        }

        // Layer 2: Confidence intervals (rcap) if requested
        if inlist("`interval'", "ci", "both") {
            if "`side'" == "two" {
                local layers `"`layers' (rcap `ci_lo' `ci_hi' `plot_pos', lcolor("0 114 178") lwidth(medthin))"'
            }
            else if "`side'" == "left" {
                // One-sided left: lower bound is -inf, use rspike from estimate down
                // and show upper cap
                local layers `"`layers' (rspike `estimate_v' `ci_hi' `plot_pos', lcolor("0 114 178") lwidth(medthin))"'
            }
            else {
                // One-sided right: upper bound is +inf
                local layers `"`layers' (rspike `ci_lo' `estimate_v' `plot_pos', lcolor("0 114 178") lwidth(medthin))"'
            }
        }

        // Layer 3: Point estimates
        local layers `"`layers' (scatter `estimate_v' `plot_pos', mcolor("0 114 178") msymbol(O) msize(medium))"'

        // Assemble graph command
        twoway `layers' ///
            , xlabel(`xlabels', valuelabel) ///
              xscale(range(`xmin' `xmax')) ///
              xtitle("Evaluation point") ///
              ytitle("`ytitle'") ///
              title(`"`title'"') ///
              subtitle(`"`subtitle'"') ///
              legend(off) ///
              `scheme_opt' `saving_opt' `name_opt' `nodraw_opt'
    }
    else if "`type'" == "heterogeneity" {
        // ===============================================================
        // HETEROGENEITY PLOT: Treatment effect along boundary
        // X-axis: evaluation point index with coordinate labels
        // Y-axis: tau_hat with confidence interval bands
        // ===============================================================
        local xmin 0.5
        local xmax = `nplot' + 0.5

        // Build xlabel with boundary coordinate labels
        local xlabels ""
        local _obs 0
        foreach _row of numlist `plotrows' {
            local _obs = `_obs' + 1
            local _b1 = string(`resmat'[`_row', 1], "%5.2f")
            local _b2 = string(`resmat'[`_row', 2], "%5.2f")
            local xlabels `"`xlabels' `_obs' "(`_b1',`_b2')""'
        }

        // Construct twoway layers
        local layers ""

        // Layer 0: zero reference line
        local layers `"`layers' (function y=0, range(`xmin' `xmax') lcolor(gs10) lpattern(dash) lwidth(thin))"'

        // Layer 1: CI band as shaded area (rarea)
        if inlist("`interval'", "ci", "both") {
            local layers `"`layers' (rarea `ci_lo' `ci_hi' `plot_pos', color("0 114 178%20") lwidth(none))"'
        }

        // Layer 2: CB band as lighter area if available
        if inlist("`interval'", "cb", "both") {
            local layers `"`layers' (rarea `cb_lo' `cb_hi' `plot_pos', color("213 94 0%15") lwidth(none))"'
        }

        // Layer 3: Connected point estimates
        local layers `"`layers' (connected `estimate_v' `plot_pos', lcolor("0 114 178") lpattern(solid) lwidth(medthick) mcolor("0 114 178") msymbol(O) msize(medium))"'

        // Legend
        local legend_items ""
        local legend_idx 2
        if inlist("`interval'", "ci", "both") {
            local legend_items `"`legend_items' `legend_idx' "`level'% CI""'
            local legend_idx = `legend_idx' + 1
        }
        if inlist("`interval'", "cb", "both") {
            local legend_items `"`legend_items' `legend_idx' "Uniform CB""'
            local legend_idx = `legend_idx' + 1
        }
        local legend_items `"`legend_items' `legend_idx' "Point estimate""'
        local legend_opt `"legend(order(`legend_items') position(6) rows(1))"'

        twoway `layers' ///
            , xlabel(`xlabels', angle(45) labsize(small)) ///
              xscale(range(`xmin' `xmax')) ///
              xtitle("Boundary point (b{subscript:1}, b{subscript:2})") ///
              ytitle("`ytitle'") ///
              title(`"`title'"') ///
              subtitle(`"`subtitle'"') ///
              `legend_opt' ///
              `scheme_opt' `saving_opt' `name_opt' `nodraw_opt'
    }
    else {
        // ===============================================================
        // HEAT MAP (bubble/scatter with color coding)
        // ===============================================================

        // Compute absolute estimate for marker sizing
        tempvar abs_est sign_est
        quietly {
            gen double `abs_est' = abs(`estimate_v')
            gen byte `sign_est' = cond(`estimate_v' >= 0, 1, 0)
        }

        // Separate positive and negative for color coding
        tempvar est_pos est_neg b1_pos b2_pos b1_neg b2_neg
        quietly {
            gen double `est_pos' = `abs_est' if `sign_est' == 1
            gen double `est_neg' = `abs_est' if `sign_est' == 0
            gen double `b1_pos' = `b1v' if `sign_est' == 1
            gen double `b2_pos' = `b2v' if `sign_est' == 1
            gen double `b1_neg' = `b1v' if `sign_est' == 0
            gen double `b2_neg' = `b2v' if `sign_est' == 0
        }

        // Positive effects: blue (#0072B2 = "0 114 178")
        // Negative effects: orange (#D55E00 = "213 94 0")
        local layers ""
        local layers `"`layers' (scatter `b2_pos' `b1_pos' [aweight=`est_pos'], mcolor("0 114 178%70") msymbol(O) msize(*2))"'
        local layers `"`layers' (scatter `b2_neg' `b1_neg' [aweight=`est_neg'], mcolor("213 94 0%70") msymbol(O) msize(*2))"'

        // Legend for heat map
        local legend_opt `"legend(order(1 "Positive effect" 2 "Negative effect") position(6) rows(1))"'

        twoway `layers' ///
            , xtitle("Boundary coordinate (x1)") ///
              ytitle("Boundary coordinate (x2)") ///
              title(`"`title'"') ///
              subtitle(`"`subtitle'"') ///
              `legend_opt' ///
              `scheme_opt' `saving_opt' `name_opt' `nodraw_opt'
    }

    restore
end

*! rd2d_plot version 1.1.0 24jun2026
*! Authors: Cattaneo, Titiunik, Yu
