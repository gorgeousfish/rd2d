*! version 1.1.0 30jun2026
* rd2d_bwsens: Bandwidth sensitivity analysis for rd2d and rd2d_dist
* Re-estimates the model over a grid of bandwidth multipliers to assess
* robustness of point estimates, standard errors, and confidence intervals
* to the bandwidth choice.

program define rd2d_bwsens, rclass
    version 16.0
    syntax [, GRID(numlist >0) SUBset(numlist integer >0) PLOT ///
             SAVing(string) QUIET FORmat(string)]

    // =========================================================================
    // 1. Validation: require rd2d or rd2d_dist estimation results
    // =========================================================================
    if "`e(cmd)'" != "rd2d" & "`e(cmd)'" != "rd2d_dist" {
        di as err "rd2d_bwsens requires estimation results from {bf:rd2d} or {bf:rd2d_dist}"
        exit 301
    }

    local cmd "`e(cmd)'"
    local neval = e(neval)

    // Check for stored variable names (requires updated rd2d/rd2d_dist)
    if "`e(depvar)'" == "" {
        di as err "rd2d_bwsens requires e(depvar); please re-run `cmd' with the latest version"
        exit 301
    }

    // =========================================================================
    // 2. Parse options
    // =========================================================================
    // Default grid of bandwidth multipliers
    if "`grid'" == "" {
        local grid "0.5 0.75 0.9 1.0 1.1 1.25 1.5 2.0"
    }
    local ngrid : word count `grid'

    // Display format
    if "`format'" == "" local format "%9.4f"

    // Subset validation
    local eval_indices ""
    if "`subset'" == "" {
        forvalues i = 1/`neval' {
            local eval_indices "`eval_indices' `i'"
        }
    }
    else {
        foreach idx of numlist `subset' {
            if `idx' < 1 | `idx' > `neval' {
                di as err "subset() values must be between 1 and `neval'"
                exit 198
            }
        }
        local eval_indices "`subset'"
    }
    local ndisp : word count `eval_indices'

    // =========================================================================
    // 3. Extract current estimation settings
    // =========================================================================
    local depvar   "`e(depvar)'"
    local kernel   "`e(kernel)'"
    local vce      "`e(vce)'"
    local rbc      "`e(rbc)'"
    local side     "`e(side)'"
    local level    = e(level)
    local p_order  = e(p)
    local q_order  = e(q)
    local cluster  "`e(cluster)'"
    local masspts  "`e(masspoints_opt)'"
    local scaleregul = e(scaleregul)

    // Extract baseline bandwidth matrix
    tempname bws_base
    matrix `bws_base' = e(bws)

    // Command-specific settings
    if "`cmd'" == "rd2d" {
        local x1var    "`e(x1var)'"
        local x2var    "`e(x2var)'"
        local tvar     "`e(tvar)'"
        local atstring `"`e(atstring)'"'
        local ktype    "`e(ktype)'"
        local bwselect "`e(bwselect)'"
        local method   "`e(method)'"
        local deriv    "`e(deriv)'"
        local tangvec  "`e(tangvec)'"
        local stdvars  "`e(stdvars)'"
    }
    else {
        // rd2d_dist
        local dvars    "`e(dvars)'"
        local bwselect "`e(bwselect)'"
        local kink     "`e(kink)'"
    }

    // =========================================================================
    // 4. Build baseline bandwidth values for scaling
    // =========================================================================
    // For rd2d: bws has columns b1 b2 h01 h02 h11 h12 Nh0 Nh1
    // For rd2d_dist: bws has columns b1 b2 h0 h1 Nh0 Nh1
    if "`cmd'" == "rd2d" {
        // Extract baseline h values: 4 per eval point (h01 h02 h11 h12)
        local base_h ""
        forvalues j = 1/`neval' {
            local h01 = `bws_base'[`j', 3]
            local h02 = `bws_base'[`j', 4]
            local h11 = `bws_base'[`j', 5]
            local h12 = `bws_base'[`j', 6]
            local base_h "`base_h' `h01' `h02' `h11' `h12'"
        }
    }
    else {
        // Extract baseline h values: 2 per eval point (h0 h1)
        local base_h ""
        forvalues j = 1/`neval' {
            local h0 = `bws_base'[`j', 3]
            local h1 = `bws_base'[`j', 4]
            local base_h "`base_h' `h0' `h1'"
        }
    }

    // =========================================================================
    // 5. Reconstruct common options string
    // =========================================================================
    local common_opts "p(`p_order') q(`q_order') kernel(`kernel') vce(`vce')"
    local common_opts "`common_opts' level(`level') side(`side') rbc(`rbc')"
    local common_opts "`common_opts' masspoints(`masspts')"
    if "`cluster'" != "" {
        local common_opts "`common_opts' cluster(`cluster')"
    }

    if "`cmd'" == "rd2d" {
        local common_opts "`common_opts' ktype(`ktype') scaleregul(`scaleregul')"
        if "`deriv'" != "" & "`deriv'" != "0 0" {
            local common_opts "`common_opts' deriv(`deriv')"
        }
        if "`tangvec'" != "" {
            local common_opts "`common_opts' tangvec(`tangvec')"
        }
    }
    else {
        local common_opts "`common_opts' scaleregul(`scaleregul')"
        if "`kink'" == "on" {
            local common_opts "`common_opts' kink(on)"
        }
    }

    // =========================================================================
    // 6. Iterate over grid and re-estimate
    // =========================================================================
    // Result storage: ngrid rows x (ndisp * 4 + 1) columns
    // For each grid point: [multiplier, {Est_q, Se_q, CI_lo, CI_hi} x ndisp]
    tempname sens_results
    matrix `sens_results' = J(`ngrid', 1 + 4 * `ndisp', .)

    local gi = 0
    foreach mult of numlist `grid' {
        local gi = `gi' + 1
        matrix `sens_results'[`gi', 1] = `mult'

        // Build scaled bandwidth string
        local h_scaled ""
        if "`cmd'" == "rd2d" {
            local nhvals = 4 * `neval'
            forvalues i = 1/`nhvals' {
                local bh : word `i' of `base_h'
                local sv = `mult' * `bh'
                local h_scaled "`h_scaled' `sv'"
            }
        }
        else {
            local nhvals = 2 * `neval'
            forvalues i = 1/`nhvals' {
                local bh : word `i' of `base_h'
                local sv = `mult' * `bh'
                local h_scaled "`h_scaled' `sv'"
            }
        }

        // Reconstruct and run the command quietly
        if "`quiet'" == "" & "`mult'" == "1" {
            // Baseline multiplier - note for user
        }

        capture {
            if "`cmd'" == "rd2d" {
                quietly rd2d `depvar' `x1var' `x2var' `tvar', ///
                    at(`atstring') h(`h_scaled') `common_opts'
            }
            else {
                quietly rd2d_dist `depvar' `dvars', ///
                    h(`h_scaled') `common_opts'
            }
        }

        if _rc {
            if "`quiet'" == "" {
                di as txt "note: estimation failed at multiplier `mult' (rc=" _rc ")"
            }
            continue
        }

        // Extract results for displayed evaluation points
        tempname tmp_results
        matrix `tmp_results' = e(results)
        local ci = 2
        foreach idx of numlist `eval_indices' {
            matrix `sens_results'[`gi', `ci']     = `tmp_results'[`idx', 5]
            matrix `sens_results'[`gi', `ci' + 1] = `tmp_results'[`idx', 6]
            matrix `sens_results'[`gi', `ci' + 2] = `tmp_results'[`idx', 9]
            matrix `sens_results'[`gi', `ci' + 3] = `tmp_results'[`idx', 10]
            local ci = `ci' + 4
        }
    }

    // =========================================================================
    // 7. Restore original estimation (re-run at multiplier 1.0)
    // =========================================================================
    capture {
        if "`cmd'" == "rd2d" {
            quietly rd2d `depvar' `x1var' `x2var' `tvar', ///
                at(`atstring') h(`base_h') `common_opts'
        }
        else {
            quietly rd2d_dist `depvar' `dvars', ///
                h(`base_h') `common_opts'
        }
    }

    // =========================================================================
    // 8. Display formatted table
    // =========================================================================
    if "`quiet'" == "" {
        di ""
        di as txt "Bandwidth Sensitivity Analysis"
        di as txt "Command: {bf:`cmd'}"
        di as txt "Baseline bandwidth source: `e(bwsource)'"
        di as txt "Grid multipliers: `grid'"
        di ""

        // Display one sub-table per evaluation point
        local di_idx = 0
        foreach idx of numlist `eval_indices' {
            local di_idx = `di_idx' + 1

            // Header for this evaluation point
            if `ndisp' > 1 | "`subset'" != "" {
                if "`cmd'" == "rd2d" {
                    local b1v = `bws_base'[`idx', 1]
                    local b2v = `bws_base'[`idx', 2]
                    di as txt "--- Evaluation point `idx' (b1=" ///
                        string(`b1v', "%7.3f") ", b2=" string(`b2v', "%7.3f") ") ---"
                }
                else {
                    di as txt "--- Evaluation point `idx' ---"
                }
            }

            local hdr_line "=================================================================="
            di as txt "`hdr_line'"
            di as txt "  Mult." _col(12) "    h" _col(24) "   Est.q" _col(36) "   Se.q" _col(48) "  `=string(`level', "%3.0f")'% CI"
            di as txt "`hdr_line'"

            forvalues gi = 1/`ngrid' {
                local mv = `sens_results'[`gi', 1]
                local col_start = 2 + 4 * (`di_idx' - 1)
                local est_v  = `sens_results'[`gi', `col_start']
                local se_v   = `sens_results'[`gi', `col_start' + 1]
                local ci_lo  = `sens_results'[`gi', `col_start' + 2]
                local ci_hi  = `sens_results'[`gi', `col_start' + 3]

                // Compute displayed bandwidth (average of dimensions)
                if "`cmd'" == "rd2d" {
                    local h01 = `mv' * `bws_base'[`idx', 3]
                    local h02 = `mv' * `bws_base'[`idx', 4]
                    local h_disp = (`h01' + `h02') / 2
                }
                else {
                    local h0 = `mv' * `bws_base'[`idx', 3]
                    local h_disp = `h0'
                }

                if `est_v' >= . {
                    di as res " " %5.2f `mv' _col(12) `format' `h_disp' _col(24) "       ." _col(36) "       ." _col(48) "  [., .]"
                }
                else {
                    local lo_str = string(`ci_lo', "`format'")
                    local hi_str = string(`ci_hi', "`format'")
                    // Mark baseline
                    local star ""
                    if abs(`mv' - 1.0) < 1e-8 {
                        local star " *"
                    }
                    di as res " " %5.2f `mv' _col(12) `format' `h_disp' _col(24) `format' `est_v' _col(36) `format' `se_v' _col(48) "  [`lo_str', `hi_str']`star'"
                }
            }

            di as txt "`hdr_line'"
            di as txt "  * denotes baseline bandwidth (multiplier = 1.0)"
            di ""
        }
    }

    // =========================================================================
    // 9. Optional: sensitivity plot
    // =========================================================================
    if "`plot'" != "" {
        // Create temporary dataset for plotting
        preserve
        quietly {
            drop _all
            set obs `ngrid'

            // Generate multiplier variable
            generate double _mult = .
            forvalues gi = 1/`ngrid' {
                replace _mult = `sens_results'[`gi', 1] in `gi'
            }

            // Generate estimate and CI for first displayed point
            generate double _est = .
            generate double _ci_lo = .
            generate double _ci_hi = .
            forvalues gi = 1/`ngrid' {
                replace _est   = `sens_results'[`gi', 2] in `gi'
                replace _ci_lo = `sens_results'[`gi', 4] in `gi'
                replace _ci_hi = `sens_results'[`gi', 5] in `gi'
            }
        }

        local plot_title "Bandwidth Sensitivity: `cmd'"
        local plot_opts "xtitle(Bandwidth multiplier) ytitle(Treatment effect)"
        local plot_opts "`plot_opts' title(`plot_title')"
        local plot_opts "`plot_opts' legend(order(1 \"Point estimate\" 2 \"`=string(`level', "%3.0f")'% CI\"))"
        local plot_opts "`plot_opts' xline(1, lpattern(dash) lcolor(gs8))"

        twoway (connected _est _mult, msymbol(O) mcolor(navy) lcolor(navy)) ///
               (rcap _ci_hi _ci_lo _mult, lcolor(cranberry)) ///
               , `plot_opts'

        if "`saving'" != "" {
            graph export "`saving'", replace
        }

        restore
    }

    // =========================================================================
    // 10. Return results
    // =========================================================================
    // Build column names
    local cnames "multiplier"
    foreach idx of numlist `eval_indices' {
        local cnames "`cnames' Est_`idx' Se_`idx' CI_lo_`idx' CI_hi_`idx'"
    }
    matrix colnames `sens_results' = `cnames'

    // Row names from multipliers
    local rnames ""
    foreach mult of numlist `grid' {
        local rnames "`rnames' m`=subinstr("`mult'", ".", "_", .)'"
    }
    matrix rownames `sens_results' = `rnames'

    return matrix sens_results = `sens_results'
    return local grid "`grid'"
    return scalar ngrid = `ngrid'
    return scalar ndisp = `ndisp'
    return local cmd "`cmd'"
end
