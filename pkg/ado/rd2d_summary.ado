*! version 1.1.0 24jun2026
* rd2d_summary: Post-estimation summary for rd2d and rd2d_dist
* Provides structured, publication-quality display of estimation results.

program define rd2d_summary
    version 16.0
    syntax [, OUTput(string) CBuniform SUBset(numlist integer >0) ///
             AATE(numlist)]

    // =========================================================================
    // 1. Validation: require rd2d or rd2d_dist estimation results
    // =========================================================================
    if "`e(cmd)'" != "rd2d" & "`e(cmd)'" != "rd2d_dist" {
        di as err "rd2d_summary requires estimation results from rd2d or rd2d_dist"
        exit 301
    }

    local cmd "`e(cmd)'"
    local neval = e(neval)

    // output() option: default "main"
    if "`output'" == "" local output "main"
    local output = lower("`output'")
    if "`output'" != "main" & "`output'" != "bw" {
        di as err "output() must be {bf:main} or {bf:bw}; got '`output''"
        exit 198
    }

    // cbuniform: requires cbands estimation
    if "`cbuniform'" != "" & "`e(cbands)'" != "on" {
        di as err "cbuniform requires estimation with cbands option"
        exit 198
    }

    // subset: validate range
    local indices ""
    if "`subset'" == "" {
        forvalues i = 1/`neval' {
            local indices "`indices' `i'"
        }
    }
    else {
        foreach idx of numlist `subset' {
            if `idx' < 1 | `idx' > `neval' {
                di as err "subset() values must be between 1 and `neval'"
                exit 198
            }
        }
        local indices "`subset'"
    }
    local ndisp : word count `indices'

    // aate: validate length matches displayed points
    if "`aate'" != "" {
        local nwt : word count `aate'
        if `nwt' != `ndisp' {
            di as err "aate() must have `ndisp' weights (one per displayed evaluation point); got `nwt'"
            exit 198
        }
    }

    // =========================================================================
    // 2. Extract stored results
    // =========================================================================
    tempname results bws
    matrix `results' = e(results)
    matrix `bws' = e(bws)
    local level = e(level)

    // Determine if b1/b2 columns are all missing (rd2d_dist case)
    local has_bdy 1
    if "`cmd'" == "rd2d_dist" {
        local has_bdy 0
        forvalues i = 1/`neval' {
            if `results'[`i', 1] < . | `results'[`i', 2] < . {
                local has_bdy 1
                continue, break
            }
        }
    }

    // =========================================================================
    // 3. Display header (meta-information)
    // =========================================================================
    di ""
    di as txt "`cmd'"
    di ""

    // Total observations
    di as txt "Number of Obs." _col(40) as res %10.0f e(N)

    // BW type
    if "`e(bwsource)'" == "user" {
        di as txt "BW type" _col(40) as res %10s "user"
    }
    else {
        local bwdisp "`e(bwselect)'"
        if "`cmd'" == "rd2d" {
            local bwdisp "`e(bwselect)'-`e(method)'"
        }
        di as txt "BW type" _col(40) as res %10s "`bwdisp'"
    }

    // Kernel
    if "`cmd'" == "rd2d" {
        local kerndisp "`e(kernel)'-`e(ktype)'"
    }
    else {
        local kerndisp "`e(kernel)'"
    }
    di as txt "Kernel" _col(40) as res %10s "`kerndisp'"

    // VCE
    di as txt "VCE method" _col(40) as res %10s "`e(vce)'"

    // Masspoints
    di as txt "Masspoints" _col(40) as res %10s "`e(masspoints_opt)'"

    // Kink (rd2d_dist only)
    if "`cmd'" == "rd2d_dist" {
        di as txt "Kink" _col(40) as res %10s "`e(kink)'"
    }

    di ""

    // Control/Treated breakdown
    di as txt _col(35) "Control" _col(48) "Treated"
    di as txt "Number of Obs." _col(35) as res %7.0f e(N0) _col(48) as res %7.0f e(N1)
    di as txt "Order est. (p)" _col(35) as res %7.0f e(p) _col(48) as res %7.0f e(p)
    di as txt "Order rbc. (q)" _col(35) as res %7.0f e(q) _col(48) as res %7.0f e(q)

    di ""

    // =========================================================================
    // 4. Main table output
    // =========================================================================
    if "`output'" == "main" {
        // --- Determine CI/CB label ---
        local ci_label "`=string(`level', "%3.0f")'% CI"
        if "`cbuniform'" != "" {
            local ci_label "`=string(`level', "%3.0f")'% Unif. CB"
        }

        // --- Table header ---
        if `has_bdy' {
            local hdr_line "===================================================================="
            di as txt "`hdr_line'"
            di as txt "  ID" _col(9) "     b1" _col(18) "     b2" _col(27) "   Est." _col(36) "      z" _col(45) "  P>|z|" _col(55) "    `ci_label'"
            di as txt "`hdr_line'"
        }
        else {
            local hdr_line "========================================================"
            di as txt "`hdr_line'"
            di as txt "  ID" _col(9) "   Est." _col(18) "      z" _col(27) "  P>|z|" _col(37) "    `ci_label'"
            di as txt "`hdr_line'"
        }

        // --- Table body ---
        foreach idx of numlist `indices' {
            local est  = `results'[`idx', 5]
            local zval = `results'[`idx', 7]
            local pval = `results'[`idx', 8]

            // CI or CB
            if "`cbuniform'" != "" {
                local ci_lo = `results'[`idx', 11]
                local ci_hi = `results'[`idx', 12]
            }
            else {
                local ci_lo = `results'[`idx', 9]
                local ci_hi = `results'[`idx', 10]
            }

            // Format interval bounds (handle inf)
            local lo_str = string(`ci_lo', "%7.4f")
            local hi_str = string(`ci_hi', "%7.4f")
            if `ci_lo' <= -c(maxdouble) + 1 | `ci_lo' >= . {
                local lo_str "   -inf"
            }
            if `ci_hi' >= c(maxdouble) - 1 | `ci_hi' >= . {
                local hi_str "    inf"
            }

            local interval "[`lo_str', `hi_str']"

            if `has_bdy' {
                local b1val = `results'[`idx', 1]
                local b2val = `results'[`idx', 2]
                di as res " " %3.0f `idx' " " %8.3f `b1val' " " %8.3f `b2val' " " %7.4f `est' " " %8.4f `zval' " " %8.4f `pval' "   `interval'"
            }
            else {
                di as res " " %3.0f `idx' " " %7.4f `est' " " %8.4f `zval' " " %8.4f `pval' "   `interval'"
            }
        }

        // --- Table footer ---
        if `has_bdy' {
            di as txt "===================================================================="
        }
        else {
            di as txt "========================================================"
        }

        // =================================================================
        // 5. AATE aggregation (if requested)
        // =================================================================
        if "`aate'" != "" {
            // Check cov_q availability
            capture confirm matrix e(cov_q)
            if _rc {
                di as err "aate() requires estimation with cbands option (e(cov_q) not found)"
                exit 198
            }

            tempname covq covq_sub wvec tau_vec
            matrix `covq' = e(cov_q)

            // Build weight vector and estimates vector
            matrix `wvec' = J(1, `ndisp', 0)
            matrix `tau_vec' = J(`ndisp', 1, 0)
            local wsum = 0
            local wi = 1
            foreach w of numlist `aate' {
                matrix `wvec'[1, `wi'] = `w'
                local wsum = `wsum' + `w'
                local wi = `wi' + 1
            }
            // Normalize weights
            if `wsum' <= 0 {
                di as err "aate() weights must sum to a positive value"
                exit 198
            }
            forvalues i = 1/`ndisp' {
                matrix `wvec'[1, `i'] = `wvec'[1, `i'] / `wsum'
            }

            // Extract subset estimate values and covariance sub-matrix
            matrix `covq_sub' = J(`ndisp', `ndisp', 0)
            local ri = 1
            foreach ridx of numlist `indices' {
                matrix `tau_vec'[`ri', 1] = `results'[`ridx', 5]
                local ci = 1
                foreach cidx of numlist `indices' {
                    matrix `covq_sub'[`ri', `ci'] = `covq'[`ridx', `cidx']
                    local ci = `ci' + 1
                }
                local ri = `ri' + 1
            }

            // AATE = w' * tau
            local aate_val = 0
            forvalues i = 1/`ndisp' {
                local aate_val = `aate_val' + `wvec'[1, `i'] * `tau_vec'[`i', 1]
            }

            // Var(AATE) = w' * Cov * w
            tempname wcol wcov_prod
            matrix `wcol' = `wvec''
            matrix `wcov_prod' = `wvec' * `covq_sub' * `wcol'
            local aate_var = `wcov_prod'[1, 1]
            if `aate_var' < 0 {
                local aate_var = 0
            }
            local aate_se = sqrt(`aate_var')

            // z, p, CI
            if `aate_se' > 0 {
                local aate_z = `aate_val' / `aate_se'
            }
            else {
                local aate_z = .
            }
            if `aate_z' < . {
                local aate_p = 2 * (1 - normal(abs(`aate_z')))
            }
            else {
                local aate_p = .
            }
            local crit = invnormal((`level' + 100) / 200)
            local aate_ci_lo = `aate_val' - `crit' * `aate_se'
            local aate_ci_hi = `aate_val' + `crit' * `aate_se'

            // Display AATE
            di as txt "--------------------------------------------------------------------"
            di as txt "AATE (Aggregated Average Treatment Effect)"
            di as txt "--------------------------------------------------------------------"

            // Display weights
            di as txt "  Weights:" _c
            forvalues i = 1/`ndisp' {
                di as res "  " %6.4f `wvec'[1, `i'] _c
            }
            di ""

            di as txt "  AATE:" _col(20) as res %12.4f `aate_val'
            di as txt "  Std. Err.:" _col(20) as res %12.4f `aate_se'
            di as txt "  z:" _col(20) as res %12.4f `aate_z'
            di as txt "  P>|z|:" _col(20) as res %12.4f `aate_p'
            di as txt "  `=string(`level', "%3.0f")'% CI:" _col(20) as res "  [" %7.4f `aate_ci_lo' ", " %7.4f `aate_ci_hi' "]"
            di as txt "--------------------------------------------------------------------"
        }
    }

    // =========================================================================
    // 6. Bandwidth table output
    // =========================================================================
    if "`output'" == "bw" {
        if "`cmd'" == "rd2d" {
            // rd2d bws: b1 b2 h01 h02 h11 h12 Nh0 Nh1
            local bw_line "==================================================================================="
            di as txt "`bw_line'"
            di as txt _col(15) "Bdy Points" _col(35) "BW Control" _col(55) "BW Treatment" _col(72) "Eff. N"
            di as txt "  ID" _col(9) "     b1" _col(18) "     b2" _col(27) "    h01" _col(36) "    h02" _col(45) "    h11" _col(54) "    h12" _col(63) "    Nh0" _col(72) "    Nh1"
            di as txt "`bw_line'"

            foreach idx of numlist `indices' {
                local b1v = `bws'[`idx', 1]
                local b2v = `bws'[`idx', 2]
                local h01 = `bws'[`idx', 3]
                local h02 = `bws'[`idx', 4]
                local h11 = `bws'[`idx', 5]
                local h12 = `bws'[`idx', 6]
                local nh0 = `bws'[`idx', 7]
                local nh1 = `bws'[`idx', 8]
                di as res " " %3.0f `idx' " " %8.3f `b1v' " " %8.3f `b2v' " " %8.3f `h01' " " %8.3f `h02' " " %8.3f `h11' " " %8.3f `h12' " " %7.0f `nh0' " " %7.0f `nh1'
            }

            di as txt "`bw_line'"
        }
        else {
            // rd2d_dist bws: b1 b2 h0 h1 Nh0 Nh1
            if `has_bdy' {
                local bw_line "================================================================="
                di as txt "`bw_line'"
                di as txt "  ID" _col(9) "     b1" _col(18) "     b2" _col(27) "     h0" _col(36) "     h1" _col(45) "    Nh0" _col(54) "    Nh1"
                di as txt "`bw_line'"

                foreach idx of numlist `indices' {
                    local b1v = `bws'[`idx', 1]
                    local b2v = `bws'[`idx', 2]
                    local h0  = `bws'[`idx', 3]
                    local h1  = `bws'[`idx', 4]
                    local nh0 = `bws'[`idx', 5]
                    local nh1 = `bws'[`idx', 6]
                    di as res " " %3.0f `idx' " " %8.3f `b1v' " " %8.3f `b2v' " " %8.3f `h0' " " %8.3f `h1' " " %7.0f `nh0' " " %7.0f `nh1'
                }

                di as txt "`bw_line'"
            }
            else {
                local bw_line "================================================="
                di as txt "`bw_line'"
                di as txt "  ID" _col(9) "     h0" _col(18) "     h1" _col(27) "    Nh0" _col(36) "    Nh1"
                di as txt "`bw_line'"

                foreach idx of numlist `indices' {
                    local h0  = `bws'[`idx', 3]
                    local h1  = `bws'[`idx', 4]
                    local nh0 = `bws'[`idx', 5]
                    local nh1 = `bws'[`idx', 6]
                    di as res " " %3.0f `idx' " " %8.3f `h0' " " %8.3f `h1' " " %7.0f `nh0' " " %7.0f `nh1'
                }

                di as txt "`bw_line'"
            }
        }
    }
end
