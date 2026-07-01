*! version 1.1.0 24jun2026
program define rd2d_aggregate, rclass
    version 16.0
    syntax , [Method(string) Weights(numlist) SUBset(numlist integer >0) ///
             ESTimate(string) LEVEL(real 95)]

    // ─── Verify method() supplied ────────────────────────────────────────
    if "`method'" == "" {
        di as err "option method() required"
        di as err "  available methods: wbate, aate, lbate"
        exit 198
    }

    // ─── Verify estimation context ───────────────────────────────────────
    if "`e(cmd)'" != "rd2d" & "`e(cmd)'" != "rd2d_dist" {
        di as err "rd2d_aggregate requires estimation results from {bf:rd2d} or {bf:rd2d_dist}"
        exit 198
    }

    // ─── Method validation ───────────────────────────────────────────────
    local method = lower("`method'")
    if !inlist("`method'", "wbate", "aate", "lbate") {
        di as err "method() must be wbate, aate, or lbate; got '`method''"
        di as err "  wbate: weighted boundary average treatment effect"
        di as err "  aate:  aggregated average treatment effect (equal weights)"
        di as err "  lbate: largest boundary average treatment effect"
        exit 198
    }

    // ─── LBATE requires cbands ───────────────────────────────────────────
    if "`method'" == "lbate" {
        if "`e(cbands)'" != "on" {
            di as err "method(lbate) requires estimation with the {bf:cbands} option"
            di as err "re-run {bf:`e(cmd)'} with the {bf:cbands} option"
            exit 198
        }
    }

    // ─── WBATE/AATE require covariance matrix ────────────────────────────
    if "`method'" == "wbate" | "`method'" == "aate" {
        capture confirm matrix e(cov_q)
        if _rc {
            di as err "method(`method') requires the covariance matrix e(cov_q)"
            di as err "re-run {bf:`e(cmd)'} with the {bf:cbands} option"
            exit 198
        }
    }

    // ─── Estimate order validation ───────────────────────────────────────
    local est "q"
    if "`estimate'" != "" {
        local est = lower("`estimate'")
        if !inlist("`est'", "q", "p") {
            di as err "estimate() must be q or p; got '`estimate''"
            exit 198
        }
    }

    // ─── Level validation ────────────────────────────────────────────────
    if `level' <= 0 | `level' >= 100 {
        di as err "level() must be between 0 and 100 (exclusive); got `level'"
        exit 198
    }

    // ─── Subset processing ───────────────────────────────────────────────
    local neval = e(neval)
    local J = 0
    local indices ""

    if "`subset'" == "" {
        local J = `neval'
        forvalues i = 1/`neval' {
            local indices "`indices' `i'"
        }
    }
    else {
        foreach idx of numlist `subset' {
            if `idx' < 1 | `idx' > `neval' {
                di as err "subset() values must be between 1 and `neval'; got `idx'"
                exit 198
            }
            local J = `J' + 1
        }
        local indices "`subset'"
    }

    if `J' < 1 {
        di as err "no evaluation points selected"
        exit 198
    }

    // ─── WBATE weights validation ────────────────────────────────────────
    local wsum = 0
    if "`method'" == "wbate" {
        if "`weights'" == "" {
            di as err "method(wbate) requires the {bf:weights()} option"
            exit 198
        }
        local nw = 0
        foreach w of numlist `weights' {
            local nw = `nw' + 1
        }
        if `nw' != `J' {
            di as err "weights() must have `J' values (one per evaluation point); got `nw'"
            exit 198
        }
        foreach w of numlist `weights' {
            if `w' < 0 {
                di as err "weights() values must be non-negative; got `w'"
                exit 198
            }
        }
        foreach w of numlist `weights' {
            local wsum = `wsum' + `w'
        }
        if `wsum' <= 0 {
            di as err "weights() must have positive sum; got `wsum'"
            exit 198
        }
    }

    // ─── Column indices for e(results) ───────────────────────────────────
    // Layout: b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper ...
    if "`est'" == "q" {
        local col_est = 5
        local col_se  = 6
    }
    else {
        local col_est = 3
        local col_se  = 4
    }

    // ─── Extract results matrix ──────────────────────────────────────────
    tempname results
    mat `results' = e(results)

    // ═══════════════════════════════════════════════════════════════════════
    // WBATE / AATE computation
    // ═══════════════════════════════════════════════════════════════════════
    if "`method'" == "wbate" | "`method'" == "aate" {

        // Build weight vector (normalized to sum=1)
        tempname wvec
        mat `wvec' = J(`J', 1, .)
        if "`method'" == "aate" {
            forvalues i = 1/`J' {
                mat `wvec'[`i', 1] = 1 / `J'
            }
        }
        else {
            local wi = 0
            foreach w of numlist `weights' {
                local wi = `wi' + 1
                mat `wvec'[`wi', 1] = `w' / `wsum'
            }
        }

        // Extract estimation vector
        tempname est_vec
        mat `est_vec' = J(`J', 1, .)
        local ri = 0
        foreach ridx of numlist `indices' {
            local ri = `ri' + 1
            mat `est_vec'[`ri', 1] = `results'[`ridx', `col_est']
        }

        // Point estimate: tau = w' * est
        local tau_agg = 0
        forvalues i = 1/`J' {
            local tau_agg = `tau_agg' + `wvec'[`i', 1] * `est_vec'[`i', 1]
        }

        // Variance estimation (only for q-order with cov_q)
        local has_variance = 0
        if "`est'" == "q" {
            local has_variance = 1
            tempname covq covq_sub tmp_vec
            mat `covq' = e(cov_q)

            // Extract covariance submatrix
            mat `covq_sub' = J(`J', `J', .)
            local ri = 0
            foreach ridx of numlist `indices' {
                local ri = `ri' + 1
                local ci = 0
                foreach cidx of numlist `indices' {
                    local ci = `ci' + 1
                    mat `covq_sub'[`ri', `ci'] = `covq'[`ridx', `cidx']
                }
            }

            // Variance = w' * Cov * w
            mat `tmp_vec' = `covq_sub' * `wvec'
            local var_agg = 0
            forvalues i = 1/`J' {
                local var_agg = `var_agg' + `wvec'[`i', 1] * `tmp_vec'[`i', 1]
            }

            // Guard against numerical negatives
            if `var_agg' < 0 {
                local var_agg = 0
            }

            // Inference
            local se_agg = sqrt(`var_agg')
            if `se_agg' > 0 {
                local z_agg = `tau_agg' / `se_agg'
                local p_agg = 2 * (1 - normal(abs(`z_agg')))
            }
            else {
                local z_agg = .
                local p_agg = .
            }
            local z_alpha = invnormal((`level' + 100) / 200)
            local ci_lower = `tau_agg' - `z_alpha' * `se_agg'
            local ci_upper = `tau_agg' + `z_alpha' * `se_agg'
        }

        // ─── Display ─────────────────────────────────────────────────────
        di ""
        if "`method'" == "wbate" {
            di as txt "Weighted Boundary Average Treatment Effect (WBATE)"
        }
        else {
            di as txt "Aggregated Average Treatment Effect (AATE)"
        }
        di as txt "{hline 68}"
        di as txt "  Method:       " as res "`=upper("`method'")'"
        di as txt "  Eval. points: " as res "`J'"
        if "`method'" == "aate" {
            di as txt "  Weights:      " as res "equal (1/`J')"
        }
        else {
            di as txt "  Weights:      " as res "user-specified (normalized)"
        }
        di as txt "  Estimate:     " as res "`=upper("`est'")'-order"

        if `has_variance' {
            di as txt "{hline 68}"
            di as txt _col(3) "Estimate" _col(15) "Std. Err." _col(28) ///
                "z" _col(35) "P>|z|" _col(50) "[`level'% CI]"
            di as txt "{hline 68}"
            di as res _col(3) %10.6g `tau_agg' _col(15) %10.6g `se_agg' ///
                _col(27) %8.4f `z_agg' _col(36) %6.4f `p_agg' _col(47) ///
                "[" %9.6g `ci_lower' ", " %9.6g `ci_upper' "]"
            di as txt "{hline 68}"
        }
        else {
            di as txt "{hline 68}"
            di as txt _col(3) "Estimate"
            di as txt "{hline 68}"
            di as res _col(3) %10.6g `tau_agg'
            di as txt "{hline 68}"
            di as txt _col(3) "Note: Variance not available for P-order estimates"
            di as txt _col(3) "(covariance matrix e(cov_q) is Q-order only)."
        }

        // ─── Return values ───────────────────────────────────────────────
        return scalar estimate = `tau_agg'
        return scalar J = `J'
        return scalar level = `level'
        return local method "`method'"
        return local estimate_order "`est'"
        if `has_variance' {
            return scalar se = `se_agg'
            return scalar z = `z_agg'
            return scalar p = `p_agg'
            return scalar ci_lower = `ci_lower'
            return scalar ci_upper = `ci_upper'
            return scalar variance = `var_agg'
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // LBATE computation
    // ═══════════════════════════════════════════════════════════════════════
    if "`method'" == "lbate" {

        // LBATE CI defined only via q-order uniform band (Theorem 6)
        if "`est'" != "q" {
            di as err "method(lbate) requires q-order estimates;"
            di as err "LBATE confidence intervals are defined only via the"
            di as err "uniform band critical value from q-order inference."
            di as err "Remove the estimate(p) option or re-estimate with cbands."
            exit 198
        }

        // Find maximum estimate
        local tau_lbate = -c(maxdouble)
        local jmax = 1
        local ji = 0
        foreach idx of numlist `indices' {
            local ji = `ji' + 1
            local tau_j = `results'[`idx', `col_est']
            if `tau_j' > `tau_lbate' {
                local tau_lbate = `tau_j'
                local jmax = `idx'
            }
        }

        // CI using uniform band critical value (Theorem 6)
        local cb_crit = e(cb_crit)
        local ci_lower = -c(maxdouble)
        local ci_upper = -c(maxdouble)
        foreach idx of numlist `indices' {
            local tau_j = `results'[`idx', `col_est']
            local se_j = `results'[`idx', `col_se']
            local lb_j = `tau_j' - `cb_crit' * `se_j'
            local ub_j = `tau_j' + `cb_crit' * `se_j'
            if `lb_j' > `ci_lower' {
                local ci_lower = `lb_j'
            }
            if `ub_j' > `ci_upper' {
                local ci_upper = `ub_j'
            }
        }

        // ─── Display ─────────────────────────────────────────────────────
        di ""
        di as txt "Largest Boundary Average Treatment Effect (LBATE)"
        di as txt "{hline 68}"
        di as txt "  Method:       " as res "LBATE"
        di as txt "  Eval. points: " as res "`J'"
        di as txt "  Max at:       " as res "point `jmax'"
        di as txt "  CB critical:  " as res %9.4g `cb_crit'
        di as txt "  Estimate:     " as res "`=upper("`est'")'-order"
        di as txt "{hline 68}"
        di as txt _col(3) "Estimate" _col(50) "[`level'% CI]"
        di as txt "{hline 68}"
        di as res _col(3) %10.6g `tau_lbate' _col(47) ///
            "[" %9.6g `ci_lower' ", " %9.6g `ci_upper' "]"
        di as txt "{hline 68}"
        di as txt _col(3) "Note: CI based on uniform band critical value (Theorem 6,"
        di as txt _col(3) "Cattaneo-Titiunik-Yu, 2025). No z-statistic for LBATE."

        // ─── Return values ───────────────────────────────────────────────
        return scalar estimate = `tau_lbate'
        return scalar ci_lower = `ci_lower'
        return scalar ci_upper = `ci_upper'
        return scalar cb_crit = `cb_crit'
        return scalar jmax = `jmax'
        return scalar J = `J'
        return scalar level = `level'
        return local method "lbate"
        return local estimate_order "`est'"
    }

end
