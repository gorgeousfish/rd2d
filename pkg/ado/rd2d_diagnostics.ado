*! version 1.1.0 30jun2026
* rd2d_diagnostics: Post-estimation diagnostics for rd2d and rd2d_dist
* Formats and displays estimation quality diagnostics from e(diagnostics).

program define rd2d_diagnostics, rclass
    version 16.0
    syntax [, OUTput(string) SUBset(numlist integer >0)]

    // =========================================================================
    // 1. Validation: require rd2d or rd2d_dist estimation results
    // =========================================================================
    if "`e(cmd)'" != "rd2d" & "`e(cmd)'" != "rd2d_dist" {
        di as err "rd2d_diagnostics requires estimation results from rd2d or rd2d_dist"
        exit 301
    }

    local cmd "`e(cmd)'"
    local neval = e(neval)

    // output() option: default "summary"
    if "`output'" == "" local output "summary"
    local output = lower("`output'")
    if "`output'" != "summary" & "`output'" != "full" & "`output'" != "warnings" {
        di as err "output() must be {bf:summary}, {bf:full}, or {bf:warnings}; got '`output''"
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

    // =========================================================================
    // 2. Extract stored results
    // =========================================================================
    tempname diag results massinfo
    matrix `diag' = e(diagnostics)
    matrix `results' = e(results)
    capture matrix `massinfo' = e(masspoints)
    local has_mass = (_rc == 0)

    local p = e(p)
    local q = e(q)
    // Full-rank basis count for 2D polynomial of order k: (k+1)(k+2)/2
    local fullrank_p = (`p' + 1) * (`p' + 2) / 2
    local fullrank_q = (`q' + 1) * (`q' + 2) / 2

    // CB diagnostics (may not exist)
    local cb_psd = 0
    local cb_mineig = .
    local cb_crit = .
    capture {
        local cb_psd = e(cb_psd_adjusted)
        local cb_mineig = e(cb_min_eigen)
        local cb_crit = e(cb_crit)
    }

    // =========================================================================
    // 3. Compute aggregate diagnostics across all requested points
    // =========================================================================
    local max_cond = 0
    local n_fallback = 0
    local min_effN = .
    local has_warnings = 0

    foreach idx of numlist `indices' {
        // Condition numbers: columns 4, 7, 10, 13
        forvalues col = 4(3)13 {
            local cval = `diag'[`idx', `col']
            if `cval' < . & `cval' > `max_cond' {
                local max_cond = `cval'
            }
        }
        // Fallback flags: columns 5, 8, 11, 14
        forvalues col = 5(3)14 {
            local fval = `diag'[`idx', `col']
            if `fval' == 1 {
                local n_fallback = `n_fallback' + 1
                continue, break
            }
        }
        // Effective sample sizes: results columns 17, 18
        local nh0 = `results'[`idx', 17]
        local nh1 = `results'[`idx', 18]
        if `nh0' < . & `nh0' < `min_effN' {
            local min_effN = `nh0'
        }
        if `nh1' < . & `nh1' < `min_effN' {
            local min_effN = `nh1'
        }
    }

    // Determine if any warnings exist
    if `max_cond' > 1e6 | `n_fallback' > 0 {
        local has_warnings = 1
    }
    if `min_effN' < 30 & `min_effN' < . {
        local has_warnings = 1
    }
    if `has_mass' {
        foreach idx of numlist `indices' {
            local mratio = `massinfo'[`idx', 4]
            if `mratio' > 0.2 & `mratio' < . {
                local has_warnings = 1
            }
        }
    }

    // =========================================================================
    // 4. Display header
    // =========================================================================
    di ""
    di as txt "{bf:rd2d diagnostics}" _col(40) "after " as res "`cmd'"
    di ""
    di as txt "Polynomial order (p)" _col(40) as res %5.0f `p'
    di as txt "Bias-correction order (q)" _col(40) as res %5.0f `q'
    di as txt "Full-rank basis (p)" _col(40) as res %5.0f `fullrank_p'
    di as txt "Full-rank basis (q)" _col(40) as res %5.0f `fullrank_q'
    di as txt "Number of Obs." _col(40) as res %10.0f e(N)
    di as txt "  Control" _col(40) as res %10.0f e(N0)
    di as txt "  Treated" _col(40) as res %10.0f e(N1)
    di ""

    // CB PSD info if available
    if `cb_crit' < . {
        di as txt "CB critical value" _col(40) as res %10.4f `cb_crit'
        if `cb_psd' == 1 {
            di as txt "CB PSD adjustment" _col(40) as res "{bf:applied}"
            di as txt "CB min eigenvalue" _col(40) as res %10.2e `cb_mineig'
        }
        else {
            di as txt "CB PSD adjustment" _col(40) as res "none"
        }
        di ""
    }

    // =========================================================================
    // 5. Summary mode
    // =========================================================================
    if "`output'" == "summary" {
        local hdr_line "========================================================================================"
        di as txt "`hdr_line'"
        di as txt "  ID" _col(8) "  Nh0" _col(15) "  Nh1" _col(22) " cond_p0" _col(32) " cond_p1" _col(42) " cond_q0" _col(52) " cond_q1" _col(62) "fallback" _col(72) " mass_r" _col(80) "status"
        di as txt "`hdr_line'"

        foreach idx of numlist `indices' {
            local nh0 = `results'[`idx', 17]
            local nh1 = `results'[`idx', 18]
            local cp0 = `diag'[`idx', 4]
            local cp1 = `diag'[`idx', 7]
            local cq0 = `diag'[`idx', 10]
            local cq1 = `diag'[`idx', 13]

            // Determine fallback status
            local fb_str "invsym"
            local fp0 = `diag'[`idx', 5]
            local fp1 = `diag'[`idx', 8]
            local fq0 = `diag'[`idx', 11]
            local fq1 = `diag'[`idx', 14]
            if `fp0' == 1 | `fp1' == 1 | `fq0' == 1 | `fq1' == 1 {
                local fb_str "pinv"
            }

            // Mass ratio
            local mr_str "    ."
            if `has_mass' {
                local mratio = `massinfo'[`idx', 4]
                if `mratio' < . {
                    local mr_str = string(`mratio', "%6.3f")
                }
            }

            // Status flag
            local status_str "ok"
            if `cp0' > 1e6 | `cp1' > 1e6 | `cq0' > 1e6 | `cq1' > 1e6 {
                local status_str "WARN"
            }
            if "`fb_str'" == "pinv" {
                local status_str "WARN"
            }
            if (`nh0' < 30 & `nh0' < .) | (`nh1' < 30 & `nh1' < .) {
                local status_str "WARN"
            }
            if `has_mass' {
                local mratio = `massinfo'[`idx', 4]
                if `mratio' > 0.2 & `mratio' < . {
                    local status_str "WARN"
                }
            }

            // Format condition numbers
            local cp0_str = cond(`cp0' < ., string(`cp0', "%9.1f"), "       .")
            local cp1_str = cond(`cp1' < ., string(`cp1', "%9.1f"), "       .")
            local cq0_str = cond(`cq0' < ., string(`cq0', "%9.1f"), "       .")
            local cq1_str = cond(`cq1' < ., string(`cq1', "%9.1f"), "       .")

            if "`status_str'" == "WARN" {
                di as err " " %3.0f `idx' _col(8) %5.0f `nh0' _col(15) %5.0f `nh1' _col(22) "`cp0_str'" _col(32) "`cp1_str'" _col(42) "`cq0_str'" _col(52) "`cq1_str'" _col(62) "`fb_str'" _col(72) "`mr_str'" _col(80) "`status_str'"
            }
            else {
                di as res " " %3.0f `idx' _col(8) %5.0f `nh0' _col(15) %5.0f `nh1' _col(22) "`cp0_str'" _col(32) "`cp1_str'" _col(42) "`cq0_str'" _col(52) "`cq1_str'" _col(62) "`fb_str'" _col(72) "`mr_str'" _col(80) "`status_str'"
            }
        }

        di as txt "`hdr_line'"
    }

    // =========================================================================
    // 6. Full mode
    // =========================================================================
    if "`output'" == "full" {
        foreach idx of numlist `indices' {
            local b1 = `diag'[`idx', 1]
            local b2 = `diag'[`idx', 2]
            di as txt "{hline 68}"
            di as txt "Evaluation point `idx'" _col(30) "b1 = " as res %8.4f `b1' as txt "  b2 = " as res %8.4f `b2'
            di as txt "{hline 68}"

            // Effective N
            local nh0 = `results'[`idx', 17]
            local nh1 = `results'[`idx', 18]
            di as txt "  Effective N:" _col(25) "Control = " as res %7.0f `nh0' as txt "   Treated = " as res %7.0f `nh1'

            // Masspoints
            if `has_mass' {
                local m_all = `massinfo'[`idx', 1]
                local m0 = `massinfo'[`idx', 2]
                local m1 = `massinfo'[`idx', 3]
                local mratio = `massinfo'[`idx', 4]
                di as txt "  Masspoints:" _col(25) "Total = " as res %5.0f `m_all' as txt "  M0 = " as res %5.0f `m0' as txt "  M1 = " as res %5.0f `m1' as txt "  ratio = " as res %5.3f `mratio'
            }

            di ""
            di as txt "  {ul:Order p (estimation)}" _col(40) "full-rank = `fullrank_p'"
            // Control side p
            local r0p = `diag'[`idx', 3]
            local c0p = `diag'[`idx', 4]
            local f0p = `diag'[`idx', 5]
            local fb0p_str = cond(`f0p' == 1, "pinv", "invsym")
            di as txt "    Control:" _col(20) "rank = " as res %3.0f `r0p' as txt "  cond = " as res %12.2f `c0p' as txt "  method = " as res "`fb0p_str'"

            // Treated side p
            local r1p = `diag'[`idx', 6]
            local c1p = `diag'[`idx', 7]
            local f1p = `diag'[`idx', 8]
            local fb1p_str = cond(`f1p' == 1, "pinv", "invsym")
            di as txt "    Treated:" _col(20) "rank = " as res %3.0f `r1p' as txt "  cond = " as res %12.2f `c1p' as txt "  method = " as res "`fb1p_str'"

            di ""
            di as txt "  {ul:Order q (bias-correction)}" _col(40) "full-rank = `fullrank_q'"
            // Control side q
            local r0q = `diag'[`idx', 9]
            local c0q = `diag'[`idx', 10]
            local f0q = `diag'[`idx', 11]
            local fb0q_str = cond(`f0q' == 1, "pinv", "invsym")
            di as txt "    Control:" _col(20) "rank = " as res %3.0f `r0q' as txt "  cond = " as res %12.2f `c0q' as txt "  method = " as res "`fb0q_str'"

            // Treated side q
            local r1q = `diag'[`idx', 12]
            local c1q = `diag'[`idx', 13]
            local f1q = `diag'[`idx', 14]
            local fb1q_str = cond(`f1q' == 1, "pinv", "invsym")
            di as txt "    Treated:" _col(20) "rank = " as res %3.0f `r1q' as txt "  cond = " as res %12.2f `c1q' as txt "  method = " as res "`fb1q_str'"

            di ""
        }
    }

    // =========================================================================
    // 7. Warnings mode
    // =========================================================================
    if "`output'" == "warnings" {
        local any_shown = 0

        foreach idx of numlist `indices' {
            local point_warnings ""

            // Check condition numbers
            local cp0 = `diag'[`idx', 4]
            local cp1 = `diag'[`idx', 7]
            local cq0 = `diag'[`idx', 10]
            local cq1 = `diag'[`idx', 13]
            local maxc = 0
            foreach cv in `cp0' `cp1' `cq0' `cq1' {
                if `cv' < . & `cv' > `maxc' {
                    local maxc = `cv'
                }
            }
            if `maxc' > 1e6 {
                local point_warnings "`point_warnings' cond"
            }

            // Check fallback
            local fp0 = `diag'[`idx', 5]
            local fp1 = `diag'[`idx', 8]
            local fq0 = `diag'[`idx', 11]
            local fq1 = `diag'[`idx', 14]
            if `fp0' == 1 | `fp1' == 1 | `fq0' == 1 | `fq1' == 1 {
                local point_warnings "`point_warnings' fallback"
            }

            // Check effective N
            local nh0 = `results'[`idx', 17]
            local nh1 = `results'[`idx', 18]
            if (`nh0' < 30 & `nh0' < .) | (`nh1' < 30 & `nh1' < .) {
                local point_warnings "`point_warnings' lowN"
            }

            // Check mass ratio
            if `has_mass' {
                local mratio = `massinfo'[`idx', 4]
                if `mratio' > 0.2 & `mratio' < . {
                    local point_warnings "`point_warnings' mass"
                }
            }

            // Display if warnings exist
            if "`point_warnings'" != "" {
                if `any_shown' == 0 {
                    di as txt "{hline 68}"
                    di as txt "{bf:Diagnostic warnings}"
                    di as txt "{hline 68}"
                }
                local any_shown = 1
                local b1 = `diag'[`idx', 1]
                local b2 = `diag'[`idx', 2]
                di ""
                di as err "  Point `idx' (b1=" %7.4f `b1' ", b2=" %7.4f `b2' "):"

                if strpos("`point_warnings'", "cond") {
                    di as txt "    {c -} Condition number > 1e6 (max = " as res string(`maxc', "%12.2e") as txt ")"
                    di as txt "      Suggestion: consider {bf:stdvars} option or larger bandwidth"
                }
                if strpos("`point_warnings'", "fallback") {
                    di as txt "    {c -} Rank-deficient design matrix (pinv fallback used)"
                    di as txt "      Suggestion: increase bandwidth or reduce polynomial order"
                }
                if strpos("`point_warnings'", "lowN") {
                    di as txt "    {c -} Very few effective observations (Nh0=" as res %4.0f `nh0' as txt ", Nh1=" as res %4.0f `nh1' as txt ")"
                    di as txt "      Suggestion: increase bandwidth"
                }
                if strpos("`point_warnings'", "mass") {
                    local mratio = `massinfo'[`idx', 4]
                    di as txt "    {c -} High mass-point ratio (" as res string(`mratio', "%5.3f") as txt " > 0.2)"
                    di as txt "      Suggestion: consider {bf:masspoints(adjust)} option"
                }
            }
        }

        if `any_shown' == 0 {
            di as txt "{hline 68}"
            di as txt "{bf:No diagnostic warnings detected.} All evaluation points pass."
            di as txt "{hline 68}"
        }
        else {
            di ""
            di as txt "{hline 68}"
        }
    }

    // =========================================================================
    // 8. Return values
    // =========================================================================
    return scalar max_cond = `max_cond'
    return scalar n_fallback = `n_fallback'
    return scalar min_effN = `min_effN'
    return scalar has_warnings = `has_warnings'
end
