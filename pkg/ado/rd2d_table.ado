*! version 1.1.0 24jun2026
* rd2d_table: Post-estimation command for esttab/estout integration
* Repacks rd2d or rd2d_dist results into eclass (b, V) for tabulation tools.

program define rd2d_table, eclass
    version 16.0
    syntax [using/] [, REPlace FORmat(string) ESTimate(string) ///
           SUBset(numlist integer >0) TEX CSV ///
           NOStar LEVEL(real -1) STATs(string) ///
           MLAbels(string) TITLE(string)]

    // =========================================================================
    // 1. Validation: require rd2d or rd2d_dist estimation results
    // =========================================================================
    if "`e(cmd)'" != "rd2d" & "`e(cmd)'" != "rd2d_dist" {
        di as err "rd2d_table requires estimation results from {bf:rd2d} or {bf:rd2d_dist}"
        exit 301
    }

    // =========================================================================
    // 2. Capture estimation context before ereturn post clears it
    // =========================================================================
    local cmd "`e(cmd)'"
    local neval = e(neval)
    local depvar "`e(depvar)'"
    local vce_type "`e(vce)'"
    local N_total = e(N)
    local p_order = e(p)
    local q_order = e(q)
    local stored_level = e(level)

    // Capture optional scalars that may exist
    capture local N0_val = e(N0)
    if _rc local N0_val = .
    capture local N1_val = e(N1)
    if _rc local N1_val = .

    // Extract results matrix
    tempname results
    matrix `results' = e(results)

    // Attempt to capture covariance matrix
    local has_cov_q 0
    local has_cov_p 0
    tempname covq covp
    capture confirm matrix e(cov_q)
    if !_rc {
        matrix `covq' = e(cov_q)
        local has_cov_q 1
    }
    capture confirm matrix e(cov_p)
    if !_rc {
        matrix `covp' = e(cov_p)
        local has_cov_p 1
    }

    // =========================================================================
    // 3. Parse options
    // =========================================================================

    // estimate(): default "q" (bias-corrected)
    if "`estimate'" == "" local estimate "q"
    local estimate = lower("`estimate'")
    if "`estimate'" != "q" & "`estimate'" != "p" {
        di as err "estimate() must be {bf:q} (bias-corrected) or {bf:p} (conventional); got '`estimate''"
        exit 198
    }

    // format(): default %9.4f
    if "`format'" == "" local format "%9.4f"

    // level(): inherit from estimation if not specified
    if `level' == -1 {
        local level = `stored_level'
    }
    if `level' <= 0 | `level' >= 100 {
        di as err "level() must be between 0 and 100"
        exit 198
    }

    // =========================================================================
    // 4. Determine evaluation point indices
    // =========================================================================
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
    local nuse : word count `indices'

    // =========================================================================
    // 5. Determine matrix columns for selected estimate order
    // =========================================================================
    // e(results) columns for rd2d:
    //   1:b1  2:b2  3:Est_p  4:Se_p  5:Est_q  6:Se_q
    //   7:z   8:pval  9:CI_l  10:CI_u  11:CB_l  12:CB_u
    //   13:h01  14:h02  15:h11  16:h12  17:Nh0  18:Nh1
    if "`estimate'" == "q" {
        local col_est = 5
        local col_se  = 6
    }
    else {
        local col_est = 3
        local col_se  = 4
    }

    // =========================================================================
    // 6. Build coefficient vector and variance matrix
    // =========================================================================
    tempname b V

    matrix `b' = J(1, `nuse', 0)
    matrix `V' = J(`nuse', `nuse', 0)

    // Determine if boundary coordinates are available
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

    // Fill b and diagonal of V; construct coefficient names
    local names ""
    local j = 1
    foreach idx of numlist `indices' {
        matrix `b'[1, `j'] = `results'[`idx', `col_est']
        matrix `V'[`j', `j'] = `results'[`idx', `col_se']^2

        // Create readable coefficient name from boundary coordinates
        if `has_bdy' {
            local b1_raw = `results'[`idx', 1]
            local b2_raw = `results'[`idx', 2]
            local b1_str = string(`b1_raw', "%9.3g")
            local b2_str = string(`b2_raw', "%9.3g")
            // Clean up spaces for valid Stata names
            local b1_str = subinstr("`b1_str'", " ", "", .)
            local b2_str = subinstr("`b2_str'", " ", "", .)
            local b1_str = subinstr("`b1_str'", "-", "n", .)
            local b2_str = subinstr("`b2_str'", "-", "n", .)
            local b1_str = subinstr("`b1_str'", ".", "d", .)
            local b2_str = subinstr("`b2_str'", ".", "d", .)
            local cname "tau_`b1_str'_`b2_str'"
        }
        else {
            local cname "tau_`idx'"
        }
        local names "`names' `cname'"

        local j = `j' + 1
    }

    // Use full covariance matrix if available
    if "`estimate'" == "q" & `has_cov_q' {
        // Extract submatrix for selected indices
        local ri = 1
        foreach ridx of numlist `indices' {
            local ci = 1
            foreach cidx of numlist `indices' {
                matrix `V'[`ri', `ci'] = `covq'[`ridx', `cidx']
                local ci = `ci' + 1
            }
            local ri = `ri' + 1
        }
    }
    else if "`estimate'" == "p" & `has_cov_p' {
        local ri = 1
        foreach ridx of numlist `indices' {
            local ci = 1
            foreach cidx of numlist `indices' {
                matrix `V'[`ri', `ci'] = `covp'[`ridx', `cidx']
                local ci = `ci' + 1
            }
            local ri = `ri' + 1
        }
    }

    // Apply matrix names
    matrix colnames `b' = `names'
    matrix rownames `V' = `names'
    matrix colnames `V' = `names'

    // =========================================================================
    // 7. Post as eclass estimation
    // =========================================================================
    ereturn post `b' `V', obs(`N_total')

    ereturn local cmd "rd2d_table"
    ereturn local cmd_source "`cmd'"
    ereturn local depvar "`depvar'"
    ereturn local vce "`vce_type'"
    ereturn local estimate_order "`estimate'"
    ereturn scalar neval = `nuse'
    ereturn scalar p = `p_order'
    ereturn scalar q = `q_order'
    ereturn scalar level = `level'

    if `N0_val' < . {
        ereturn scalar N0 = `N0_val'
    }
    if `N1_val' < . {
        ereturn scalar N1 = `N1_val'
    }

    // Store the original results submatrix for reference
    tempname results_sub
    matrix `results_sub' = J(`nuse', colsof(`results'), 0)
    local j = 1
    foreach idx of numlist `indices' {
        forvalues c = 1/`=colsof(`results')' {
            matrix `results_sub'[`j', `c'] = `results'[`idx', `c']
        }
        local j = `j' + 1
    }
    ereturn matrix results = `results_sub'

    // =========================================================================
    // 8. Export or display
    // =========================================================================
    if "`using'" != "" {
        // Check esttab availability
        capture which esttab
        if _rc {
            di as err "esttab/estout not found; install via: {stata ssc install estout}"
            di as err "Alternatively, use {bf:rd2d_table} without {bf:using} to display results"
            exit 199
        }

        // Build esttab options
        local esttab_opts `"cells(b(fmt(`format')) se(par fmt(`format')))"'
        if "`nostar'" != "" {
            local esttab_opts `"`esttab_opts' nostar"'
        }
        if "`title'" != "" {
            local esttab_opts `"`esttab_opts' title("`title'")"'
        }

        // File format
        if "`tex'" != "" {
            local esttab_opts `"`esttab_opts' booktabs"'
            esttab . using "`using'", `esttab_opts' `replace'
        }
        else if "`csv'" != "" {
            esttab . using "`using'", `esttab_opts' csv `replace'
        }
        else {
            // Default: infer from extension
            local ext = substr("`using'", -4, .)
            if "`ext'" == ".tex" {
                local esttab_opts `"`esttab_opts' booktabs"'
            }
            else if "`ext'" == ".csv" {
                local esttab_opts `"`esttab_opts' csv"'
            }
            esttab . using "`using'", `esttab_opts' `replace'
        }

        di as txt "  -> table written to `using'"
    }
    else {
        // Display using ereturn display
        di ""
        di as txt "{hline 70}"
        if "`estimate'" == "q" {
            di as txt "rd2d Table: Bias-corrected estimates (order q = `q_order')"
        }
        else {
            di as txt "rd2d Table: Conventional estimates (order p = `p_order')"
        }
        di as txt "Source estimation: `cmd' | VCE: `vce_type' | " ///
            "Dep. var: `depvar'"
        di as txt "{hline 70}"
        ereturn display, level(`level')
    }
end
