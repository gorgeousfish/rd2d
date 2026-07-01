*! version 1.1.0 24jun2026
program define rd2d_dist, eclass
    version 16.0
    if `"`1'"' == "version" | `"`1'"' == ",version" {
        di as txt "rd2d_dist version 1.1.0 (24 June 2026)"
        exit
    }
    syntax varlist(min=2 numeric) [if] [in], ///
        [H(string asis) P(integer 1) Q(string) KINK(string) KERnel(string) ///
         LEVEL(real 95) SIDE(string) BWSELect(string) VCE(string) RBC(string) ///
         BWCHeck(string) MASSPoints(string) SCALEregul(real 1) CQT(real 0.5) ///
         CBANDs NOCBANDs REPP(integer 1000) CLuster(varname) ///
         MP(string) SR(real -1) BWS(string)]

    gettoken yvar dvars : varlist
    local neval : word count `dvars'
    if (`neval' == 0) {
        di as err "at least one distance variable is required after the outcome variable"
        exit 198
    }

    // --- cbands default logic (v1.2.0: default ON, nocbands to disable) ---
    local cbands_on = ("`nocbands'" == "")

    local clustername "`cluster'"
    local clusterwork "`cluster'"
    marksample touse
    markout `touse' `varlist'
    if ("`clustername'" != "") {
        markout `touse' `clustername', strok
        tempvar clusterid
        quietly egen long `clusterid' = group(`clustername') if `touse'
        local clusterwork "`clusterid'"
    }

    // --- Option alias resolution (4.4) ---
    if ("`mp'" != "") {
        if ("`masspoints'" != "") {
            di as err "cannot specify both mp() and masspoints()"
            exit 198
        }
        local masspoints "`mp'"
    }
    if (`sr' != -1) {
        if (`scaleregul' != 1) {
            di as err "cannot specify both sr() and scaleregul()"
            exit 198
        }
        local scaleregul = `sr'
    }
    if ("`bws'" != "") {
        if ("`bwselect'" != "") {
            di as err "cannot specify both bws() and bwselect()"
            exit 198
        }
        local bwselect "`bws'"
    }

    if (`p' < 0) {
        di as err "p() must be a positive integer (typically 1, 2, or 3); got `p'"
        di as err "  p=1 (local linear) is recommended for most applications"
        exit 198
    }
    local qeff -1
    if ("`q'" != "") {
        capture confirm integer number `q'
        if (_rc | `q' < 0 | `q' >= .) {
            di as err "q() must be a nonnegative integer >= p (currently p=`p'); got `q'"
            di as err "  q=p+1 is the default for bias correction"
            exit 198
        }
        local qeff = `q'
    }
    if (`level' <= 0 | `level' >= 100) {
        di as err "level() must be between 0 and 100 (e.g., 95 for 95% confidence)"
        exit 198
    }
    if (`repp' < 1) {
        di as err "repp() must be a positive integer"
        exit 198
    }
    if (`scaleregul' >= . | `scaleregul' < 0) {
        di as err "scaleregul() must be a non-negative number (default: 1 for distance)"
        exit 198
    }
    if (`cqt' <= 0 | `cqt' >= 1) {
        di as err "cqt() must be between 0 and 1"
        exit 198
    }

    local kernel = lower("`kernel'")
    if ("`kernel'" == "") local kernel "triangular"
    if inlist("`kernel'", "uni", "unif") local kernel "uniform"
    if inlist("`kernel'", "tri", "triag") local kernel "triangular"
    if inlist("`kernel'", "epa", "epan") local kernel "epanechnikov"
    if inlist("`kernel'", "gau") local kernel "gaussian"
    if !inlist("`kernel'", "uniform", "triangular", "epanechnikov", "gaussian") {
        di as err "kernel() must be one of: {bf:triangular}, {bf:epanechnikov}, {bf:uniform}, {bf:gaussian}"
        di as err "  (abbreviations: tri, epa, uni, gau); got '`kernel''"
        exit 198
    }

    // bwselect options:
    //   mserd   = MSE-optimal common bandwidth (minimizes combined MSE across sides)
    //   imserd  = integrated MSE-optimal common bandwidth (averages MSE across eval points)
    //   msetwo  = MSE-optimal separate bandwidths for control and treated sides
    //   imsetwo = integrated MSE-optimal separate bandwidths
    local bwselect = lower("`bwselect'")
    if ("`bwselect'" == "") local bwselect "mserd"
    if !inlist("`bwselect'", "mserd", "imserd", "msetwo", "imsetwo") {
        di as err "bwselect() must be one of: {bf:mserd}, {bf:imserd}, {bf:msetwo}, {bf:imsetwo}"
        di as err "  mserd = MSE-optimal common bandwidth; msetwo = side-specific bandwidths"
        exit 198
    }

    local kink = lower("`kink'")
    if ("`kink'" == "") local kink "off"
    if !inlist("`kink'", "off", "on") {
        di as err "kink() must be {bf:on} or {bf:off}; kink(on) uses undersmoothing instead of RBC"
        exit 198
    }

    local vce = lower("`vce'")
    if ("`vce'" == "") local vce "hc1"
    if !inlist("`vce'", "hc0", "hc1", "hc2", "hc3") {
        di as err "vce() must be one of: {bf:hc0}, {bf:hc1}, {bf:hc2}, {bf:hc3}"
        di as err "  (hc1 is recommended for small-sample adjustment)"
        exit 198
    }
    if ("`clustername'" != "" & !inlist("`vce'", "hc0", "hc1")) {
        di as txt "{bf:note}: vce(`vce') is not available with cluster(); using hc1 instead"
        di as txt "  Reason: leverage-based HC2/HC3 corrections require observation-level structure"
        local vce "hc1"
    }

    local rbc = lower("`rbc'")
    if ("`rbc'" == "") local rbc "on"
    if !inlist("`rbc'", "on", "off") {
        di as err "rbc() must be {bf:on} (robust bias-corrected) or {bf:off} (conventional)"
        exit 198
    }

    local side = lower("`side'")
    if ("`side'" == "") local side "two"
    if !inlist("`side'", "two", "left", "right") {
        di as err "side() must be {bf:two} (two-sided), {bf:left}, or {bf:right}"
        exit 198
    }

    local masspoints = lower("`masspoints'")
    if ("`masspoints'" == "") local masspoints "check"
    if !inlist("`masspoints'", "check", "adjust", "off") {
        di as err "masspoints() must be one of: {bf:check}, {bf:adjust}, {bf:off}"
        di as err "  check = warn if mass points detected; adjust = adapt bandwidth"
        exit 198
    }

    // h() bandwidth parsing: accepts 1 (common), 2 (side-specific), or 2*neval values
    local hcount 0
    local hlist ""
    if ("`h'" != "") {
        local hlist "`h'"
        local hcount : word count `hlist'
        forvalues i = 1/`hcount' {
            local hv : word `i' of `hlist'
            capture confirm number `hv'
            if (_rc) {
                di as err "h() bandwidths must be positive numbers; element `i' is not a number"
                exit 198
            }
            if (`hv' >= . | `hv' <= 0) {
                di as err "h() bandwidths must be finite positive numbers; element `i' has value `hv'"
                exit 198
            }
        }
    }

    local bwcheck_raw "`bwcheck'"
    if ("`bwcheck_raw'" == "") {
        local bwcheck = 50 + `p' + 1
    }
    else {
        capture confirm integer number `bwcheck_raw'
        if (_rc | `bwcheck_raw' < 0) {
            di as err "bwcheck() must be a nonnegative integer"
            exit 198
        }
        local bwcheck = `bwcheck_raw'
    }
    quietly count if `touse'
    local N = r(N)
    if (`N' == 0) {
        di as err "no observations"
        exit 2000
    }
    if (`bwcheck' > 0 & `N' < `bwcheck') {
        di as err "Insufficient observations (N=`N') for bandwidth estimation. Minimum required: `bwcheck' (bwcheck option). Consider reducing bwcheck() or providing manual h()."
        exit 2001
    }

    // Distance variable sign convention: D<0 = control side, D>=0 = treated side.
    // All distance variables must share the same sign pattern (same obs on same side).
    if (`neval' > 1) {
        local firstd : word 1 of `dvars'
        forvalues j = 2/`neval' {
            local dvar : word `j' of `dvars'
            quietly count if `touse' & ((`dvar' >= 0) != (`firstd' >= 0))
            if (r(N) > 0) {
                di as err "distance columns must use a common treatment-side sign pattern (D<0=control, D>=0=treated)"
                di as err "  All distance variables must classify the same observations as control/treated"
                exit 198
            }
        }
    }

    if ("`rbc'" == "off") local qeff = `p'
    else if (`qeff' < 0) local qeff = cond("`kink'" == "on", `p', `p' + 1)
    if (`qeff' < `p') {
        di as err "q() must be an integer >= p (currently p=`p'); got q=`qeff'"
        di as err "  q=p+1 is the default for bias correction"
        exit 198
    }
    // kink(on): boundary kink degrades bias from h^{p+1} to h, making
    // standard RBC (q=p+1) invalid. Use undersmoothing (q=p) instead.
    if ("`kink'" == "on" & "`rbc'" == "on" & `qeff' > `p') {
        di as txt "warning: q() reset to p() because kink(on) uses undersmoothing, not RBC."
        local qeff = `p'
    }

    local side_min = `qeff' + 1
    forvalues j = 1/`neval' {
        local dvar : word `j' of `dvars'
        quietly count if `touse' & `dvar' < 0
        local N0j = r(N)
        quietly count if `touse' & `dvar' >= 0
        local N1j = r(N)
        if (`j' == 1) {
            local N0 = `N0j'
            local N1 = `N1j'
        }
        if (`N0j' < `side_min' | `N1j' < `side_min') {
            local failside = cond(`N0j' < `side_min', "control", "treated")
            local failnh = cond(`N0j' < `side_min', `N0j', `N1j')
            di as err "Insufficient observations on `failside' side (Nh=`failnh', need `side_min') at distance column `j'"
            di as err "  Suggestion: increase bandwidth with h() or reduce polynomial order p()"
            exit 2001
        }
    }

    tempname bws massinfo results resultsA0 resultsA1 bwsout hpmat hqmat

    if (`hcount' > 0) {
        local expected_h = 2 * `neval'
        if !inlist(`hcount', 1, 2, `expected_h') {
            di as err "h() requires 1, 2, or 2*neval values; got `hcount' for `neval' distance variables"
            exit 198
        }
        matrix `bws' = J(`neval', 6, .)
        // massinfo matrix columns: [M, M0, M1, mass_ratio]
        // M = total unique distance values, M0/M1 = unique per side
        // mass_ratio = 1 - M/N; values >= 0.2 trigger warnings
        matrix `massinfo' = J(`neval', 4, .)
        forvalues j = 1/`neval' {
            local dvar : word `j' of `dvars'
            if (`hcount' == 1) {
                local h0v : word 1 of `hlist'
                local h1v : word 1 of `hlist'
            }
            else if (`hcount' == 2) {
                local h0v : word 1 of `hlist'
                local h1v : word 2 of `hlist'
            }
            else {
                local h0i = 2 * (`j' - 1) + 1
                local h1i = 2 * (`j' - 1) + 2
                local h0v : word `h0i' of `hlist'
                local h1v : word `h1i' of `hlist'
            }
            matrix `bws'[`j', 3] = `h0v'
            matrix `bws'[`j', 4] = `h1v'
            if ("`masspoints'" != "off") {
                tempvar mass_all mass_0 mass_1
                quietly egen long `mass_all' = group(`dvar') if `touse'
                quietly summarize `mass_all' if `touse', meanonly
                local M = r(max)
                quietly egen long `mass_0' = group(`dvar') if `touse' & `dvar' < 0
                quietly summarize `mass_0' if `touse' & `dvar' < 0, meanonly
                local M0 = r(max)
                quietly egen long `mass_1' = group(`dvar') if `touse' & `dvar' >= 0
                quietly summarize `mass_1' if `touse' & `dvar' >= 0, meanonly
                local M1 = r(max)
                local mass = 1 - `M' / `N'
                matrix `massinfo'[`j', 1] = `M'
                matrix `massinfo'[`j', 2] = `M0'
                matrix `massinfo'[`j', 3] = `M1'
                matrix `massinfo'[`j', 4] = `mass'
            }
            else {
                quietly count if `touse' & `dvar' < 0
                local N0j = r(N)
                quietly count if `touse' & `dvar' >= 0
                local N1j = r(N)
                matrix `massinfo'[`j', 1] = `N'
                matrix `massinfo'[`j', 2] = `N0j'
                matrix `massinfo'[`j', 3] = `N1j'
                matrix `massinfo'[`j', 4] = 0
            }
        }
        _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')
        local bwsource "user"
        local bwselect "user provided"
    }
    else {
        local clusteropt ""
        if ("`clusterwork'" != "") local clusteropt "cluster(`clusterwork')"
        quietly rdbw2d_dist `yvar' `dvars' if `touse', p(`p') kernel(`kernel') ///
            bwselect(`bwselect') kink(`kink') vce(`vce') bwcheck(`bwcheck') ///
            masspoints(`masspoints') scaleregul(`scaleregul') cqt(`cqt') `clusteropt'
        matrix `bws' = r(bws)
        matrix `massinfo' = r(masspoints)
        _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')
        local bwsource "automatic"
    }

    // --- Bandwidth zero/negative check (4.3) ---
    forvalues j = 1/`neval' {
        if (`bws'[`j', 3] <= 0 | `bws'[`j', 4] <= 0) {
            di as err "computed bandwidth is zero or negative at point j=`j'; check data density near this evaluation point"
            exit 498
        }
    }

    tempname diagnostics
    matrix `results' = J(`neval', 18, .)
    matrix `resultsA0' = J(`neval', 9, .)
    matrix `resultsA1' = J(`neval', 9, .)
    matrix `bwsout' = J(`neval', 6, .)
    matrix `hpmat' = J(`neval', 2, .)
    matrix `hqmat' = J(`neval', 2, .)
    matrix `diagnostics' = J(`neval', 14, .)
    local anyfallback 0

    tempname mu0p mu1p se0p se1p mu0q mu1q se0q se1q h0 h1 h0rbc h1rbc
    tempname taup sep tauq seq zval pvalue crit cil ciu nh0 nh1 nh0q nh1q
    tempname fit0 fit1 fit0q fit1q
    tempname r0p c0p f0p r1p c1p f1p r0q c0q f0q r1q c1q f1q

    forvalues j = 1/`neval' {
        local dvar : word `j' of `dvars'
        scalar `h0' = `bws'[`j', 3]
        scalar `h1' = `bws'[`j', 4]
        scalar `h0rbc' = `h0'
        scalar `h1rbc' = `h1'

        if ("`rbc'" == "on" & "`kink'" == "on" & "`bwsource'" == "automatic") {
            if inlist("`bwselect'", "mserd", "imserd") {
                scalar `h0rbc' = `h0' * `massinfo'[`j', 1]^(-1/3) / `massinfo'[`j', 1]^(-1/4)
                scalar `h1rbc' = `h1' * `massinfo'[`j', 1]^(-1/3) / `massinfo'[`j', 1]^(-1/4)
            }
            else {
                scalar `h0rbc' = `h0' * `massinfo'[`j', 2]^(-1/3) / `massinfo'[`j', 2]^(-1/4)
                scalar `h1rbc' = `h1' * `massinfo'[`j', 3]^(-1/3) / `massinfo'[`j', 3]^(-1/4)
            }
        }

        mata: st_matrix("`fit0'", _rd2d_dist_fit_raw("`yvar'", "`dvar'", "`touse'", ///
            "`clusterwork'", `=`h0'', `p', "control", "`kernel'", "`vce'"))
        scalar `mu0p' = `fit0'[1,1]
        scalar `se0p' = `fit0'[1,2]
        scalar `nh0' = `fit0'[1,3]
        scalar `r0p' = `fit0'[1,4]
        scalar `c0p' = `fit0'[1,5]
        scalar `f0p' = `fit0'[1,6]

        mata: st_matrix("`fit1'", _rd2d_dist_fit_raw("`yvar'", "`dvar'", "`touse'", ///
            "`clusterwork'", `=`h1'', `p', "treated", "`kernel'", "`vce'"))
        scalar `mu1p' = `fit1'[1,1]
        scalar `se1p' = `fit1'[1,2]
        scalar `nh1' = `fit1'[1,3]
        scalar `r1p' = `fit1'[1,4]
        scalar `c1p' = `fit1'[1,5]
        scalar `f1p' = `fit1'[1,6]

        mata: st_matrix("`fit0q'", _rd2d_dist_fit_raw("`yvar'", "`dvar'", "`touse'", ///
            "`clusterwork'", `=`h0rbc'', `qeff', "control", "`kernel'", "`vce'"))
        scalar `mu0q' = `fit0q'[1,1]
        scalar `se0q' = `fit0q'[1,2]
        scalar `nh0q' = `fit0q'[1,3]
        scalar `r0q' = `fit0q'[1,4]
        scalar `c0q' = `fit0q'[1,5]
        scalar `f0q' = `fit0q'[1,6]

        mata: st_matrix("`fit1q'", _rd2d_dist_fit_raw("`yvar'", "`dvar'", "`touse'", ///
            "`clusterwork'", `=`h1rbc'', `qeff', "treated", "`kernel'", "`vce'"))
        scalar `mu1q' = `fit1q'[1,1]
        scalar `se1q' = `fit1q'[1,2]
        scalar `nh1q' = `fit1q'[1,3]
        scalar `r1q' = `fit1q'[1,4]
        scalar `c1q' = `fit1q'[1,5]
        scalar `f1q' = `fit1q'[1,6]

        if (`f0p' != 0 | `f1p' != 0 | `f0q' != 0 | `f1q' != 0) {
            local anyfallback 1
            if (`f0p' != 0) di as txt "{bf:note}: design matrix rank-deficient at point j=`j' control-side p-fit" _newline ///
                "  rank=" %4.0f `r0p' " cond=" %9.3e `c0p' _newline ///
                "  Using generalized inverse (pinv). Consider: larger h(), stdvars, or lower p()"
            if (`f1p' != 0) di as txt "{bf:note}: design matrix rank-deficient at point j=`j' treated-side p-fit" _newline ///
                "  rank=" %4.0f `r1p' " cond=" %9.3e `c1p' _newline ///
                "  Using generalized inverse (pinv). Consider: larger h(), stdvars, or lower p()"
            if (`f0q' != 0) di as txt "{bf:note}: design matrix rank-deficient at point j=`j' control-side q-fit" _newline ///
                "  rank=" %4.0f `r0q' " cond=" %9.3e `c0q' _newline ///
                "  Using generalized inverse (pinv). Consider: larger h(), stdvars, or lower p()"
            if (`f1q' != 0) di as txt "{bf:note}: design matrix rank-deficient at point j=`j' treated-side q-fit" _newline ///
                "  rank=" %4.0f `r1q' " cond=" %9.3e `c1q' _newline ///
                "  Using generalized inverse (pinv). Consider: larger h(), stdvars, or lower p()"
        }

        scalar `taup' = `mu1p' - `mu0p'
        scalar `sep' = sqrt(`se0p'^2 + `se1p'^2)
        scalar `tauq' = `mu1q' - `mu0q'
        scalar `seq' = sqrt(`se0q'^2 + `se1q'^2)
        if (`seq' > 0 & `seq' < .) {
            scalar `zval' = `tauq' / `seq'
            scalar `pvalue' = 2 * normal(-abs(`zval'))
        }
        else {
            scalar `zval' = .
            scalar `pvalue' = .
        }

        if ("`side'" == "two") {
            scalar `crit' = invnormal((100 + `level') / 200)
            scalar `cil' = `tauq' - `crit' * `seq'
            scalar `ciu' = `tauq' + `crit' * `seq'
        }
        else if ("`side'" == "left") {
            scalar `crit' = invnormal(`level' / 100)
            scalar `cil' = -c(maxdouble)
            scalar `ciu' = `tauq' + `crit' * `seq'
        }
        else {
            scalar `crit' = invnormal(`level' / 100)
            scalar `cil' = `tauq' - `crit' * `seq'
            scalar `ciu' = c(maxdouble)
        }

        matrix `results'[`j', 1] = `bws'[`j', 1]
        matrix `results'[`j', 2] = `bws'[`j', 2]
        matrix `results'[`j', 3] = `taup'
        matrix `results'[`j', 4] = `sep'
        matrix `results'[`j', 5] = `tauq'
        matrix `results'[`j', 6] = `seq'
        matrix `results'[`j', 7] = `zval'
        matrix `results'[`j', 8] = `pvalue'
        matrix `results'[`j', 9] = `cil'
        matrix `results'[`j', 10] = `ciu'
        matrix `results'[`j', 13] = `h0'
        matrix `results'[`j', 14] = `h1'
        matrix `results'[`j', 15] = `h0rbc'
        matrix `results'[`j', 16] = `h1rbc'
        matrix `results'[`j', 17] = `nh0'
        matrix `results'[`j', 18] = `nh1'

        matrix `resultsA0'[`j', 1] = `bws'[`j', 1]
        matrix `resultsA0'[`j', 2] = `bws'[`j', 2]
        matrix `resultsA0'[`j', 3] = `mu0p'
        matrix `resultsA0'[`j', 4] = `se0p'
        matrix `resultsA0'[`j', 5] = `mu0q'
        matrix `resultsA0'[`j', 6] = `se0q'
        matrix `resultsA0'[`j', 7] = `h0'
        matrix `resultsA0'[`j', 8] = `h0rbc'
        matrix `resultsA0'[`j', 9] = `nh0'

        matrix `resultsA1'[`j', 1] = `bws'[`j', 1]
        matrix `resultsA1'[`j', 2] = `bws'[`j', 2]
        matrix `resultsA1'[`j', 3] = `mu1p'
        matrix `resultsA1'[`j', 4] = `se1p'
        matrix `resultsA1'[`j', 5] = `mu1q'
        matrix `resultsA1'[`j', 6] = `se1q'
        matrix `resultsA1'[`j', 7] = `h1'
        matrix `resultsA1'[`j', 8] = `h1rbc'
        matrix `resultsA1'[`j', 9] = `nh1'

        matrix `bwsout'[`j', 1] = `bws'[`j', 1]
        matrix `bwsout'[`j', 2] = `bws'[`j', 2]
        matrix `bwsout'[`j', 3] = `h0'
        matrix `bwsout'[`j', 4] = `h1'
        matrix `bwsout'[`j', 5] = `nh0'
        matrix `bwsout'[`j', 6] = `nh1'
        matrix `hpmat'[`j', 1] = `h0'
        matrix `hpmat'[`j', 2] = `h1'
        matrix `hqmat'[`j', 1] = `h0rbc'
        matrix `hqmat'[`j', 2] = `h1rbc'
        matrix `diagnostics'[`j', 1] = `bws'[`j', 1]
        matrix `diagnostics'[`j', 2] = `bws'[`j', 2]
        matrix `diagnostics'[`j', 3] = `r0p'
        matrix `diagnostics'[`j', 4] = `c0p'
        matrix `diagnostics'[`j', 5] = `f0p'
        matrix `diagnostics'[`j', 6] = `r1p'
        matrix `diagnostics'[`j', 7] = `c1p'
        matrix `diagnostics'[`j', 8] = `f1p'
        matrix `diagnostics'[`j', 9] = `r0q'
        matrix `diagnostics'[`j', 10] = `c0q'
        matrix `diagnostics'[`j', 11] = `f0q'
        matrix `diagnostics'[`j', 12] = `r1q'
        matrix `diagnostics'[`j', 13] = `c1q'
        matrix `diagnostics'[`j', 14] = `f1q'
    }

    tempname covp corrp covq corrq cbcrit cbpsd cbmineig
    scalar `cbcrit' = .
    scalar `cbpsd' = 0
    scalar `cbmineig' = .
    matrix `covq' = J(`neval', `neval', .)
    matrix `corrq' = J(`neval', `neval', .)
    forvalues j = 1/`neval' {
        matrix `covq'[`j', `j'] = `results'[`j', 6]^2
        matrix `corrq'[`j', `j'] = cond(`results'[`j', 6] > 0 & `results'[`j', 6] < ., 1, .)
    }
    if ("`clusterwork'" != "") {
        mata: _rd2d_dist_cov_q("`yvar'", "`dvars'", "`touse'", "`hpmat'", `p', ///
            "`kernel'", "`vce'", "`clusterwork'", "`covp'", "`corrp'")
        forvalues j = 1/`neval' {
            matrix `results'[`j', 4] = sqrt(`covp'[`j', `j'])
        }
    }
    if (`neval' > 1 | `cbands_on' | "`clusterwork'" != "") {
        mata: _rd2d_dist_cov_q("`yvar'", "`dvars'", "`touse'", "`hqmat'", `qeff', ///
            "`kernel'", "`vce'", "`clusterwork'", "`covq'", "`corrq'")
    }
    forvalues j = 1/`neval' {
        matrix `results'[`j', 6] = sqrt(`covq'[`j', `j'])
        if (`results'[`j', 6] > 0 & `results'[`j', 6] < .) {
            matrix `results'[`j', 7] = `results'[`j', 5] / `results'[`j', 6]
            matrix `results'[`j', 8] = 2 * normal(-abs(`results'[`j', 7]))
        }
        else {
            matrix `results'[`j', 7] = .
            matrix `results'[`j', 8] = .
        }
        if (`results'[`j', 15] == `results'[`j', 13] & ///
            `results'[`j', 16] == `results'[`j', 14] & `qeff' == `p') {
            matrix `results'[`j', 4] = `results'[`j', 6]
        }
        if ("`side'" == "two") {
            matrix `results'[`j', 9] = `results'[`j', 5] - invnormal((100 + `level') / 200) * `results'[`j', 6]
            matrix `results'[`j', 10] = `results'[`j', 5] + invnormal((100 + `level') / 200) * `results'[`j', 6]
        }
        else if ("`side'" == "left") {
            matrix `results'[`j', 9] = -c(maxdouble)
            matrix `results'[`j', 10] = `results'[`j', 5] + invnormal(`level' / 100) * `results'[`j', 6]
        }
        else {
            matrix `results'[`j', 9] = `results'[`j', 5] - invnormal(`level' / 100) * `results'[`j', 6]
            matrix `results'[`j', 10] = c(maxdouble)
        }
    }
    if (`cbands_on') {
        forvalues j = 1/`neval' {
            if (`results'[`j', 6] <= 0 | `results'[`j', 6] >= .) {
                di as err "unable to compute cbands with nonpositive standard error"
                exit 498
            }
        }
    }
    if (`cbands_on') {
        mata: _rd2d_dist_cb_from_corr("`corrq'", `repp', `level', "`side'", ///
            "`cbcrit'", "`cbpsd'", "`cbmineig'")
        forvalues j = 1/`neval' {
            if ("`side'" == "two") {
                matrix `results'[`j', 11] = `results'[`j', 5] - `cbcrit' * `results'[`j', 6]
                matrix `results'[`j', 12] = `results'[`j', 5] + `cbcrit' * `results'[`j', 6]
            }
            else if ("`side'" == "left") {
                matrix `results'[`j', 11] = -c(maxdouble)
                matrix `results'[`j', 12] = `results'[`j', 5] + `cbcrit' * `results'[`j', 6]
            }
            else {
                matrix `results'[`j', 11] = `results'[`j', 5] - `cbcrit' * `results'[`j', 6]
                matrix `results'[`j', 12] = c(maxdouble)
            }
        }
    }

    matrix colnames `results' = b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h0 h1 h0_rbc h1_rbc Nh0 Nh1
    matrix colnames `resultsA0' = b1 b2 Est_p Se_p Est_q Se_q h0 h0_rbc Nh0
    matrix colnames `resultsA1' = b1 b2 Est_p Se_p Est_q Se_q h1 h1_rbc Nh1
    matrix colnames `bwsout' = b1 b2 h0 h1 Nh0 Nh1
    matrix colnames `diagnostics' = b1 b2 rank_p0 cond_p0 fb_p0 rank_p1 cond_p1 fb_p1 rank_q0 cond_q0 fb_q0 rank_q1 cond_q1 fb_q1
    matrix colnames `massinfo' = M M0 M1 mass
    matrix rownames `results' = `dvars'
    matrix rownames `resultsA0' = `dvars'
    matrix rownames `resultsA1' = `dvars'
    matrix rownames `bwsout' = `dvars'
    matrix rownames `diagnostics' = `dvars'
    matrix rownames `massinfo' = `dvars'
    matrix colnames `covq' = `dvars'
    matrix rownames `covq' = `dvars'
    matrix colnames `corrq' = `dvars'
    matrix rownames `corrq' = `dvars'

    tempname eb eV
    matrix `eb' = J(1, `neval', .)
    forvalues j = 1/`neval' {
        matrix `eb'[1, `j'] = `results'[`j', 5]
    }
    matrix colnames `eb' = `dvars'
    matrix `eV' = `covq'

    ereturn clear
    ereturn post `eb' `eV'
    matrix colnames `results' = b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h0 h1 h0rbc h1rbc Nh0 Nh1
    matrix colnames `resultsA0' = b1 b2 Est_p Se_p Est_q Se_q h0 h0rbc Nh0
    matrix colnames `resultsA1' = b1 b2 Est_p Se_p Est_q Se_q h1 h1rbc Nh1
    matrix colnames `bwsout' = b1 b2 h0 h1 Nh0 Nh1
    matrix colnames `diagnostics' = b1 b2 rank_p0 cond_p0 fb_p0 rank_p1 cond_p1 fb_p1 rank_q0 cond_q0 fb_q0 rank_q1 cond_q1 fb_q1
    matrix colnames `massinfo' = M M0 M1 mass_ratio
    ereturn matrix results = `results'
    ereturn matrix results_A0 = `resultsA0'
    ereturn matrix results_A1 = `resultsA1'
    ereturn matrix bws = `bwsout'
    ereturn matrix diagnostics = `diagnostics'
    ereturn matrix masspoints = `massinfo'
    ereturn matrix cov_q = `covq'
    ereturn matrix corr_q = `corrq'
    ereturn scalar N = `N'
    ereturn scalar N0 = `N0'
    ereturn scalar N1 = `N1'
    ereturn scalar neval = `neval'
    ereturn scalar p = `p'
    ereturn scalar q = `qeff'
    ereturn scalar level = `level'
    ereturn scalar repp = `repp'
    ereturn scalar cb_crit = `cbcrit'
    ereturn scalar cb_psd_adjusted = `cbpsd'
    ereturn scalar cb_min_eigen = `cbmineig'
    ereturn scalar bwcheck = `bwcheck'
    ereturn scalar scaleregul = `scaleregul'
    ereturn scalar cqt = `cqt'
    ereturn local cmd "rd2d_dist"
        ereturn local depvar "`yvar'"
    ereturn local kernel "`kernel'"
    ereturn local bwselect "`bwselect'"
    ereturn local kink "`kink'"
    ereturn local vce "`vce'"
    ereturn local rbc "`rbc'"
    ereturn local side "`side'"
    ereturn local masspoints_opt "`masspoints'"
    ereturn local bwsource "`bwsource'"
    ereturn local cbands = cond(`cbands_on', "on", "off")
    ereturn local clustered = cond("`clustername'" != "", "on", "off")
    ereturn local cluster "`clustername'"
    ereturn local fallback = cond(`anyfallback', "pinv", "invsym")
    ereturn local version "1.1.0"
    ereturn local depvar "`yvar'"
    ereturn local dvars "`dvars'"

    // --- cbands multi-point hint (4.2) ---
    if (`neval' > 1 & !`cbands_on') {
        di as txt "note: consider removing {cmd:nocbands} option for uniform inference over `neval' evaluation points"
    }

    // ===================================================================
    // TABLE DISPLAY SECTION
    // Layout modes: ultra (<50), compact (<63), narrow (<79), normal (>=79)
    // With cbands: shows Est.q, SE.q, CI.lo, CB.lo, CB.hi
    // Without cbands: shows Est.p, SE.p, Est.q, SE.q, CI.lo, CI.hi
    // ===================================================================
    * --- Table layout parameters ---
    local line_width = min(79, c(linesize))
    local hline_rule `"di as txt \"{hline `line_width'}\""'
    local ultra_table = c(linesize) < 50
    local compact_table = c(linesize) < 63
    local narrow_table = c(linesize) < 79
    di as txt _newline "Distance RD estimation"
    di as txt "{hline `line_width'}"
    if (`ultra_table' | `compact_table') {
        di as txt "  Eval points: " as res %9.0f `neval'
        di as txt "  Observations: " as res %9.0f `N'
        di as txt "  VCE: " as res "`vce'" as txt "  RBC: " as res "`rbc'"
        di as txt "  Side: " as res "`side'" as txt "  Kink: " as res "`kink'"
        di as txt "  Bandwidth: " as res "`bwsource'"
        di as txt "  Kernel: " as res "`kernel'"
    }
    else {
        di as txt "  Evaluation points: " as res %9.0f `neval' ///
            as txt "    Observations: " as res %9.0f `N'
        di as txt "  VCE: " as res "`vce'" ///
            as txt "    RBC: " as res "`rbc'" ///
            as txt "    Side: " as res "`side'"
        di as txt "  Bandwidth source: " as res "`bwsource'" ///
            as txt "    Kernel: " as res "`kernel'" ///
            as txt "    Kink: " as res "`kink'"
    }
    di as txt "{hline `line_width'}"
    if (`ultra_table') {
        if (`cbands_on') {
            di as txt %4s "Pt" _col(7) %8s "Est" _col(17) %8s "SE"
            di as txt %4s "" _col(7) %7s "CI<" _col(16) %7s "CI>" ///
                _col(25) %7s "CB<" _col(34) %7s "CB>"
        }
        else {
            di as txt %4s "Pt" _col(7) %8s "Ep" _col(17) %8s "SEp" ///
                _col(27) %8s "Eq"
            di as txt %4s "" _col(7) %7s "SEq" _col(16) %7s "CI<" ///
                _col(25) %7s "CI>"
        }
    }
    else if (`compact_table') {
        if (`cbands_on') {
            di as txt %6s "Point" _col(8) %6s "Est.q" _col(15) %6s "SE.q" ///
                _col(22) %6s "CI.lo" _col(29) %6s "CI.hi" ///
                _col(36) %6s "CB.lo" _col(43) %6s "CB.hi"
        }
        else {
            di as txt %6s "Point" _col(8) %6s "Est.p" _col(15) %6s "SE.p" ///
                _col(22) %6s "Est.q" _col(29) %6s "SE.q" ///
                _col(36) %6s "CI.lo" _col(43) %6s "CI.hi"
        }
    }
    else if (`narrow_table') {
        if (`cbands_on') {
            di as txt %8s "Point" _col(10) %8s "Est.q" _col(19) %8s "SE.q" ///
                _col(28) %8s "CI.low" _col(37) %8s "CI.high" ///
                _col(46) %8s "CB.low" _col(55) %8s "CB.high"
        }
        else {
            di as txt %8s "Point" _col(10) %8s "Est.p" _col(19) %8s "SE.p" ///
                _col(28) %8s "Est.q" _col(37) %8s "SE.q" ///
                _col(46) %8s "CI.low" _col(55) %8s "CI.high"
        }
    }
    else {
        if (`cbands_on') {
            di as txt %10s "Point" _col(12) %10s "Est.q" _col(23) %9s "SE.q" ///
                _col(33) %9s "CI.low" _col(43) %9s "CI.high" ///
                _col(53) %9s "CB.low" _col(63) %9s "CB.high"
        }
        else {
            di as txt %10s "Point" _col(12) %11s "Est.p" _col(24) %11s "SE.p" ///
                _col(36) %11s "Est.q" _col(48) %11s "SE.q" ///
                _col(60) %9s "CI.low" _col(70) %9s "CI.high"
        }
    }
    di as txt "{hline `line_width'}"
    tempname display_results
    matrix `display_results' = e(results)
    forvalues j = 1/`neval' {
        local rname : word `j' of `dvars'
        local suffix : display %02.0f `j'
        local label_width = cond(`ultra_table', 4, cond(`compact_table', 6, cond(`narrow_table', 8, 10)))
        local prefix_len = max(0, `label_width' - strlen("`suffix'") - 1)
        local dname = cond(strlen("`rname'") > `label_width', substr("`rname'", 1, `prefix_len') + "~`suffix'", "`rname'")
        local estp = el(`display_results', `j', 3)
        local sep = el(`display_results', `j', 4)
        local estq = el(`display_results', `j', 5)
        local seq = el(`display_results', `j', 6)
        local cil = el(`display_results', `j', 9)
        local ciu = el(`display_results', `j', 10)
        if (`ultra_table' & `cbands_on') {
            local cbl = el(`display_results', `j', 11)
            local cbu = el(`display_results', `j', 12)
            local cil_u : display %7.3g `cil'
            local ciu_u : display %7.3g `ciu'
            local cbl_u : display %7.3g `cbl'
            local cbu_u : display %7.3g `cbu'
            if (`cil' <= -c(maxdouble) / 2) local cil_u "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_u "inf"
            if (`cbl' <= -c(maxdouble) / 2) local cbl_u "-inf"
            if (`cbu' >= c(maxdouble) / 2) local cbu_u "inf"
            di as res %4s "`dname'" _col(7) %8.3g `estq' ///
                _col(17) %8.3g `seq'
            di as res %4s "" _col(7) %7s "`cil_u'" ///
                _col(16) %7s "`ciu_u'" _col(25) %7s "`cbl_u'" ///
                _col(34) %7s "`cbu_u'"
        }
        else if (`compact_table' & `cbands_on') {
            local cbl = el(`display_results', `j', 11)
            local cbu = el(`display_results', `j', 12)
            local cil_c : display %6.3g `cil'
            local ciu_c : display %6.3g `ciu'
            local cbl_c : display %6.3g `cbl'
            local cbu_c : display %6.3g `cbu'
            if (`cil' <= -c(maxdouble) / 2) local cil_c "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_c "inf"
            if (`cbl' <= -c(maxdouble) / 2) local cbl_c "-inf"
            if (`cbu' >= c(maxdouble) / 2) local cbu_c "inf"
            di as res %6s "`dname'" _col(8) %6.3g `estq' ///
                _col(15) %6.3g `seq' _col(22) %6s "`cil_c'" ///
                _col(29) %6s "`ciu_c'" _col(36) %6s "`cbl_c'" ///
                _col(43) %6s "`cbu_c'"
        }
        else if (`narrow_table' & `cbands_on') {
            local cbl = el(`display_results', `j', 11)
            local cbu = el(`display_results', `j', 12)
            local cil_n : display %8.4g `cil'
            local ciu_n : display %8.4g `ciu'
            local cbl_n : display %8.4g `cbl'
            local cbu_n : display %8.4g `cbu'
            if (`cil' <= -c(maxdouble) / 2) local cil_n "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_n "inf"
            if (`cbl' <= -c(maxdouble) / 2) local cbl_n "-inf"
            if (`cbu' >= c(maxdouble) / 2) local cbu_n "inf"
            di as res %8s "`dname'" _col(10) %8.4g `estq' ///
                _col(19) %8.4g `seq' _col(28) %8s "`cil_n'" ///
                _col(37) %8s "`ciu_n'" _col(46) %8s "`cbl_n'" ///
                _col(55) %8s "`cbu_n'"
        }
        else if (`cbands_on') {
            local cbl = el(`display_results', `j', 11)
            local cbu = el(`display_results', `j', 12)
            local cil_w : display %9.4g `cil'
            local ciu_w : display %9.4g `ciu'
            local cbl_w : display %9.4g `cbl'
            local cbu_w : display %9.4g `cbu'
            if (`cil' <= -c(maxdouble) / 2) local cil_w "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_w "inf"
            if (`cbl' <= -c(maxdouble) / 2) local cbl_w "-inf"
            if (`cbu' >= c(maxdouble) / 2) local cbu_w "inf"
            di as res %10s "`dname'" _col(12) %10.4g `estq' ///
                _col(23) %9.4g `seq' _col(33) %9s "`cil_w'" ///
                _col(43) %9s "`ciu_w'" _col(53) %9s "`cbl_w'" ///
                _col(63) %9s "`cbu_w'"
        }
        else if (`ultra_table') {
            local cil_u : display %7.3g `cil'
            local ciu_u : display %7.3g `ciu'
            if (`cil' <= -c(maxdouble) / 2) local cil_u "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_u "inf"
            di as res %4s "`dname'" _col(7) %8.3g `estp' ///
                _col(17) %8.3g `sep' _col(27) %8.3g `estq'
            di as res %4s "" _col(7) %7.3g `seq' ///
                _col(16) %7s "`cil_u'" _col(25) %7s "`ciu_u'"
        }
        else if (`compact_table') {
            local cil_c : display %6.3g `cil'
            local ciu_c : display %6.3g `ciu'
            if (`cil' <= -c(maxdouble) / 2) local cil_c "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_c "inf"
            di as res %6s "`dname'" _col(8) %6.3g `estp' ///
                _col(15) %6.3g `sep' ///
                _col(22) %6.3g `estq' ///
                _col(29) %6.3g `seq' ///
                _col(36) %6s "`cil_c'" ///
                _col(43) %6s "`ciu_c'"
        }
        else if (`narrow_table') {
            local cil_n : display %8.4g `cil'
            local ciu_n : display %8.4g `ciu'
            if (`cil' <= -c(maxdouble) / 2) local cil_n "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_n "inf"
            di as res %8s "`dname'" _col(10) %8.4g `estp' ///
                _col(19) %8.4g `sep' ///
                _col(28) %8.4g `estq' ///
                _col(37) %8.4g `seq' ///
                _col(46) %8s "`cil_n'" ///
                _col(55) %8s "`ciu_n'"
        }
        else {
            local cil_w : display %9.4g `cil'
            local ciu_w : display %9.4g `ciu'
            if (`cil' <= -c(maxdouble) / 2) local cil_w "-inf"
            if (`ciu' >= c(maxdouble) / 2) local ciu_w "inf"
            di as res %10s "`dname'" _col(12) %11.4g `estp' ///
                _col(24) %11.4g `sep' ///
                _col(36) %11.4g `estq' ///
                _col(48) %11.4g `seq' ///
                _col(60) %9s "`cil_w'" ///
                _col(70) %9s "`ciu_w'"
        }
    }
    di as txt "{hline `line_width'}"
    if (`cbands_on') {
        if (`ultra_table') {
            di as txt "  Uniform bands: " as res "on"
            di as txt "  Critical value: " as res %9.4g e(cb_crit)
            di as txt "  Repetitions: " as res %9.0f `repp'
        }
        else {
            di as txt "  Uniform bands: " as res "on" ///
                as txt "    Critical value: " as res %9.4g e(cb_crit) ///
                as txt "    Repetitions: " as res %9.0f `repp'
        }
    }
    else {
        di as txt "  Uniform bands: " as res "off"
    }
end

mata:
real scalar _rd2d_dist_infnorm(real rowvector x)
{
    return(max(abs(x)))
}

real colvector _rd2d_dist_rstyle_kernel(real colvector u, string scalar kernel)
{
    real colvector w

    if (kernel == "uniform") {
        w = 0.5 :* (abs(u) :<= 1)
    }
    else if (kernel == "triangular") {
        w = (1 :- abs(u)) :* (abs(u) :<= 1)
    }
    else if (kernel == "epanechnikov") {
        w = 0.75 :* (1 :- u:^2) :* (abs(u) :<= 1)
    }
    else {
        w = exp(-0.5 :* u:^2) :/ sqrt(2 * pi())
    }

    return(w)
}

real rowvector _rd2d_dist_bread_diag(real matrix bread)
{
    real matrix A, eigvec
    real rowvector eigval, pos
    real scalar k, maxeig, minpos, tol, r, cond, fallback

    A = (bread + bread') / 2
    k = cols(A)
    symeigensystem(A, eigvec, eigval)
    maxeig = max(abs(eigval))
    if (!(maxeig > 0)) return((0, ., 1))

    tol = 1e-12 * maxeig
    pos = select(eigval, eigval :> tol)
    r = length(pos)
    if (r > 0) {
        minpos = min(pos)
        cond = maxeig / minpos
    }
    else {
        cond = .
    }
    fallback = (r < k | !(cond < 1e12))

    return((r, cond, fallback))
}

real matrix _rd2d_dist_bread_inverse(real matrix bread, real rowvector diag)
{
    real matrix A

    A = (bread + bread') / 2
    if (diag[3] != 0) return(pinv(A))
    return(invsym(A))
}

real rowvector _rd2d_dist_fit_raw(
    string scalar yname,
    string scalar dname,
    string scalar tousename,
    string scalar cname,
    real scalar h,
    real scalar p,
    string scalar side,
    string scalar kernel,
    string scalar vce)
{
    real colvector y, d, C, x, w, ind, resd, hii, clusters, cidx
    real matrix X, bread, ibread, beta, meat, V, sqrtwX, Xwr, score
    real rowvector diag
    real scalar k, eN, mu, se, wvce, j, hascluster, g, factor

    st_view(y = ., ., yname, tousename)
    st_view(d = ., ., dname, tousename)
    hascluster = (cname != "")
    if (hascluster) st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)

    // Select side: control = D<0, treated = D>=0
    if (side == "control") ind = (d :< 0)
    else ind = (d :>= 0)

    y = select(y, ind)
    d = select(d, ind)
    C = select(C, ind)
    x = abs(d)
    w = _rd2d_dist_rstyle_kernel(x :/ h, kernel) :/ (h^2)
    ind = (w :> 0)
    y = select(y, ind)
    x = select(x, ind)
    w = select(w, ind)
    C = select(C, ind)

    eN = rows(y)
    k = p + 1
    if (eN < k) {
        errprintf("not enough observations for rd2d_dist fit\n")
        _error(2001)
    }

    X = J(eN, k, 1)
    for (j = 2; j <= k; j++) {
        X[, j] = x :^(j - 1)
    }

    bread = quadcross(X, w, X)
    diag = _rd2d_dist_bread_diag(bread)
    ibread = _rd2d_dist_bread_inverse(bread, diag)
    beta = ibread * quadcross(X, w, y)
    mu = beta[1,1]

    resd = y - X * beta
    sqrtwX = sqrt(w) :* X

    if (vce == "hc0") {
        wvce = 1
    }
    else if (vce == "hc1") {
        if (eN <= k) {
            errprintf("not enough effective observations for hc1 degrees-of-freedom correction\n")
            _error(2001)
        }
        wvce = (hascluster ? 1 : sqrt(eN / (eN - k)))
    }
    else if (vce == "hc2" | vce == "hc3") {
        hii = rowsum((sqrtwX * ibread) :* sqrtwX)
        if (min(1 :- hii) <= 0) {
            errprintf("leverage is too high for hc2/hc3 variance calculation\n")
            _error(2001)
        }
        if (vce == "hc2") {
            resd = resd :* sqrt(1 :/ (1 :- hii))
        }
        else {
            resd = resd :* 1 :/ (1 :- hii)
        }
    }
    else {
        errprintf("unsupported vce option\n")
        _error(198)
    }

    if (vce == "hc0" | vce == "hc1") {
        resd = resd :* wvce
    }

    if (!hascluster) {
        Xwr = X :* (w :* resd)
        meat = quadcross(Xwr, Xwr)
    }
    else {
        clusters = uniqrows(sort(C, 1))
        g = rows(clusters)
        if (g <= 1) {
            errprintf("cluster() must contain at least two clusters inside the bandwidth on the %s side\n", side)
            _error(2001)
        }
        factor = 1
        if (vce == "hc1") factor = ((eN - 1) / (eN - k)) * (g / (g - 1))
        meat = J(k, k, 0)
        for (j = 1; j <= g; j++) {
            cidx = (C :== clusters[j])
            score = quadcross(select(X, cidx), select(w, cidx) :* select(resd, cidx))
            meat = meat + score * score'
        }
        meat = meat * factor
    }
    V = ibread * meat * ibread
    se = sqrt(V[1,1])

    return((mu, se, eN, diag[1], diag[2], diag[3]))
}

real colvector _rd2d_dist_side_influence(
    real colvector yfull,
    real colvector dfull,
    real colvector Cfull,
    real scalar h,
    real scalar p,
    string scalar side,
    string scalar kernel,
    string scalar vce,
    real scalar hascluster)
{
    real colvector psi, y, d, C, x, w, active, idx, resid, hii, clusters
    real matrix X, bread, ibread, beta, sqrtwX, scores
    real rowvector diag
    real scalar nfull, n, k, j, g, factor

    nfull = rows(yfull)
    psi = J(nfull, 1, 0)
    if (side == "control") active = (dfull :< 0)
    else active = (dfull :>= 0)

    idx = selectindex(active)
    y = select(yfull, active)
    d = select(dfull, active)
    C = select(Cfull, active)
    x = abs(d)
    w = _rd2d_dist_rstyle_kernel(x :/ h, kernel) :/ (h^2)
    active = (w :> 0)
    idx = select(idx, active)
    y = select(y, active)
    C = select(C, active)
    x = select(x, active)
    w = select(w, active)

    n = rows(y)
    k = p + 1
    if (n < k) {
        errprintf("not enough observations for rd2d_dist covariance fit\n")
        _error(2001)
    }

    X = J(n, k, 1)
    for (j = 2; j <= k; j++) {
        X[, j] = x :^(j - 1)
    }

    bread = quadcross(X, w, X)
    diag = _rd2d_dist_bread_diag(bread)
    ibread = _rd2d_dist_bread_inverse(bread, diag)
    beta = ibread * quadcross(X, w, y)
    resid = y - X * beta
    sqrtwX = sqrt(w) :* X

    if (vce == "hc1") {
        if (n <= k) {
            errprintf("not enough effective observations for hc1 covariance correction\n")
            _error(2001)
        }
        if (!hascluster) resid = resid :* sqrt(n / (n - k))
    }
    else if (vce == "hc2" | vce == "hc3") {
        hii = rowsum((sqrtwX * ibread) :* sqrtwX)
        if (min(1 :- hii) <= 0) {
            errprintf("leverage is too high for hc2/hc3 covariance calculation\n")
            _error(2001)
        }
        if (vce == "hc2") {
            resid = resid :* sqrt(1 :/ (1 :- hii))
        }
        else {
            resid = resid :* 1 :/ (1 :- hii)
        }
    }

    scores = (X :* (w :* resid)) * ibread'
    if (hascluster) {
        clusters = uniqrows(sort(C, 1))
        g = rows(clusters)
        if (g <= 1) {
            errprintf("cluster() must contain at least two clusters inside the bandwidth on the %s side\n", side)
            _error(2001)
        }
        factor = 1
        if (vce == "hc1") factor = ((n - 1) / (n - k)) * (g / (g - 1))
        scores = scores :* sqrt(factor)
    }
    psi[idx] = scores[, 1]
    return(psi)
}

real matrix _rd2d_dist_cluster_crossprod(real matrix psi, real colvector C)
{
    real colvector clusters, cidx
    real matrix cps
    real scalar j, g

    clusters = uniqrows(sort(C, 1))
    g = rows(clusters)
    if (g <= 1) {
        errprintf("cluster() must contain at least two clusters\n")
        _error(2001)
    }
    cps = J(g, cols(psi), 0)
    for (j = 1; j <= g; j++) {
        cidx = (C :== clusters[j])
        cps[j,] = colsum(select(psi, cidx))
    }
    return(quadcross(cps, cps))
}

void _rd2d_dist_cov_q(
    string scalar yname,
    string scalar dnames,
    string scalar tousename,
    string scalar hqname,
    real scalar p,
    string scalar kernel,
    string scalar vce,
    string scalar cname,
    string scalar covname,
    string scalar corrname)
{
    real colvector y, d, C, se
    real matrix D, hq, psit, psic, psi, cov, corr
    real scalar j, l, k, denom, hascluster

    st_view(y = ., ., yname, tousename)
    st_view(D = ., ., tokens(dnames), tousename)
    hascluster = (cname != "")
    if (hascluster) st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)
    hq = st_matrix(hqname)
    k = cols(D)
    psit = J(rows(y), k, 0)
    psic = J(rows(y), k, 0)

    for (j = 1; j <= k; j++) {
        d = D[, j]
        psit[, j] = _rd2d_dist_side_influence(y, d, C, hq[j, 2], p, "treated", kernel, vce, hascluster)
        psic[, j] = _rd2d_dist_side_influence(y, d, C, hq[j, 1], p, "control", kernel, vce, hascluster)
    }

    psi = psit - psic
    if (!hascluster) {
        cov = quadcross(psi, psi)
    }
    else {
        cov = _rd2d_dist_cluster_crossprod(psi, C)
    }
    se = sqrt(diagonal(cov))
    corr = J(k, k, .)
    for (j = 1; j <= k; j++) {
        for (l = 1; l <= k; l++) {
            denom = se[j] * se[l]
            if (denom > 0 & denom < .) corr[j, l] = cov[j, l] / denom
        }
    }
    corr = (corr + corr') / 2
    for (j = 1; j <= k; j++) {
        if (se[j] > 0 & se[j] < .) corr[j, j] = 1
    }

    st_matrix(covname, cov)
    st_matrix(corrname, corr)
}

real scalar _rd2d_dist_cb_quantile(real colvector x, real scalar q)
{
    real scalar n, h, j
    real colvector sx

    sx = sort(x, 1)
    n = rows(sx)
    if (n == 0) return(.)
    if (n == 1) return(sx[1])
    if (q <= 0) return(sx[1])
    if (q >= 1) return(sx[n])

    h = n * q
    j = floor(h)
    if (abs(h - j) < 1e-12) {
        if (j < 1) return(sx[1])
        if (j >= n) return(sx[n])
        return((sx[j] + sx[j + 1]) / 2)
    }

    return(sx[ceil(h)])
}

void _rd2d_dist_cb_from_corr(
    string scalar corrname,
    real scalar rep,
    real scalar level,
    string scalar side,
    string scalar critname,
    string scalar adjustedname,
    string scalar mineigname)
{
    real matrix corr, eigvec, L, sim
    real rowvector eigval
    real colvector tvec
    real scalar k, adjusted, mineig, tol

    corr = st_matrix(corrname)
    corr = (corr + corr') / 2
    k = rows(corr)
    if (k == 0 | k != cols(corr) | any(corr :>= .) | any(diagonal(corr) :<= 0)) {
        errprintf("invalid correlation matrix for cbands simulation\n")
        _error(498)
    }
    symeigensystem(corr, eigvec, eigval)
    mineig = min(eigval)
    tol = 1e-12 * max((1, max(abs(eigval))))
    adjusted = any(eigval :< -tol)
    eigval = eigval :* (eigval :> tol)
    L = eigvec * diag(sqrt(eigval'))
    sim = invnormal(runiform(rep, k)) * L'

    if (side == "two") {
        tvec = rowmax(abs(sim))
    }
    else {
        tvec = rowmax(sim)
    }

    st_numscalar(critname, _rd2d_dist_cb_quantile(tvec, level / 100))
    st_numscalar(adjustedname, adjusted)
    st_numscalar(mineigname, mineig)
}
end
