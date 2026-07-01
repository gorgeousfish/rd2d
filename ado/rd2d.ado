*! version 1.1.0 24jun2026
program define rd2d, eclass
    version 16.0
    if "`1'" == "version" | "`1'" == ",version" {
        di as txt "rd2d version 1.1.0 (24 June 2026)"
        exit
    }
    syntax varlist(min=4 max=4 numeric) [if] [in], AT(string asis) ///
        [H(string asis) DERiv(string asis) TANGvec(string asis) P(integer 1) ///
        Q(string) KERnel(string) KTYPE(string) VCE(string) BWCHeck(string) ///
         MASSPoints(string) SCALEregul(real 3) SCALEbiascrct(real 1) STDVars RAWVARS ///
         LEVEL(real 95) REPP(integer 1000) SIDE(string) BWSELect(string) METHOD(string) RBC(string) ///
         CBANDs NOCBANDs CLuster(varname) ///
         MP(string) SR(real -1) BWS(string) KT(string)]

    // --- cbands default logic (v1.2.0: default ON, nocbands to disable) ---
    local cbands_on = ("`nocbands'" == "")

    gettoken yvar rest : varlist
    gettoken x1var rest : rest
    gettoken x2var tvar : rest
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
        if (`scaleregul' != 3) {
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
    if ("`kt'" != "") {
        if ("`ktype'" != "") {
            di as err "cannot specify both kt() and ktype()"
            exit 198
        }
        local ktype "`kt'"
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
        di as err "scaleregul() must be a non-negative number (default: 3 for location); got `scaleregul'"
        exit 198
    }
    if (`scalebiascrct' >= . | `scalebiascrct' < 0) {
        di as err "scalebiascrct() must be finite and nonnegative"
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

    local ktype = lower("`ktype'")
    if ("`ktype'" == "") local ktype "prod"
    if !inlist("`ktype'", "prod", "rad") {
        di as err "ktype() must be {bf:prod} (product kernel) or {bf:rad} (radial kernel)"
        exit 198
    }

    local bwselect = lower("`bwselect'")
    if ("`bwselect'" == "") local bwselect "mserd"
    if !inlist("`bwselect'", "mserd", "imserd", "msetwo", "imsetwo") {
        di as err "bwselect() must be one of: {bf:mserd}, {bf:imserd}, {bf:msetwo}, {bf:imsetwo}"
        di as err "  mserd = MSE-optimal common bandwidth; msetwo = side-specific bandwidths"
        exit 198
    }

    local method = lower("`method'")
    if ("`method'" == "") local method "dpi"
    if !inlist("`method'", "dpi", "rot") {
        di as err "method() must be {bf:dpi} (data-driven plug-in) or {bf:rot} (rule of thumb)"
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

    // --- stdvars / rawvars logic (4.1) ---
    if ("`stdvars'" != "" & "`rawvars'" != "") {
        di as err "cannot specify both stdvars and rawvars"
        exit 198
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

    local bwcheck_raw "`bwcheck'"
    if ("`bwcheck_raw'" == "") {
        local bwcheck = 50 + `p' + 1
    }
    else {
        capture confirm integer number `bwcheck_raw'
        if (_rc | `bwcheck_raw' < 0) {
            di as err "bwcheck() must be a nonnegative integer; got '`bwcheck_raw''"
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
    quietly count if `touse' & !inlist(`tvar', 0, 1)
    if (r(N) > 0) {
        di as err "treatment indicator must contain only 0 and 1"
        exit 198
    }
    if (`bwcheck' > 0 & `N' < `bwcheck') {
        di as err "Insufficient observations (N=`N') for bandwidth estimation. Minimum required: `bwcheck' (bwcheck option). Consider reducing bwcheck() or providing manual h()."
        exit 2001
    }

    local atlist "`at'"
    local atcount : word count `atlist'
    if (`atcount' < 2 | mod(`atcount', 2) != 0) {
        di as err "at() requires pairs of coordinates: at(b1 b2) or at(b1 b2 b1 b2 ...)"
        di as err "  Example: at(0 0) for single point, at(0 0 1 1) for two points"
        exit 198
    }
    forvalues i = 1/`atcount' {
        local av : word `i' of `atlist'
        capture confirm number `av'
        if (_rc) {
            di as err "at() must contain numeric b1 b2 evaluation-point pairs"
            exit 198
        }
        if (`av' >= .) {
            di as err "at() must contain finite b1 b2 evaluation-point pairs"
            exit 198
        }
    }
    local neval = `atcount' / 2

    local d1 0
    local d2 0
    if ("`deriv'" != "" & "`tangvec'" == "") {
        local derivlist "`deriv'"
        local derivcount : word count `derivlist'
        if (`derivcount' != 2) {
            di as err "deriv() must contain exactly two nonnegative integers"
            exit 198
        }
        local d1 : word 1 of `derivlist'
        local d2 : word 2 of `derivlist'
        capture confirm integer number `d1'
        if (_rc) {
            di as err "deriv() must contain exactly two nonnegative integers"
            exit 198
        }
        capture confirm integer number `d2'
        if (_rc) {
            di as err "deriv() must contain exactly two nonnegative integers"
            exit 198
        }
        if (`d1' < 0 | `d2' < 0 | `d1' + `d2' > `p') {
            di as err "deriv() components sum (d1+d2=`=`d1'+`d2'') exceeds polynomial order p(`p'); reduce deriv() or increase p()"
            exit 198
        }
    }

    local has_tangvec = 0
    local tv1 .
    local tv2 .
    if ("`tangvec'" != "") {
        if (`neval' > 1) {
            di as err "tangvec() is currently supported only with one at() point"
            exit 198
        }
        local tvlist "`tangvec'"
        local tvcount : word count `tvlist'
        if (`tvcount' != 2) {
            di as err "tangvec() must contain exactly two numbers for the single at() point"
            exit 198
        }
        local tv1 : word 1 of `tvlist'
        local tv2 : word 2 of `tvlist'
        capture confirm number `tv1'
        if (_rc) {
            di as err "tangvec() must contain exactly two numbers"
            exit 198
        }
        capture confirm number `tv2'
        if (_rc) {
            di as err "tangvec() must contain exactly two numbers"
            exit 198
        }
        if (`p' < 1) {
            di as err "tangvec() requires p() of at least 1"
            exit 198
        }
        if (`tv1' >= . | `tv2' >= . | (`tv1' == 0 & `tv2' == 0)) {
            di as err "tangvec() must contain a nonzero finite direction vector"
            exit 198
        }
        local has_tangvec = 1
        if ("`deriv'" != "") di as txt "warning: tangvec() provided; deriv() is ignored."
    }
    local derivsum = `d1' + `d2'
    if ("`tangvec'" != "") local derivsum = 1

    if ("`rbc'" == "off") local qeff = `p'
    else if (`qeff' < 0) local qeff = `p' + 1
    if (`qeff' < `p') {
        di as err "q() must be an integer >= p (currently p=`p'); got q=`qeff'"
        di as err "  q=p+1 is the default for bias correction"
        exit 198
    }

    local qbasis = (`qeff' + 1) * (`qeff' + 2) / 2
    quietly count if `touse' & `tvar' == 0
    local N0 = r(N)
    quietly count if `touse' & `tvar' == 1
    local N1 = r(N)
    if (`N0' < `qbasis' | `N1' < `qbasis') {
        local failside = cond(`N0' < `qbasis', "control", "treated")
        local failnh = cond(`N0' < `qbasis', `N0', `N1')
        di as err "Insufficient observations on `failside' side (Nh=`failnh', need `qbasis')"
        di as err "  Suggestion: increase bandwidth with h() or reduce polynomial order p()"
        exit 2001
    }

    tempname bws massinfo bwconsts bwbounds bwsout results resultsA0 resultsA1 hqmat
    if ("`h'" != "") {
        local hlist "`h'"
        local hcount : word count `hlist'
        if !inlist(`hcount', 1, 4, 4 * `neval') {
            di as err "h() requires 1, 4, or 4*neval values; got `hcount' for `neval' evaluation points"
            exit 198
        }
        forvalues i = 1/`hcount' {
            local hv : word `i' of `hlist'
            capture confirm number `hv'
            if (_rc) {
                di as err "h() bandwidths must be positive numbers"
                exit 198
            }
            if (`hv' >= . | `hv' <= 0) {
                di as err "h() bandwidths must be finite positive numbers"
                exit 198
            }
        }
        matrix `bws' = J(`neval', 6, .)
        forvalues j = 1/`neval' {
            local ai = 2 * `j' - 1
            local bi = 2 * `j'
            local b1 : word `ai' of `atlist'
            local b2 : word `bi' of `atlist'
            matrix `bws'[`j', 1] = `b1'
            matrix `bws'[`j', 2] = `b2'
            if (`hcount' == 1) {
                local hv : word 1 of `hlist'
                matrix `bws'[`j', 3] = `hv'
                matrix `bws'[`j', 4] = `hv'
                matrix `bws'[`j', 5] = `hv'
                matrix `bws'[`j', 6] = `hv'
            }
            else if (`hcount' == 4) {
                forvalues i = 1/4 {
                    local hv : word `i' of `hlist'
                    matrix `bws'[`j', `i' + 2] = `hv'
                }
            }
            else {
                forvalues i = 1/4 {
                    local hi = 4 * (`j' - 1) + `i'
                    local hv : word `hi' of `hlist'
                    matrix `bws'[`j', `i' + 2] = `hv'
                }
            }
        }
        // massinfo matrix: 4 columns = [M, M0, M1, mass_ratio]
        //   M  = number of unique (x1, x2) points overall
        //   M0 = unique points on control side (t==0)
        //   M1 = unique points on treated side (t==1)
        //   mass_ratio = 1 - M/N; proportion of duplicated locations
        // When h() is provided, mass-point detection uses the full sample
        // (not the bandwidth-local sample) because bandwidths are user-supplied.
        matrix `massinfo' = J(`neval', 4, .)
        if ("`masspoints'" != "off") {
            tempvar mass_all mass_0 mass_1
            quietly egen long `mass_all' = group(`x1var' `x2var') if `touse'
            quietly summarize `mass_all' if `touse', meanonly
            local M = r(max)
            quietly egen long `mass_0' = group(`x1var' `x2var') if `touse' & `tvar' == 0
            quietly summarize `mass_0' if `touse' & `tvar' == 0, meanonly
            local M0 = r(max)
            quietly egen long `mass_1' = group(`x1var' `x2var') if `touse' & `tvar' == 1
            quietly summarize `mass_1' if `touse' & `tvar' == 1, meanonly
            local M1 = r(max)
            local mass = 1 - `M' / `N'
            forvalues j = 1/`neval' {
                matrix `massinfo'[`j', 1] = `M'
                matrix `massinfo'[`j', 2] = `M0'
                matrix `massinfo'[`j', 3] = `M1'
                matrix `massinfo'[`j', 4] = `mass'
            }
            _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')
        }
        else {
            // masspoints="off": store N/N0/N1 as if every observation is unique
            // (M=N implies mass_ratio=0), so downstream code has a uniform schema.
            forvalues j = 1/`neval' {
                matrix `massinfo'[`j', 1] = `N'
                matrix `massinfo'[`j', 2] = `N0'
                matrix `massinfo'[`j', 3] = `N1'
                matrix `massinfo'[`j', 4] = 0
            }
        }
        local bwselect "user provided"
        local bwsource "user"
        // stdvars in manual h() path: only on if user explicitly says stdvars
        local stdvars_active = ("`stdvars'" != "")
        local stdvars_source = cond("`stdvars'" != "" | "`rawvars'" != "", "user", "")
    }
    else {
        // stdvars default: ON for automatic bandwidth unless rawvars specified
        local stdvars_active = 1
        if ("`rawvars'" != "") local stdvars_active = 0
        if ("`stdvars'" != "") local stdvars_active = 1
        local stdopt ""
        if (`stdvars_active') local stdopt "stdvars"
        local stdflag = `stdvars_active'
        local clusteropt ""
        if ("`clusterwork'" != "") local clusteropt "cluster(`clusterwork')"
        local targetopt "deriv(`d1' `d2')"
        if ("`tangvec'" != "") local targetopt "tangvec(`tangvec')"
        matrix `bws' = J(`neval', 6, .)
        matrix `massinfo' = J(`neval', 4, .)
        matrix `bwconsts' = J(`neval', 10, .)
        matrix `bwbounds' = J(`neval', 4, .)
        // Initialize DPI multi-point cache (invalidated; populated on first rdbw2d call)
        mata: _rdbw2d_cache_valid = 0
        forvalues j = 1/`neval' {
            local ai = 2 * `j' - 1
            local bi = 2 * `j'
            local b1 : word `ai' of `atlist'
            local b2 : word `bi' of `atlist'
            quietly rdbw2d `yvar' `x1var' `x2var' `tvar' if `touse', at(`b1' `b2') ///
                p(`p') kernel(`kernel') ktype(`ktype') bwselect(`bwselect') ///
                method(`method') vce(`vce') bwcheck(`bwcheck') masspoints(`masspoints') ///
                scaleregul(`scaleregul') scalebiascrct(`scalebiascrct') `stdopt' `clusteropt' ///
                `targetopt'
            tempname bwone massone mseone boundone
            matrix `bwone' = r(bws)
            matrix `massone' = r(masspoints)
            matrix `mseone' = r(mseconsts)
            mata: _rd2d_loc_bwcheck_bounds("`x1var'", "`x2var'", "`tvar'", "`touse'", ///
                `b1', `b2', "`ktype'", `bwcheck', "`masspoints'", `stdflag', "`boundone'")
            forvalues c = 1/6 {
                matrix `bws'[`j', `c'] = `bwone'[1, `c']
            }
            forvalues c = 1/4 {
                matrix `massinfo'[`j', `c'] = `massone'[1, `c']
            }
            forvalues c = 1/10 {
                matrix `bwconsts'[`j', `c'] = `mseone'[1, `c']
            }
            forvalues c = 1/4 {
                matrix `bwbounds'[`j', `c'] = `boundone'[1, `c']
            }
        }
        // Clear DPI multi-point cache (release memory)
        mata: _rdbw2d_cache_valid = 0
        mata: _rdbw2d_cache_x1_0 = _rdbw2d_cache_x2_0 = _rdbw2d_cache_y_0 = _rdbw2d_cache_C_0 = J(0, 1, .)
        mata: _rdbw2d_cache_x1_1 = _rdbw2d_cache_x2_1 = _rdbw2d_cache_y_1 = _rdbw2d_cache_C_1 = J(0, 1, .)
        mata: _rdbw2d_cache_unique0 = _rdbw2d_cache_unique1 = J(0, 0, .)
        // Mass-point warning for automatic-bandwidth path: rdbw2d fills
        // massinfo per evaluation point from its local sample.
        _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')
        if inlist("`bwselect'", "imserd", "imsetwo") {
            tempname sd1 sd2 Vbar Bbar Vbar0 Vbar1 Bbar0 Bbar1 himse himse0 himse1 den
            scalar `sd1' = 1
            scalar `sd2' = 1
            if (`stdvars_active') {
                quietly summarize `x1var' if `touse'
                scalar `sd1' = r(sd)
                quietly summarize `x2var' if `touse'
                scalar `sd2' = r(sd)
            }
            if ("`bwselect'" == "imserd") {
                scalar `Vbar' = 0
                scalar `Bbar' = 0
                forvalues j = 1/`neval' {
                    scalar `Vbar' = `Vbar' + (`bwconsts'[`j', 5] + `bwconsts'[`j', 6]) / `neval'
                    scalar `Bbar' = `Bbar' + ((`bwconsts'[`j', 3] + `scalebiascrct' * `bwconsts'[`j', 7] - ///
                        `bwconsts'[`j', 4] - `scalebiascrct' * `bwconsts'[`j', 8])^2 + ///
                        `scaleregul' * `bwconsts'[`j', 9] + `scaleregul' * `bwconsts'[`j', 10]) / `neval'
                }
                scalar `den' = (2 * `p' + 2 - 2 * `derivsum') * `Bbar'
                if (!(`Vbar' > 0) | !(`den' > 0) | missing(`Vbar') | missing(`den')) {
                    di as err "integrated bandwidth constants could not be calculated"
                    exit 498
                }
                scalar `himse' = ((2 + 2 * `derivsum') * `Vbar' / `den')^(1 / (2 * `p' + 4))
                if (!(`himse' > 0) | missing(`himse')) {
                    di as err "integrated bandwidth could not be calculated"
                    exit 498
                }
                tempname hj hlo hhi
                forvalues j = 1/`neval' {
                    scalar `hj' = `himse'
                    if (`bwcheck' > 0) {
                        scalar `hlo' = max(`bwbounds'[`j', 1], `bwbounds'[`j', 2])
                        scalar `hhi' = max(`bwbounds'[`j', 3], `bwbounds'[`j', 4])
                        scalar `hj' = min(max(`hj', `hlo'), `hhi')
                    }
                    matrix `bws'[`j', 3] = `hj' * `sd1'
                    matrix `bws'[`j', 4] = `hj' * `sd2'
                    matrix `bws'[`j', 5] = `hj' * `sd1'
                    matrix `bws'[`j', 6] = `hj' * `sd2'
                }
            }
            else {
                scalar `Vbar0' = 0
                scalar `Vbar1' = 0
                scalar `Bbar0' = 0
                scalar `Bbar1' = 0
                forvalues j = 1/`neval' {
                    scalar `Vbar0' = `Vbar0' + `bwconsts'[`j', 5] / `neval'
                    scalar `Vbar1' = `Vbar1' + `bwconsts'[`j', 6] / `neval'
                    scalar `Bbar0' = `Bbar0' + ((`bwconsts'[`j', 3] + `scalebiascrct' * `bwconsts'[`j', 7])^2 + ///
                        `scaleregul' * `bwconsts'[`j', 9]) / `neval'
                    scalar `Bbar1' = `Bbar1' + ((`bwconsts'[`j', 4] + `scalebiascrct' * `bwconsts'[`j', 8])^2 + ///
                        `scaleregul' * `bwconsts'[`j', 10]) / `neval'
                }
                scalar `den' = (2 * `p' + 2 - 2 * `derivsum') * `Bbar0'
                if (!(`Vbar0' > 0) | !(`den' > 0) | missing(`Vbar0') | missing(`den')) {
                    di as err "integrated control bandwidth constants could not be calculated"
                    exit 498
                }
                scalar `himse0' = ((2 + 2 * `derivsum') * `Vbar0' / `den')^(1 / (2 * `p' + 4))
                scalar `den' = (2 * `p' + 2 - 2 * `derivsum') * `Bbar1'
                if (!(`Vbar1' > 0) | !(`den' > 0) | missing(`Vbar1') | missing(`den')) {
                    di as err "integrated treated bandwidth constants could not be calculated"
                    exit 498
                }
                scalar `himse1' = ((2 + 2 * `derivsum') * `Vbar1' / `den')^(1 / (2 * `p' + 4))
                if (!(`himse0' > 0) | !(`himse1' > 0) | missing(`himse0') | missing(`himse1')) {
                    di as err "integrated bandwidth could not be calculated"
                    exit 498
                }
                tempname h0j h1j
                forvalues j = 1/`neval' {
                    scalar `h0j' = `himse0'
                    scalar `h1j' = `himse1'
                    if (`bwcheck' > 0) {
                        scalar `h0j' = min(max(`h0j', `bwbounds'[`j', 1]), `bwbounds'[`j', 3])
                        scalar `h1j' = min(max(`h1j', `bwbounds'[`j', 2]), `bwbounds'[`j', 4])
                    }
                    matrix `bws'[`j', 3] = `h0j' * `sd1'
                    matrix `bws'[`j', 4] = `h0j' * `sd2'
                    matrix `bws'[`j', 5] = `h1j' * `sd1'
                    matrix `bws'[`j', 6] = `h1j' * `sd2'
                }
            }
        }
        local bwsource "automatic"
    }

    // --- Bandwidth zero/negative check (4.3) ---
    forvalues j = 1/`neval' {
        if (`bws'[`j', 3] <= 0 | `bws'[`j', 4] <= 0 | `bws'[`j', 5] <= 0 | `bws'[`j', 6] <= 0) {
            di as err "computed bandwidth is zero or negative at point j=`j'; check data density near this evaluation point"
            exit 498
        }
    }

    tempname targetp targetq
    local fitcluster ""
    if ("`clusterwork'" != "") local fitcluster "cluster(`clusterwork')"
    quietly _rd2d_build_target, p(`p') deriv(`d1' `d2') tangvec(`tangvec')
    matrix `targetp' = r(target)
    quietly _rd2d_build_target, p(`qeff') deriv(`d1' `d2') tangvec(`tangvec')
    matrix `targetq' = r(target)

    tempname covp corrp covq corrq cbcrit cbpsd cbmineig diagnostics
    matrix `results' = J(`neval', 18, .)
    matrix `resultsA0' = J(`neval', 9, .)
    matrix `resultsA1' = J(`neval', 9, .)
    matrix `bwsout' = J(`neval', 8, .)
    matrix `hqmat' = J(`neval', 4, .)
    matrix `diagnostics' = J(`neval', 14, .)
    matrix `covq' = J(`neval', `neval', 0)
    matrix `corrq' = J(`neval', `neval', .)
    scalar `cbcrit' = .
    scalar `cbpsd' = 0
    scalar `cbmineig' = .
    local rownames ""
    local anyfallback 0

    forvalues j = 1/`neval' {
        if (`j' == 1) local rownames "at`j'"
        else local rownames "`rownames' at`j'"
        local b1 = `bws'[`j', 1]
        local b2 = `bws'[`j', 2]
        tempname b0p V0p b1p V1p b0q V0q b1q V1q
        tempname mu0p mu1p se0p se1p mu0q mu1q se0q se1q
        tempname taup sep tauq seq zval pvalue crit cil ciu
        tempname h01 h02 h11 h12 nh0 nh1 tmp
        tempname r0p c0p f0p r1p c1p f1p r0q c0q f0q r1q c1q f1q

        scalar `h01' = `bws'[`j', 3]
        scalar `h02' = `bws'[`j', 4]
        scalar `h11' = `bws'[`j', 5]
        scalar `h12' = `bws'[`j', 6]

        _rd2d_loc_fit `yvar' `x1var' `x2var' `tvar' if `touse', at(`b1' `b2') ///
            hx(`=`h01'') hy(`=`h02'') p(`p') side(control) kernel(`kernel') ///
            ktype(`ktype') vce(`vce') `fitcluster'
        matrix `b0p' = r(b)
        matrix `V0p' = r(V)
        scalar `nh0' = round(r(N_h))
        scalar `r0p' = r(rank)
        scalar `c0p' = r(condition)
        scalar `f0p' = r(fallback)

        _rd2d_loc_fit `yvar' `x1var' `x2var' `tvar' if `touse', at(`b1' `b2') ///
            hx(`=`h11'') hy(`=`h12'') p(`p') side(treated) kernel(`kernel') ///
            ktype(`ktype') vce(`vce') `fitcluster'
        matrix `b1p' = r(b)
        matrix `V1p' = r(V)
        scalar `nh1' = round(r(N_h))
        scalar `r1p' = r(rank)
        scalar `c1p' = r(condition)
        scalar `f1p' = r(fallback)

        _rd2d_loc_fit `yvar' `x1var' `x2var' `tvar' if `touse', at(`b1' `b2') ///
            hx(`=`h01'') hy(`=`h02'') p(`qeff') side(control) kernel(`kernel') ///
            ktype(`ktype') vce(`vce') `fitcluster'
        matrix `b0q' = r(b)
        matrix `V0q' = r(V)
        scalar `r0q' = r(rank)
        scalar `c0q' = r(condition)
        scalar `f0q' = r(fallback)

        _rd2d_loc_fit `yvar' `x1var' `x2var' `tvar' if `touse', at(`b1' `b2') ///
            hx(`=`h11'') hy(`=`h12'') p(`qeff') side(treated) kernel(`kernel') ///
            ktype(`ktype') vce(`vce') `fitcluster'
        matrix `b1q' = r(b)
        matrix `V1q' = r(V)
        scalar `r1q' = r(rank)
        scalar `c1q' = r(condition)
        scalar `f1q' = r(fallback)

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

        matrix `tmp' = `targetp' * `b0p''
        scalar `mu0p' = `tmp'[1,1]
        matrix `tmp' = `targetp' * `b1p''
        scalar `mu1p' = `tmp'[1,1]
        matrix `tmp' = `targetp' * `V0p' * `targetp''
        scalar `se0p' = sqrt(max(0, `tmp'[1,1]))
        matrix `tmp' = `targetp' * `V1p' * `targetp''
        scalar `se1p' = sqrt(max(0, `tmp'[1,1]))

        matrix `tmp' = `targetq' * `b0q''
        scalar `mu0q' = `tmp'[1,1]
        matrix `tmp' = `targetq' * `b1q''
        scalar `mu1q' = `tmp'[1,1]
        matrix `tmp' = `targetq' * `V0q' * `targetq''
        scalar `se0q' = sqrt(max(0, `tmp'[1,1]))
        matrix `tmp' = `targetq' * `V1q' * `targetq''
        scalar `se1q' = sqrt(max(0, `tmp'[1,1]))

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

        matrix `results'[`j', 1] = `b1'
        matrix `results'[`j', 2] = `b2'
        matrix `results'[`j', 3] = `taup'
        matrix `results'[`j', 4] = `sep'
        matrix `results'[`j', 5] = `tauq'
        matrix `results'[`j', 6] = `seq'
        matrix `results'[`j', 7] = `zval'
        matrix `results'[`j', 8] = `pvalue'
        matrix `results'[`j', 9] = `cil'
        matrix `results'[`j', 10] = `ciu'
        matrix `results'[`j', 13] = `h01'
        matrix `results'[`j', 14] = `h02'
        matrix `results'[`j', 15] = `h11'
        matrix `results'[`j', 16] = `h12'
        matrix `results'[`j', 17] = `nh0'
        matrix `results'[`j', 18] = `nh1'

        matrix `resultsA0'[`j', 1] = `b1'
        matrix `resultsA0'[`j', 2] = `b2'
        matrix `resultsA0'[`j', 3] = `mu0p'
        matrix `resultsA0'[`j', 4] = `se0p'
        matrix `resultsA0'[`j', 5] = `mu0q'
        matrix `resultsA0'[`j', 6] = `se0q'
        matrix `resultsA0'[`j', 7] = `h01'
        matrix `resultsA0'[`j', 8] = `h02'
        matrix `resultsA0'[`j', 9] = `nh0'

        matrix `resultsA1'[`j', 1] = `b1'
        matrix `resultsA1'[`j', 2] = `b2'
        matrix `resultsA1'[`j', 3] = `mu1p'
        matrix `resultsA1'[`j', 4] = `se1p'
        matrix `resultsA1'[`j', 5] = `mu1q'
        matrix `resultsA1'[`j', 6] = `se1q'
        matrix `resultsA1'[`j', 7] = `h11'
        matrix `resultsA1'[`j', 8] = `h12'
        matrix `resultsA1'[`j', 9] = `nh1'

        matrix `bwsout'[`j', 1] = `b1'
        matrix `bwsout'[`j', 2] = `b2'
        matrix `bwsout'[`j', 3] = `h01'
        matrix `bwsout'[`j', 4] = `h02'
        matrix `bwsout'[`j', 5] = `h11'
        matrix `bwsout'[`j', 6] = `h12'
        matrix `bwsout'[`j', 7] = `nh0'
        matrix `bwsout'[`j', 8] = `nh1'
        matrix `hqmat'[`j', 1] = `h01'
        matrix `hqmat'[`j', 2] = `h02'
        matrix `hqmat'[`j', 3] = `h11'
        matrix `hqmat'[`j', 4] = `h12'
        matrix `covq'[`j', `j'] = `seq'^2
        matrix `corrq'[`j', `j'] = cond(`seq' > 0 & `seq' < ., 1, .)
        matrix `diagnostics'[`j', 1] = `b1'
        matrix `diagnostics'[`j', 2] = `b2'
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

    matrix colnames `results' = b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h01 h02 h11 h12 Nh0 Nh1
    matrix colnames `resultsA0' = b1 b2 Est_p Se_p Est_q Se_q h01 h02 Nh0
    matrix colnames `resultsA1' = b1 b2 Est_p Se_p Est_q Se_q h11 h12 Nh1
    matrix colnames `bwsout' = b1 b2 h01 h02 h11 h12 Nh0 Nh1
    matrix colnames `diagnostics' = b1 b2 rank_p0 cond_p0 fb_p0 rank_p1 cond_p1 fb_p1 rank_q0 cond_q0 fb_q0 rank_q1 cond_q1 fb_q1
    matrix colnames `massinfo' = M M0 M1 mass
    matrix rownames `results' = `rownames'
    matrix rownames `resultsA0' = `rownames'
    matrix rownames `resultsA1' = `rownames'
    matrix rownames `bwsout' = `rownames'
    matrix rownames `diagnostics' = `rownames'
    matrix rownames `massinfo' = `rownames'
    matrix colnames `covq' = `rownames'
    matrix rownames `covq' = `rownames'
    matrix colnames `corrq' = `rownames'
    matrix rownames `corrq' = `rownames'

    if ("`clusterwork'" != "") {
        mata: _rd2d_loc_cov_q("`yvar'", "`x1var'", "`x2var'", "`tvar'", "`touse'", ///
            "`bws'", "`hqmat'", "`targetp'", `p', "`kernel'", "`ktype'", "`vce'", ///
            "0", "`clusterwork'", "`covp'", "`corrp'")
        forvalues j = 1/`neval' {
            matrix `results'[`j', 4] = sqrt(`covp'[`j', `j'])
        }
    }
    if (`neval' > 1 | `cbands_on' | "`clusterwork'" != "") {
        mata: _rd2d_loc_cov_q("`yvar'", "`x1var'", "`x2var'", "`tvar'", "`touse'", ///
            "`bws'", "`hqmat'", "`targetq'", `qeff', "`kernel'", "`ktype'", "`vce'", ///
            "0", "`clusterwork'", "`covq'", "`corrq'")
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
        if (`qeff' == `p' & `results'[`j', 13] == `results'[`j', 15] & ///
            `results'[`j', 14] == `results'[`j', 16]) {
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
        mata: _rd2d_loc_cb_from_corr("`corrq'", `repp', `level', "`side'", ///
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
    matrix colnames `covq' = `rownames'
    matrix rownames `covq' = `rownames'
    matrix colnames `corrq' = `rownames'
    matrix rownames `corrq' = `rownames'

    tempname eb eV
    matrix `eb' = J(1, `neval', .)
    forvalues j = 1/`neval' {
        matrix `eb'[1, `j'] = `results'[`j', 5]
    }
    matrix colnames `eb' = `rownames'
    matrix `eV' = `covq'

    ereturn clear
    ereturn post `eb' `eV'
    matrix colnames `results' = b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h01 h02 h11 h12 Nh0 Nh1
    matrix colnames `resultsA0' = b1 b2 Est_p Se_p Est_q Se_q h01 h02 Nh0
    matrix colnames `resultsA1' = b1 b2 Est_p Se_p Est_q Se_q h11 h12 Nh1
    matrix colnames `bwsout' = b1 b2 h01 h02 h11 h12 Nh0 Nh1
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
    ereturn scalar scalebiascrct = `scalebiascrct'
    ereturn scalar derivsum = `derivsum'
    ereturn local cmd "rd2d"
    ereturn local kernel "`kernel'"
    ereturn local ktype "`ktype'"
    ereturn local bwselect "`bwselect'"
    ereturn local method "`method'"
    ereturn local vce "`vce'"
    ereturn local rbc "`rbc'"
    ereturn local side "`side'"
    ereturn local masspoints_opt "`masspoints'"
    ereturn local bwsource "`bwsource'"
    ereturn local cbands = cond(`cbands_on', "on", "off")
    ereturn local clustered = cond("`clustername'" != "", "on", "off")
    ereturn local cluster "`clustername'"
    ereturn local stdvars = cond("`bwsource'" == "automatic" & `stdvars_active', "on", "off")
    ereturn scalar stdvars_flag = cond("`bwsource'" == "automatic" & `stdvars_active', 1, 0)
    ereturn local deriv "`d1' `d2'"
    ereturn local tangvec "`tangvec'"
    ereturn local fallback = cond(`anyfallback', "pinv", "invsym")
    ereturn local version "1.1.0"
    ereturn local depvar "`yvar'"
    ereturn local x1var "`x1var'"
    ereturn local x2var "`x2var'"
    ereturn local tvar "`tvar'"
    ereturn local atstring `"`at'"'

    // --- cbands multi-point hint (4.2) ---
    if (`neval' > 1 & !`cbands_on') {
        di as txt "note: consider removing {cmd:nocbands} option for uniform inference over `neval' evaluation points"
    }


    // ===================================================================
    // TABLE DISPLAY SECTION
    // Layout modes: ultra (<50), compact (<63), narrow (<79), normal (>=79)
    // With cbands: shows Est.q, SE.q, CI.lo, CI.hi, CB.lo, CB.hi
    // Without cbands: shows Est.p, SE.p, Est.q, SE.q, CI.lo, CI.hi
    // ===================================================================
    * --- Table layout parameters ---
    local line_width = min(79, c(linesize))
    local hline_rule `"di as txt \"{hline `line_width'}\""'
    if (c(linesize) < 50) {
        local layout "ultra"
        local lw 4
        local bfmt "%7.3g"
        local bw 7
        local nfmt1 "%8.3g"
        local nfmt2 "%8.3g"
    }
    else if (c(linesize) < 63) {
        local layout "compact"
        local lw 6
        local nfmt1 "%6.3g"
        local nfmt2 "%6.3g"
        local nw1 6
        local nw2 6
        local bfmt "%6.3g"
        local bw 6
        local ci_lo "CI.lo"
        local ci_hi "CI.hi"
        local cb_lo "CB.lo"
        local cb_hi "CB.hi"
    }
    else if (c(linesize) < 79) {
        local layout "narrow"
        local lw 8
        local nfmt1 "%8.4g"
        local nfmt2 "%8.4g"
        local nw1 8
        local nw2 8
        local bfmt "%8.4g"
        local bw 8
        local ci_lo "CI.low"
        local ci_hi "CI.high"
        local cb_lo "CB.low"
        local cb_hi "CB.high"
    }
    else {
        local layout "normal"
        local lw 10
        local bfmt "%9.4g"
        local bw 9
        local ci_lo "CI.low"
        local ci_hi "CI.high"
        local cb_lo "CB.low"
        local cb_hi "CB.high"
        if (`cbands_on') {
            local nfmt1 "%10.4g"
            local nfmt2 "%9.4g"
            local nw1 10
            local nw2 9
        }
        else {
            local nfmt1 "%11.4g"
            local nfmt2 "%11.4g"
            local nw1 11
            local nw2 11
        }
    }
    * --- Derived column headers, widths, and positions (non-ultra) ---
    if ("`layout'" != "ultra") {
        if (`cbands_on') {
            local h1 "Est.q"
            local h2 "SE.q"
            local h3 "`ci_lo'"
            local h4 "`ci_hi'"
            local h5 "`cb_lo'"
            local h6 "`cb_hi'"
            local cw1 `nw1'
            local cw2 `nw2'
            local cw3 `bw'
            local cw4 `bw'
            local cw5 `bw'
            local cw6 `bw'
        }
        else {
            local h1 "Est.p"
            local h2 "SE.p"
            local h3 "Est.q"
            local h4 "SE.q"
            local h5 "`ci_lo'"
            local h6 "`ci_hi'"
            local cw1 `nw1'
            local cw2 `nw1'
            local cw3 `nw1'
            local cw4 `nw1'
            local cw5 `bw'
            local cw6 `bw'
        }
        local c1 = `lw' + 2
        local c2 = `c1' + `cw1' + 1
        local c3 = `c2' + `cw2' + 1
        local c4 = `c3' + `cw3' + 1
        local c5 = `c4' + `cw4' + 1
        local c6 = `c5' + `cw5' + 1
    }

    * --- Summary header ---
    // Determine stdvars display label
    if (`stdvars_active' & "`bwsource'" == "automatic") {
        if ("`stdvars'" != "") local stdvars_label "On (user)"
        else local stdvars_label "On (default)"
    }
    else if (`stdvars_active' & "`bwsource'" == "user") {
        local stdvars_label "On (user)"
    }
    else {
        local stdvars_label "Off"
    }
    di as txt _newline "Location RD estimation"
    di as txt "{hline `line_width'}"
    if ("`layout'" == "ultra" | "`layout'" == "compact") {
        di as txt "  Eval points: " as res %9.0f `neval'
        di as txt "  Observations: " as res %9.0f `N'
        di as txt "  VCE: " as res "`vce'" as txt "  RBC: " as res "`rbc'"
        di as txt "  Side: " as res "`side'"
        di as txt "  Bandwidth: " as res "`bwsource'"
        di as txt "  Kernel: " as res "`kernel'"
        di as txt "  Std. Vars: " as res "`stdvars_label'"
    }
    else {
        di as txt "  Evaluation points: " as res %9.0f `neval' ///
            as txt "    Observations: " as res %9.0f `N'
        di as txt "  VCE: " as res "`vce'" ///
            as txt "    RBC: " as res "`rbc'" ///
            as txt "    Side: " as res "`side'"
        di as txt "  Bandwidth source: " as res "`bwsource'" ///
            as txt "    Kernel: " as res "`kernel'" ///
            as txt "    Std. Vars: " as res "`stdvars_label'"
    }
    di as txt "{hline `line_width'}"

    * --- Column headers ---
    if ("`layout'" == "ultra") {
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
    else {
        di as txt %`lw's "Point" _col(`c1') %`cw1's "`h1'" ///
            _col(`c2') %`cw2's "`h2'" _col(`c3') %`cw3's "`h3'" ///
            _col(`c4') %`cw4's "`h4'" _col(`c5') %`cw5's "`h5'" ///
            _col(`c6') %`cw6's "`h6'"
    }
    di as txt "{hline `line_width'}"

    * --- Data rows ---
    tempname display_results
    matrix `display_results' = e(results)
    forvalues j = 1/`neval' {
        local rname : word `j' of `rownames'
        local suffix : display %02.0f `j'
        local prefix_len = max(0, `lw' - strlen("`suffix'") - 1)
        local dname = cond(strlen("`rname'") > `lw', substr("`rname'", 1, `prefix_len') + "~`suffix'", "`rname'")
        local estp = el(`display_results', `j', 3)
        local sep = el(`display_results', `j', 4)
        local estq = el(`display_results', `j', 5)
        local seq = el(`display_results', `j', 6)
        local cil = el(`display_results', `j', 9)
        local ciu = el(`display_results', `j', 10)
        * Format bounds with infinity handling
        local cil_s : display `bfmt' `cil'
        local ciu_s : display `bfmt' `ciu'
        if (`cil' <= -c(maxdouble) / 2) local cil_s "-inf"
        if (`ciu' >= c(maxdouble) / 2) local ciu_s "inf"
        if (`cbands_on') {
            local cbl = el(`display_results', `j', 11)
            local cbu = el(`display_results', `j', 12)
            local cbl_s : display `bfmt' `cbl'
            local cbu_s : display `bfmt' `cbu'
            if (`cbl' <= -c(maxdouble) / 2) local cbl_s "-inf"
            if (`cbu' >= c(maxdouble) / 2) local cbu_s "inf"
        }
        * Prepare display values (v1-v6 unified for all layouts)
        if (`cbands_on') {
            local v1 : display `nfmt1' `estq'
            local v2 : display `nfmt2' `seq'
            local v3 "`cil_s'"
            local v4 "`ciu_s'"
            local v5 "`cbl_s'"
            local v6 "`cbu_s'"
        }
        else {
            local v1 : display `nfmt1' `estp'
            local v2 : display `nfmt1' `sep'
            local v3 : display `nfmt1' `estq'
            local v4 : display `nfmt1' `seq'
            local v5 "`cil_s'"
            local v6 "`ciu_s'"
        }
        if ("`layout'" == "ultra" & !`cbands_on') local v4 : display %7.3g `seq'
        * Display row
        if ("`layout'" == "ultra") {
            if (`cbands_on') {
                di as res %4s "`dname'" _col(7) %8s "`v1'" _col(17) %8s "`v2'"
                di as res %4s "" _col(7) %7s "`v3'" ///
                    _col(16) %7s "`v4'" _col(25) %7s "`v5'" ///
                    _col(34) %7s "`v6'"
            }
            else {
                di as res %4s "`dname'" _col(7) %8s "`v1'" ///
                    _col(17) %8s "`v2'" _col(27) %8s "`v3'"
                di as res %4s "" _col(7) %7s "`v4'" ///
                    _col(16) %7s "`v5'" _col(25) %7s "`v6'"
            }
        }
        else {
            di as res %`lw's "`dname'" _col(`c1') %`cw1's "`v1'" ///
                _col(`c2') %`cw2's "`v2'" _col(`c3') %`cw3's "`v3'" ///
                _col(`c4') %`cw4's "`v4'" _col(`c5') %`cw5's "`v5'" ///
                _col(`c6') %`cw6's "`v6'"
        }
    }
    di as txt "{hline `line_width'}"

    * --- Footer ---
    if (`cbands_on') {
        if ("`layout'" == "ultra") {
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

program define _rd2d_build_target, rclass
    version 16.0
    syntax, P(integer) [DERiv(numlist min=2 max=2) TANGvec(numlist min=2 max=2)]

    local k = (`p' + 1) * (`p' + 2) / 2
    tempname target
    matrix `target' = J(1, `k', 0)

    if ("`tangvec'" != "") {
        if (`p' < 1) {
            di as err "tangvec() requires p() of at least 1"
            exit 198
        }
        local tv1 : word 1 of `tangvec'
        local tv2 : word 2 of `tangvec'
        if (`tv1' >= . | `tv2' >= . | (`tv1' == 0 & `tv2' == 0)) {
            di as err "tangvec() must contain a nonzero finite direction vector"
            exit 198
        }
        matrix `target'[1, 2] = `tv1'
        matrix `target'[1, 3] = `tv2'
    }
    else {
        local d1 0
        local d2 0
        if ("`deriv'" != "") {
            local d1 : word 1 of `deriv'
            local d2 : word 2 of `deriv'
        }
        local dsum = `d1' + `d2'
        if (`dsum' > `p') {
            di as err "deriv() components must sum to at most p()"
            exit 198
        }
        if (`dsum' == 0) local idx = 1
        else local idx = ((`dsum' + 1) * `dsum' / 2) + `d2' + 1
        local scale = exp(lnfactorial(`d1') + lnfactorial(`d2'))
        matrix `target'[1, `idx'] = `scale'
    }

    return matrix target = `target'
end

mata:
// ---------------------------------------------------------------------------
// _rd2d_loc_bwcheck_bounds: Compute bandwidth bounds from bwcheck
//                           constraint for rd2d estimation
//
// Inputs:
//   x1name     : string scalar  - Stata varname for running var 1
//   x2name     : string scalar  - Stata varname for running var 2
//   tname      : string scalar  - Stata varname for treatment (0/1)
//   tousename  : string scalar  - Stata varname for sample marker
//   b1         : real scalar    - boundary point x1-coordinate
//   b2         : real scalar    - boundary point x2-coordinate
//   ktype      : string scalar  - kernel type (prod/rad)
//   bwcheck    : real scalar    - min observations in bandwidth
//   masspoints : string scalar  - masspoints handling (adjust/off)
//   stdflag    : real scalar    - 1: standardize vars; 0: no
//   outname    : string scalar  - Stata matrix name for output
//
// Output:
//   void - stores (bwmin0, bwmin1, bwmax0, bwmax1) in outname
//
// Dependencies:
//   Called by: rd2d program (Stata caller)
//   Calls:    (Mata built-ins only)
//
// Notes:
//   - bwmin = distance to bwcheck-th nearest obs per side
//   - bwmax = maximum distance per side
//   - When bwcheck<=0, stores all missing values
//   - With masspoints="adjust", uses unique points for distances
// ---------------------------------------------------------------------------
void _rd2d_loc_bwcheck_bounds(
    string scalar x1name,
    string scalar x2name,
    string scalar tname,
    string scalar tousename,
    real scalar b1,
    real scalar b2,
    string scalar ktype,
    real scalar bwcheck,
    string scalar masspoints,
    real scalar stdflag,
    string scalar outname)
{
    real colvector x1, x2, t, idx0, idx1, d0, d1, a0, a1, z0, z1, sort0, sort1
    real matrix X0, X1
    real scalar sd1, sd2, bwmin0, bwmin1, bwmax0, bwmax1

    if (bwcheck <= 0) {
        st_matrix(outname, J(1, 4, .))
        return
    }

    st_view(x1 = ., ., x1name, tousename)
    st_view(x2 = ., ., x2name, tousename)
    st_view(t = ., ., tname, tousename)

    sd1 = 1
    sd2 = 1
    if (stdflag != 0) {
        sd1 = sqrt(variance(x1))
        sd2 = sqrt(variance(x2))
        if (!(sd1 > 0) | !(sd2 > 0)) {
            errprintf("stdvars requires positive sample standard deviations\n")
            _error(198)
        }
        x1 = x1 :/ sd1
        x2 = x2 :/ sd2
        b1 = b1 / sd1
        b2 = b2 / sd2
    }

    idx0 = (t :== 0)
    idx1 = (t :!= 0)
    if (masspoints == "adjust") {
        X0 = uniqrows(select((x1, x2), idx0))
        X1 = uniqrows(select((x1, x2), idx1))
        a0 = X0[, 1] :- b1
        z0 = X0[, 2] :- b2
        a1 = X1[, 1] :- b1
        z1 = X1[, 2] :- b2
    }
    else {
        a0 = select(x1, idx0) :- b1
        z0 = select(x2, idx0) :- b2
        a1 = select(x1, idx1) :- b1
        z1 = select(x2, idx1) :- b2
    }

    if (ktype == "prod") {
        d0 = (abs(a0) :> abs(z0)) :* abs(a0) + (abs(a0) :<= abs(z0)) :* abs(z0)
        d1 = (abs(a1) :> abs(z1)) :* abs(a1) + (abs(a1) :<= abs(z1)) :* abs(z1)
    }
    else {
        d0 = sqrt(a0:^2 + z0:^2)
        d1 = sqrt(a1:^2 + z1:^2)
    }

    sort0 = sort(d0, 1)
    sort1 = sort(d1, 1)
    bwmin0 = sort0[min((bwcheck, rows(sort0)))]
    bwmin1 = sort1[min((bwcheck, rows(sort1)))]
    bwmax0 = sort0[rows(sort0)]
    bwmax1 = sort1[rows(sort1)]

    st_matrix(outname, (bwmin0, bwmin1, bwmax0, bwmax1))
}

// ---------------------------------------------------------------------------
// _rd2d_cb_quantile: Compute sample quantile using linear interp
//
// Inputs:
//   x : real colvector  - numeric vector to compute quantile of
//   q : real scalar     - quantile level in [0, 1]
//
// Output:
//   returns : real scalar  - q-th quantile of x
//
// Dependencies:
//   Called by: _rd2d_cb_critical(), _rd2d_loc_cb_from_corr()
//   Calls:    (Mata built-ins only)
//
// Notes:
//   - Returns missing if x is empty
//   - Boundary: q<=0 returns min, q>=1 returns max
//   - Used to extract critical values from simulated distributions
// ---------------------------------------------------------------------------
real scalar _rd2d_cb_quantile(real colvector x, real scalar q)
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

// ---------------------------------------------------------------------------
// _rd2d_cb_critical: Simulate pointwise critical value for confidence
//                    bands (single target)
//
// Inputs:
//   rep   : real scalar    - number of simulation draws
//   level : real scalar    - confidence level (e.g. 95)
//   side  : string scalar  - "two" (two-sided) or "one" (one-sided)
//
// Output:
//   returns : real scalar  - simulated critical value at level/100
//
// Dependencies:
//   Called by: (utility; not invoked by rd2d program currently)
//   Calls:    _rd2d_cb_quantile()
//
// Notes:
//   - Draws from standard normal, takes abs() for two-sided
//   - Uses _rd2d_cb_quantile to extract the level-percentile
//   - For multi-target bands use _rd2d_loc_cb_from_corr instead
// ---------------------------------------------------------------------------
real scalar _rd2d_cb_critical(real scalar rep, real scalar level, string scalar side)
{
    real colvector sim, tvec

    sim = invnormal(runiform(rep, 1))
    if (side == "two") {
        tvec = abs(sim)
    }
    else {
        tvec = sim
    }

    return(_rd2d_cb_quantile(tvec, level / 100))
}

// ---------------------------------------------------------------------------
// _rd2d_loc_cluster_crossprod: Clustered cross-product of influence
//                              functions for covariance estimation
//
// Inputs:
//   psi : real matrix    - N x k matrix of influence function values
//   C   : real colvector - cluster identifiers (length N)
//
// Output:
//   returns : real matrix  - k x k clustered cross-product matrix
//
// Dependencies:
//   Called by: _rd2d_loc_cov_q()
//   Calls:    (Mata built-ins only)
//
// Notes:
//   - Errors if fewer than 2 clusters
//   - Used for multi-target covariance in confidence band estimation
// ---------------------------------------------------------------------------
real matrix _rd2d_loc_cluster_crossprod(real matrix psi, real colvector C)
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

// ---------------------------------------------------------------------------
// _rd2d_loc_cov_kernel: Evaluate univariate kernel for covariance
//                       estimation (same form as bandwidth kernel)
//
// Inputs:
//   u       : real colvector  - standardized distances (u = x/h)
//   kflag   : real scalar     - kernel index (1=uniform, 2=triangular,
//                                3=epanechnikov, 4=gaussian)
//
// Output:
//   returns : real colvector  - kernel weights K(u), non-negative
//
// Dependencies:
//   Called by: _rd2d_loc_side_influence()
//   Calls:    (Mata built-ins only)
//
// Notes:
//   - Gaussian: exp(-u^2/2)/sqrt(2*pi), infinite support
//   - Separate from rdbw2d kernel for modularity
// ---------------------------------------------------------------------------
real colvector _rd2d_loc_cov_kernel(real colvector u, real scalar kflag)
{
    real colvector a, w

    a = abs(u)
    if (kflag == 1) {
        w = 0.5 :* (a :<= 1 + 1e-12)
    }
    else if (kflag == 2) {
        w = (1 :- a) :* (a :<= 1)
    }
    else if (kflag == 3) {
        w = 0.75 :* (1 :- u:^2) :* (a :<= 1)
    }
    else {
        w = exp(-0.5 :* u:^2) :/ sqrt(2 * pi())
    }

    return(w)
}

// ---------------------------------------------------------------------------
// _rd2d_loc_cov_basis: Construct 2D polynomial basis for covariance
//                      local polynomial fit
//
// Inputs:
//   x1  : real colvector  - first running variable (centered)
//   x2  : real colvector  - second running variable (centered)
//   p   : real scalar     - polynomial order (p >= 0)
//
// Output:
//   returns : real matrix  - N x ((p+1)(p+2)/2) design matrix
//
// Notes:
//   - Same structure as _rdbw2d_loc_basis: intercept + graded
//     lexicographic polynomial terms
//   - Used in _rd2d_loc_side_influence() for influence computation
// ---------------------------------------------------------------------------
real matrix _rd2d_loc_cov_basis(real colvector x1, real colvector x2, real scalar p)
{
    real scalar j, k, count, n, cols
    real matrix X

    n = rows(x1)
    cols = (p + 1) * (p + 2) / 2
    X = J(n, cols, 1)
    count = 2
    for (j = 1; j <= p; j++) {
        for (k = 0; k <= j; k++) {
            X[, count] = (x1:^(j - k)) :* (x2:^k)
            count++
        }
    }

    return(X)
}

// ---------------------------------------------------------------------------
// _rd2d_loc_bread_diag: Diagnose numerical properties of bread matrix
//                       (X'WX) for inversion strategy selection
//
// Inputs:
//   bread : real matrix  - k x k symmetric bread matrix (X'WX)
//
// Output:
//   returns : real rowvector  - (rank, condition, fallback) where:
//     rank      = number of positive eigenvalues (> tol)
//     condition = max/min positive eigenvalue ratio
//     fallback  = 1 if pinv needed, 0 if invsym safe
//
// Notes:
//   - Symmetrizes (bread+bread')/2 before eigendecomposition
//   - tol = 1e-12 * max(|eigenvalues|) for rank determination
//   - fallback=1 when rank<k or condition >= 1e12
//   - Guides _rd2d_loc_bread_inverse() strategy
// ---------------------------------------------------------------------------
real rowvector _rd2d_loc_bread_diag(real matrix bread)
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

// ---------------------------------------------------------------------------
// _rd2d_loc_bread_inverse: Invert bread matrix using diagnostics from
//                          _rd2d_loc_bread_diag
//
// Inputs:
//   bread : real matrix    - k x k symmetric bread matrix
//   diag  : real rowvector - (rank, condition, fallback) from
//                            _rd2d_loc_bread_diag()
//
// Output:
//   returns : real matrix  - k x k inverse of bread
//
// Notes:
//   - If diag[3]=1 (fallback), uses pinv (Moore-Penrose)
//   - Otherwise uses invsym (faster, positive definite assumed)
//   - Symmetrizes before inversion for numerical stability
// ---------------------------------------------------------------------------
real matrix _rd2d_loc_bread_inverse(real matrix bread, real rowvector diag)
{
    real matrix A

    A = (bread + bread') / 2
    if (diag[3] != 0) return(pinv(A))
    return(invsym(A))
}

// ---------------------------------------------------------------------------
// _rd2d_loc_side_influence: Compute per-observation influence function
//                           for one side of the RD estimate
//
// Inputs:
//   yname      : string scalar  - Stata varname for outcome
//   x1name     : string scalar  - Stata varname for running var 1
//   x2name     : string scalar  - Stata varname for running var 2
//   dname      : string scalar  - Stata varname for treatment
//   tousename  : string scalar  - Stata varname for sample marker
//   cname      : string scalar  - cluster var (empty=none)
//   at1        : real scalar    - evaluation point x1-coordinate
//   at2        : real scalar    - evaluation point x2-coordinate
//   h1         : real scalar    - bandwidth for x1 direction
//   h2         : real scalar    - bandwidth for x2 direction
//   p          : real scalar    - polynomial order
//   side       : string scalar  - "control" or "treated"
//   target     : real rowvector - target functional vector
//   kernel     : string scalar  - kernel name
//   ktype      : string scalar  - kernel type (prod/rad)
//   vce        : string scalar  - VCE type
//   stdvars    : real scalar    - 1: standardize; 0: no
//   hascluster : real scalar    - 1 if clustered
//
// Output:
//   returns : real colvector  - N x 1 influence function values
//             (zero for observations outside bandwidth)
//
// Notes:
//   - psi_i = (X'WX)^{-1} * w_i * e_i * target' for obs in window
//   - Applies HC1/HC2/HC3 residual adjustments
//   - Cluster correction scales scores by sqrt(factor)
//   - Full-sample length vector with zeros outside bandwidth
// ---------------------------------------------------------------------------
real colvector _rd2d_loc_side_influence(
    string scalar yname,
    string scalar x1name,
    string scalar x2name,
    string scalar dname,
    string scalar tousename,
    string scalar cname,
    real scalar at1,
    real scalar at2,
    real scalar h1,
    real scalar h2,
    real scalar p,
    string scalar side,
    real rowvector target,
    string scalar kernel,
    string scalar ktype,
    string scalar vce,
    real scalar stdvars,
    real scalar hascluster,
    | real scalar sd1_cached,
      real scalar sd2_cached)
{
    real colvector y, x1, x2, d, C, xc1, xc2, fituse, w, resid, hii, adj, psi, idx
    real colvector clusters, cidx
    real matrix X, bread, ibread, beta, Xw, Xwr, scores
    real rowvector diag
    real scalar nfull, n, k, hr, sd1, sd2, g, j, factor, kflag

    kflag = (kernel == "uniform") * 1 + (kernel == "triangular") * 2 + (kernel == "epanechnikov") * 3 + (kernel == "gaussian") * 4
    if (kflag == 0) kflag = 4

    st_view(y = ., ., yname, tousename)
    st_view(x1 = ., ., x1name, tousename)
    st_view(x2 = ., ., x2name, tousename)
    st_view(d = ., ., dname, tousename)
    if (hascluster) st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)

    nfull = rows(y)

    sd1 = 1
    sd2 = 1
    if (stdvars != 0) {
        if (args() >= 20 & sd1_cached > 0 & sd2_cached > 0) {
            sd1 = sd1_cached
            sd2 = sd2_cached
        }
        else {
            sd1 = sqrt(variance(x1))
            sd2 = sqrt(variance(x2))
        }
        if (!(sd1 > 0) || !(sd2 > 0)) {
            errprintf("stdvars requires positive sample standard deviations\n")
            _error(198)
        }
        x1 = x1 :/ sd1
        x2 = x2 :/ sd2
        at1 = at1 / sd1
        at2 = at2 / sd2
        h1 = h1 / sd1
        h2 = h2 / sd2
    }

    xc1 = x1 :- at1
    xc2 = x2 :- at2

    if (side == "control") fituse = (d :== 0)
    else fituse = (d :!= 0)

    if (ktype == "prod") {
        w = _rd2d_loc_cov_kernel(xc1 :/ h1, kflag) :* _rd2d_loc_cov_kernel(xc2 :/ h2, kflag) :/ (h1 * h2)
    }
    else {
        hr = sqrt(h1^2 + h2^2)
        w = _rd2d_loc_cov_kernel(sqrt(xc1:^2 + xc2:^2) :/ hr, kflag) :/ (hr^2)
    }

    fituse = fituse :& (w :> 0)
    idx = selectindex(fituse)
    y = select(y, fituse)
    xc1 = select(xc1, fituse)
    xc2 = select(xc2, fituse)
    w = select(w, fituse)
    C = select(C, fituse)

    n = rows(y)
    k = (p + 1) * (p + 2) / 2
    if (n < k) {
        errprintf("not enough observations for p(%g) on %s side\n", p, side)
        _error(2001)
    }
    if (vce == "hc1" & n <= k) {
        errprintf("not enough effective observations for hc1 covariance correction\n")
        _error(2001)
    }

    X = _rd2d_loc_cov_basis(xc1, xc2, p)
    bread = quadcross(X, w, X)
    diag = _rd2d_loc_bread_diag(bread)
    ibread = _rd2d_loc_bread_inverse(bread, diag)
    beta = ibread * quadcross(X, w, y)
    resid = y - X * beta

    if (vce == "hc1" & !hascluster) {
        resid = resid :* sqrt(n / (n - k))
    }
    else if (vce == "hc2" | vce == "hc3") {
        Xw = sqrt(w) :* X
        hii = diagonal(Xw * ibread * Xw')
        if (any(hii :>= 1)) {
            errprintf("leverage adjustment undefined because hat values reach 1\n")
            _error(498)
        }
        if (vce == "hc2") resid = resid :* sqrt(1 :/ (1 :- hii))
        else resid = resid :* (1 :/ (1 :- hii))
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

    psi = J(nfull, 1, 0)
    psi[idx] = scores * target'
    return(psi)
}

// ---------------------------------------------------------------------------
// _rd2d_loc_cov_q: Compute multi-target covariance and correlation
//                  matrices for confidence band construction
//
// Inputs:
//   yname      : string scalar  - Stata varname for outcome
//   x1name     : string scalar  - Stata varname for running var 1
//   x2name     : string scalar  - Stata varname for running var 2
//   dname      : string scalar  - Stata varname for treatment
//   tousename  : string scalar  - Stata varname for sample marker
//   bwsname    : string scalar  - Stata matrix of boundary points
//   hqname     : string scalar  - Stata matrix of bandwidths per point
//   targetname : string scalar  - Stata matrix for target vector
//   p          : real scalar    - polynomial order
//   kernel     : string scalar  - kernel name
//   ktype      : string scalar  - kernel type (prod/rad)
//   vce        : string scalar  - VCE type
//   stdflag    : string scalar  - "1"/"0" for stdvars
//   cname      : string scalar  - cluster var (empty=none)
//   covname    : string scalar  - output matrix for covariance
//   corrname   : string scalar  - output matrix for correlation
//
// Output:
//   void - stores k x k covariance and correlation in named matrices
//
// Notes:
//   - Computes influence functions per target point on both sides
//   - psi = psi_treated - psi_control
//   - Covariance: cross-product (or clustered) of psi columns
//   - Correlation: normalized by diagonal standard errors
// ---------------------------------------------------------------------------
void _rd2d_loc_cov_q(
    string scalar yname,
    string scalar x1name,
    string scalar x2name,
    string scalar dname,
    string scalar tousename,
    string scalar bwsname,
    string scalar hqname,
    string scalar targetname,
    real scalar p,
    string scalar kernel,
    string scalar ktype,
    string scalar vce,
    string scalar stdflag,
    string scalar cname,
    string scalar covname,
    string scalar corrname)
{
    real colvector y, d, C, se
    real matrix bws, hq, target, psit, psic, psi, cov, corr
    real scalar j, l, k, hascluster

    st_view(y = ., ., yname, tousename)
    st_view(d = ., ., dname, tousename)
    hascluster = (cname != "")
    if (hascluster) st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)
    bws = st_matrix(bwsname)
    hq = st_matrix(hqname)
    target = st_matrix(targetname)
    k = rows(bws)

    // Memory safety check: psit + psic = 2 * N * k * 8 bytes
    if (2 * rows(y) * k * 8 / (1024 * 1024) > 2048) {
        printf("{txt}note: covariance computation requires %g MB for %g evaluation points\n",
               2 * rows(y) * k * 8 / (1024 * 1024), k)
        printf("{txt}      consider reducing the number of evaluation points if memory is limited\n")
    }

    psit = J(rows(y), k, 0)
    psic = J(rows(y), k, 0)

    // Pre-compute stdvars constants to avoid redundant variance() calls
    real scalar sd1_pre, sd2_pre
    real colvector x1_tmp, x2_tmp
    sd1_pre = 0
    sd2_pre = 0
    if (strtoreal(stdflag) != 0) {
        st_view(x1_tmp = ., ., x1name, tousename)
        st_view(x2_tmp = ., ., x2name, tousename)
        sd1_pre = sqrt(variance(x1_tmp))
        sd2_pre = sqrt(variance(x2_tmp))
    }

    for (j = 1; j <= k; j++) {
        psit[, j] = _rd2d_loc_side_influence(yname, x1name, x2name, dname, tousename, cname, ///
            bws[j, 1], bws[j, 2], hq[j, 3], hq[j, 4], p, "treated", target, ///
            kernel, ktype, vce, strtoreal(stdflag), hascluster, sd1_pre, sd2_pre)
        psic[, j] = _rd2d_loc_side_influence(yname, x1name, x2name, dname, tousename, cname, ///
            bws[j, 1], bws[j, 2], hq[j, 1], hq[j, 2], p, "control", target, ///
            kernel, ktype, vce, strtoreal(stdflag), hascluster, sd1_pre, sd2_pre)
    }

    psi = psit - psic
    if (!hascluster) cov = quadcross(psi, psi)
    else cov = _rd2d_loc_cluster_crossprod(psi, C)
    se = sqrt(diagonal(cov))
    corr = J(k, k, .)
    for (j = 1; j <= k; j++) {
        for (l = 1; l <= k; l++) {
            if (se[j] * se[l] > 0 & se[j] * se[l] < .) corr[j, l] = cov[j, l] / (se[j] * se[l])
        }
    }
    corr = (corr + corr') / 2
    for (j = 1; j <= k; j++) {
        if (se[j] > 0 & se[j] < .) corr[j, j] = 1
    }

    st_matrix(covname, cov)
    st_matrix(corrname, corr)
}

// ---------------------------------------------------------------------------
// _rd2d_loc_cb_from_corr: Simulate critical value for simultaneous
//                         confidence bands from correlation matrix
//
// Inputs:
//   corrname     : string scalar  - Stata matrix of correlation
//   rep          : real scalar    - number of simulation draws
//   level        : real scalar    - confidence level (e.g. 95)
//   side         : string scalar  - "two" or "one"
//   critname     : string scalar  - output scalar for critical value
//   adjustedname : string scalar  - output scalar for PSD flag
//   mineigname   : string scalar  - output scalar for min eigenvalue
//
// Output:
//   void - stores critical value, PSD-adjustment flag, and min
//          eigenvalue in named Stata scalars
//
// Notes:
//   - Eigendecomposes correlation; clamps negatives to 0
//   - Draws MVN: sim = Z * L' where L uses sqrt(eigval)
//   - Two-sided: quantile of max(|sim|) across columns
//   - One-sided: quantile of max(sim) across columns
//   - Reports whether PSD adjustment was applied
// ---------------------------------------------------------------------------
void _rd2d_loc_cb_from_corr(
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

    st_numscalar(critname, _rd2d_cb_quantile(tvec, level / 100))
    st_numscalar(adjustedname, adjusted)
    st_numscalar(mineigname, mineig)
}
end
