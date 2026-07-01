*! version 1.1.0 24jun2026
program define rdbw2d, rclass
    version 16.0
    if `"`1'"' == "version" | `"`1'"' == ",version" {
        di as txt "rdbw2d version 1.1.0 (24 June 2026)"
        exit
    }
    syntax varlist(min=4 max=4 numeric) [if] [in], AT(string asis) ///
        [P(integer 1) DERiv(string asis) TANGvec(string asis) KERnel(string) KTYPE(string) BWSELect(string) METHOD(string) ///
         VCE(string) BWCHeck(string) MASSPoints(string) C(varname) ///
         CLuster(varname) ///
         SCALEregul(real 3) SCALEbiascrct(real 1) STDVars RAWVARS ///
         MP(string) SR(real -1) BWS(string) KT(string)]

    gettoken yvar rest : varlist
    gettoken x1var rest : rest
    gettoken x2var tvar : rest
    local clustername "`c'"
    if ("`cluster'" != "") {
        if ("`c'" != "" & "`c'" != "`cluster'") {
            di as err "cluster() and c() must refer to the same variable when both are specified"
            exit 198
        }
        local clustername "`cluster'"
    }
    local clusterwork "`clustername'"
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
    // --- stdvars / rawvars conflict check ---
    if ("`stdvars'" != "" & "`rawvars'" != "") {
        di as err "cannot specify both stdvars and rawvars"
        exit 198
    }

    if (`p' < 0) {
        di as err "p() must be a nonnegative integer (0, 1, 2, ...); got `p'"
        exit 198
    }
    if (`scaleregul' >= . | `scaleregul' < 0) {
        di as err "scaleregul() must be finite and nonnegative"
        exit 198
    }
    if (`scalebiascrct' >= . | `scalebiascrct' < 0) {
        di as err "scalebiascrct() must be finite and nonnegative"
        exit 198
    }

    local d1 0
    local d2 0
    local derivsum 0
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
            di as err "deriv() components must be nonnegative and sum to at most p(); got d1=`d1' d2=`d2' with p=`p'"
            exit 198
        }
        local derivsum = `d1' + `d2'
    }

    local tv1 .
    local tv2 .
    if ("`tangvec'" != "") {
        local tvlist "`tangvec'"
        local tvcount : word count `tvlist'
        if (`tvcount' != 2) {
            di as err "tangvec() must contain exactly two numbers"
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
        local d1 1
        local d2 0
        local derivsum 1
        if ("`deriv'" != "") di as txt "warning: tangvec() provided; deriv() is ignored."
    }

    local kernel = lower("`kernel'")
    if ("`kernel'" == "") local kernel "triangular"
    if inlist("`kernel'", "uni", "unif") local kernel "uniform"
    if inlist("`kernel'", "tri", "triag") local kernel "triangular"
    if inlist("`kernel'", "epa", "epan") local kernel "epanechnikov"
    if inlist("`kernel'", "gau") local kernel "gaussian"
    if !inlist("`kernel'", "uniform", "triangular", "epanechnikov", "gaussian") {
        di as err "kernel() must be uniform, triangular, epanechnikov, or gaussian (aliases: uni, tri, epa, gau); got '`kernel''"
        exit 198
    }

    local ktype = lower("`ktype'")
    if ("`ktype'" == "") local ktype "prod"
    if !inlist("`ktype'", "prod", "rad") {
        di as err "ktype() must be prod or rad"
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
        di as err "bwselect() must be mserd, imserd, msetwo, or imsetwo; got '`bwselect''"
        exit 198
    }

    // method options:
    //   dpi = direct plug-in (uses pilot bandwidth for bias estimation)
    //   rot = rule-of-thumb (uses normal-reference formula without pilot fit)
    local method = lower("`method'")
    if ("`method'" == "") local method "dpi"
    if !inlist("`method'", "dpi", "rot") {
        di as err "method() must be 'dpi' (direct plug-in) or 'rot' (rule-of-thumb); got '`method''"
        exit 198
    }

    local vce = lower("`vce'")
    if ("`vce'" == "") local vce "hc1"
    if !inlist("`vce'", "hc0", "hc1", "hc2", "hc3") {
        di as err "vce() must be hc0, hc1, hc2, or hc3"
        exit 198
    }
    if ("`clustername'" != "" & !inlist("`vce'", "hc0", "hc1")) {
        di as txt "note: vce(`vce') not available with cluster(); using vce(hc1) cluster-robust"
        local vce "hc1"
    }

    local masspoints = lower("`masspoints'")
    if ("`masspoints'" == "") local masspoints "check"
    if !inlist("`masspoints'", "check", "adjust", "off") {
        di as err "masspoints() must be check, adjust, or off"
        exit 198
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

    // at() specifies the boundary point coordinates (b1, b2)
    local atlist "`at'"
    local atcount : word count `atlist'
    if (`atcount' != 2) {
        di as err "at() requires a pair of coordinates (b1 b2); got `atcount' values"
        exit 198
    }
    local b1 : word 1 of `atlist'
    local b2 : word 2 of `atlist'
    capture confirm number `b1'
    if (_rc) {
        di as err "at() requires a pair of numeric coordinates (b1 b2); first value is not a number"
        exit 198
    }
    capture confirm number `b2'
    if (_rc) {
        di as err "at() requires a pair of numeric coordinates (b1 b2); second value is not a number"
        exit 198
    }
    if (`b1' >= . | `b2' >= .) {
        di as err "at() requires a pair of finite coordinates (b1 b2); got missing/infinite value"
        exit 198
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
    quietly count if `touse' & `tvar' == 0
    local N0 = r(N)
    quietly count if `touse' & `tvar' == 1
    local N1 = r(N)

    // stdvars: when enabled, standardize running variables by their sample SD
    // before bandwidth selection; returned bandwidths are back-transformed
    // Default: ON (stdvars is default for automatic bandwidth) unless rawvars
    local stdflag = 1
    if ("`rawvars'" != "") local stdflag = 0
    if ("`stdvars'" != "") local stdflag = 1

    local k = (`p' + 1) * (`p' + 2) / 2
    tempname target bws mseconsts massinfo
    matrix `target' = J(1, `k', 0)
    if ("`tangvec'" != "") {
        matrix `target'[1, 2] = `tv1'
        matrix `target'[1, 3] = `tv2'
    }
    else {
        if (`derivsum' == 0) local targetidx = 1
        else local targetidx = ((`derivsum' + 1) * `derivsum' / 2) + `d2' + 1
        local targetscale = exp(lnfactorial(`d1') + lnfactorial(`d2'))
        matrix `target'[1, `targetidx'] = `targetscale'
    }

    mata: _rdbw2d_loc_mata("`yvar'", "`x1var'", "`x2var'", "`tvar'", "`touse'", ///
        `b1', `b2', `p', "`kernel'", "`ktype'", "`bwselect'", "`method'", "`vce'", ///
        `bwcheck', "`masspoints'", `scaleregul', `scalebiascrct', `stdflag', "`clusterwork'", ///
        "`target'", `derivsum', "`bws'", "`mseconsts'", "`massinfo'")

    matrix colnames `bws' = b1 b2 h01 h02 h11 h12
    matrix colnames `mseconsts' = Nh.0 Nh.1 bias.0 bias.1 var.0 var.1 reg.bias.0 reg.bias.1 reg.var.0 reg.var.1
    // massinfo matrix columns: [M, M0, M1, mass_ratio]
    // M = total unique (x1,x2) pairs, M0/M1 = unique per side
    // mass_ratio = 1 - M/N; values >= 0.2 trigger warnings
    matrix colnames `massinfo' = M M0 M1 mass
    matrix rownames `bws' = 1
    matrix rownames `mseconsts' = 1
    matrix rownames `massinfo' = 1

    _rd2d_masspoints_warn, matrix(`massinfo') neval(1) masspoints(`masspoints')

    local db1 = el(`bws', 1, 1)
    local db2 = el(`bws', 1, 2)
    local dh01 = el(`bws', 1, 3)
    local dh02 = el(`bws', 1, 4)
    local dh11 = el(`bws', 1, 5)
    local dh12 = el(`bws', 1, 6)
    local stdstatus = cond(`stdflag', "on", "off")
    // Determine stdvars display label
    if (`stdflag') {
        if ("`stdvars'" != "") local stdlabel "On (user)"
        else local stdlabel "On (default)"
    }
    else {
        local stdlabel "Off"
    }
    local line_width = min(79, c(linesize))
    local ultra_table = c(linesize) < 50
    local compact_table = c(linesize) < 63
    local narrow_table = c(linesize) < 79
    di as txt _newline "Location bandwidth selection"
    di as txt "{hline `line_width'}"
    if (`ultra_table' | `compact_table') {
        di as txt "  Observations: " as res %9.0f `N'
        di as txt "  Selector: " as res "`bwselect'" as txt "  Method: " as res "`method'"
        di as txt "  Kernel: " as res "`kernel'"
        di as txt "  Ktype: " as res "`ktype'" as txt "  VCE: " as res "`vce'"
        di as txt "  Std. Vars: " as res "`stdlabel'"
    }
    else {
        di as txt "  Observations: " as res %9.0f `N' ///
            as txt "    Selector: " as res "`bwselect'" ///
            as txt "    Method: " as res "`method'"
        di as txt "  Kernel: " as res "`kernel'" ///
            as txt "    Ktype: " as res "`ktype'" ///
            as txt "    VCE: " as res "`vce'" ///
            as txt "    Std. Vars: " as res "`stdlabel'"
    }
    di as txt "{hline `line_width'}"
    if (`ultra_table') {
        di as txt %8s "b1" _col(10) %8s "b2" ///
            _col(20) %8s "h01" _col(30) %8s "h02"
        di as res %8.2g `db1' _col(10) %8.2g `db2' ///
            _col(20) %8.3g `dh01' _col(30) %8.3g `dh02'
        di as txt %8s "" _col(20) %8s "h11" _col(30) %8s "h12"
        di as res %8s "" _col(20) %8.3g `dh11' _col(30) %8.3g `dh12'
    }
    else if (`compact_table') {
        di as txt %6s "b1" _col(10) %6s "b2" _col(18) %6s "h01" ///
            _col(25) %6s "h02" _col(32) %6s "h11" _col(39) %6s "h12"
        di as txt "{hline `line_width'}"
        di as res %6.0g `db1' _col(10) %6.0g `db2' ///
            _col(18) %6.3g `dh01' _col(25) %6.3g `dh02' ///
            _col(32) %6.3g `dh11' _col(39) %6.3g `dh12'
    }
    else if (`narrow_table') {
        di as txt %8s "b1" _col(10) %8s "b2" _col(19) %8s "h01" ///
            _col(28) %8s "h02" _col(37) %8s "h11" _col(46) %8s "h12"
        di as txt "{hline `line_width'}"
        di as res %8.4g `db1' _col(10) %8.4g `db2' ///
            _col(19) %8.4g `dh01' _col(28) %8.4g `dh02' ///
            _col(37) %8.4g `dh11' _col(46) %8.4g `dh12'
    }
    else {
        di as txt %8s "b1" _col(11) %8s "b2" _col(22) %10s "h01" ///
            _col(34) %10s "h02" _col(46) %10s "h11" _col(58) %10s "h12"
        di as txt "{hline `line_width'}"
        di as res %8.4g `db1' _col(11) %8.4g `db2' ///
            _col(22) %10.4g `dh01' _col(34) %10.4g `dh02' ///
            _col(46) %10.4g `dh11' _col(58) %10.4g `dh12'
    }
    di as txt "{hline `line_width'}"

    return matrix bws = `bws'
    return matrix mseconsts = `mseconsts'
    return matrix masspoints = `massinfo'
    return scalar N = `N'
    return scalar N0 = `N0'
    return scalar N1 = `N1'
    return scalar p = `p'
    return scalar derivsum = `derivsum'
    return scalar bwcheck = `bwcheck'
    return scalar scaleregul = `scaleregul'
    return scalar scalebiascrct = `scalebiascrct'
    return local kernel "`kernel'"
    return local ktype "`ktype'"
    return local bwselect "`bwselect'"
    return local method "`method'"
    return local vce "`vce'"
    return local masspoints_opt "`masspoints'"
    return local stdvars = cond(`stdflag', "on", "off")
    return local version "1.1.0"
    if ("`tangvec'" != "") return local deriv "0 0"
    else return local deriv "`d1' `d2'"
    return local tangvec "`tangvec'"
end

mata:
// ---------------------------------------------------------------------------
// _rdbw2d_loc_kernel: Evaluate univariate kernel function
//
// Inputs:
//   u       : real colvector  - standardized distances (u = x/h)
//   kernel  : string scalar   - kernel name: uniform, triangular,
//                                epanechnikov, or gaussian
//
// Output:
//   returns : real colvector  - kernel weights K(u), non-negative
//
// Notes:
//   - Compact support: uniform/triangular/epanechnikov return 0
//     for |u|>1
//   - Gaussian has infinite support (no truncation applied)
//   - Used by _rdbw2d_loc_weights() for product/radial kernel
// ---------------------------------------------------------------------------
real colvector _rdbw2d_loc_kernel(real colvector u, string scalar kernel)
{
    real colvector a, w
    a = abs(u)
    if (kernel == "uniform") {
        w = 0.5 :* (a :<= 1 + 1e-12)
    }
    else if (kernel == "triangular") {
        w = (1 :- a) :* (a :<= 1)
    }
    else if (kernel == "epanechnikov") {
        w = 0.75 :* (1 :- u:^2) :* (a :<= 1)
    }
    else {
        w = exp(-0.5 :* u:^2) :/ sqrt(2 * pi())
    }
    return(w)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_basis: Construct 2D polynomial basis matrix
//
// Inputs:
//   x1  : real colvector  - first running variable (centered, scaled)
//   x2  : real colvector  - second running variable (centered, scaled)
//   p   : real scalar     - polynomial order (p >= 0)
//
// Output:
//   returns : real matrix  - N x ((p+1)(p+2)/2) design matrix;
//             columns: 1, x1, x2, x1^2, x1*x2, x2^2, ...
//
// Notes:
//   - Column 1 is always the intercept (all ones)
//   - Graded lexicographic ordering: degree j terms are
//     x1^(j-k) * x2^k for k=0..j
//   - Used by all local polynomial fitting routines
// ---------------------------------------------------------------------------
real matrix _rdbw2d_loc_basis(real colvector x1, real colvector x2, real scalar p)
{
    real scalar j, k, count, cols
    real matrix X
    cols = (p + 1) * (p + 2) / 2
    X = J(rows(x1), cols, 1)
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
// _rdbw2d_loc_rot: Compute rule-of-thumb pilot bandwidth
//
// Inputs:
//   x1      : real colvector  - first running variable (full sample)
//   x2      : real colvector  - second running variable (full sample)
//   kernel  : string scalar   - kernel name for constant lookup
//
// Output:
//   returns : real scalar     - scalar pilot bandwidth h_rot > 0,
//                               or missing (.) if degenerate
//
// Notes:
//   - Based on 2D normal-reference rule: uses sample covariance
//     matrix S to approximate integrated squared second derivative
//   - Formula: h_rot = ((2*R(K)) / (N * mu2(K)^2 * C_f))^(1/6)
//   - Returns missing if S is singular or N=0
//   - Provides initial bandwidth for DPI iteration
// ---------------------------------------------------------------------------
real scalar _rdbw2d_loc_rot(real colvector x1, real colvector x2, string scalar kernel)
{
    real matrix S
    real scalar mu2K, l2K, N, traceconst, tr1

    S = variance((x1, x2))
    if (rows(S) != 2 | cols(S) != 2) return(.)
    if (kernel == "epanechnikov") {
        mu2K = 1 / 6
        l2K = 4 / (3 * pi())
    }
    else if (kernel == "triangular") {
        mu2K = 3 / 20
        l2K = 3 / (2 * pi())
    }
    else if (kernel == "uniform") {
        mu2K = 1 / 4
        l2K = 1 / pi()
    }
    else {
        mu2K = 1
        l2K = 1 / (4 * pi())
    }

    N = rows(x1)
    tr1 = trace(invsym(S))
    traceconst = 1 / (2^(2 + 2) * pi()^(1) * det(cholesky(S))) * ///
        (2 * trace(invsym(S) * invsym(S)) + tr1^2)
    if (!(traceconst > 0)) return(.)
    return(((2 * l2K) / (N * mu2K * traceconst))^(1 / 6))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_H: Build diagonal bandwidth scaling matrix
//
// Inputs:
//   h       : real scalar  - bandwidth value (h > 0)
//   p       : real scalar  - polynomial order
//   inverse : real scalar  - 0: H = diag(1, h, h, h^2, ...)
//                            nonzero: H^{-1} = diag(1, 1/h, ...)
//
// Output:
//   returns : real matrix  - ((p+1)(p+2)/2) x ((p+1)(p+2)/2) diagonal
//
// Notes:
//   - Entry (j,j) corresponds to total degree of basis column j
//   - H converts coefficients between scaled (u=x/h) and original
//     coordinates: beta_orig = H^{-1} * beta_scaled
//   - Used in _rdbw2d_loc_lm_beta_cov() and _rdbw2d_loc_bwconst()
// ---------------------------------------------------------------------------
real matrix _rdbw2d_loc_H(real scalar h, real scalar p, real scalar inverse)
{
    real scalar j, k, count, cols, val
    real matrix H

    cols = (p + 1) * (p + 2) / 2
    H = I(cols)
    count = 2
    for (j = 1; j <= p; j++) {
        for (k = 0; k <= j; k++) {
            val = h^j
            if (inverse != 0) val = 1 / val
            H[count, count] = val
            count++
        }
    }
    return(H)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_weights: Compute kernel weights for 2D local polynomial
//
// Inputs:
//   x1      : real colvector  - first running variable (centered)
//   x2      : real colvector  - second running variable (centered)
//   dist    : real colvector  - precomputed distances (radial mode)
//   h       : real scalar     - bandwidth (h > 0)
//   kernel  : string scalar   - kernel name
//   ktype   : string scalar   - "prod" (product) or "rad" (radial)
//
// Output:
//   returns : real colvector  - non-negative weights, normalized by h^2
//
// Notes:
//   - Product kernel: K(x1/h)*K(x2/h)/h^2
//   - Radial kernel:  K(dist/h)/h^2
//   - Zero weights mark observations outside the bandwidth
// ---------------------------------------------------------------------------
real colvector _rdbw2d_loc_weights(real colvector x1, real colvector x2,
    real colvector dist, real scalar h, string scalar kernel, string scalar ktype)
{
    if (ktype == "prod") {
        return(_rdbw2d_loc_kernel(x1 :/ h, kernel) :* _rdbw2d_loc_kernel(x2 :/ h, kernel) :/ (h^2))
    }
    return(_rdbw2d_loc_kernel(vec(dist) :/ h, kernel) :/ (h^2))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_distance: Compute distance metric for kernel evaluation
//
// Inputs:
//   x1    : real colvector  - first coordinate (centered at boundary)
//   x2    : real colvector  - second coordinate (centered at boundary)
//   ktype : string scalar   - "rad": Euclidean; "prod": Chebyshev
//
// Output:
//   returns : real colvector  - non-negative scalar distances
//
// Notes:
//   - Radial: sqrt(x1^2 + x2^2)
//   - Product/Chebyshev: max(|x1|, |x2|)
//   - Inputs are vectorized (vec()) before computation
//   - Determines which observations fall within bandwidth h
// ---------------------------------------------------------------------------
real colvector _rdbw2d_loc_distance(real colvector x1, real colvector x2, string scalar ktype)
{
    real colvector ax1, ax2

    x1 = vec(x1)
    x2 = vec(x2)
    if (ktype == "rad") return(sqrt(x1:^2 + x2:^2))

    ax1 = abs(x1)
    ax2 = abs(x2)
    return((ax1 :> ax2) :* ax1 :+ (ax1 :<= ax2) :* ax2)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_median: Compute sample median
//
// Inputs:
//   x : real colvector  - numeric vector (may contain any reals)
//
// Output:
//   returns : real scalar  - median value; missing (.) if x is empty
//
// Notes:
//   - For odd n: middle element of sorted x
//   - For even n: average of two middle elements
//   - Used to set pilot bandwidths for bias estimation
// ---------------------------------------------------------------------------
real scalar _rdbw2d_loc_median(real colvector x)
{
    real scalar n
    real colvector sx

    sx = sort(x, 1)
    n = rows(sx)
    if (n == 0) return(.)
    if (mod(n, 2) == 1) return(sx[(n + 1) / 2])
    return((sx[n / 2] + sx[n / 2 + 1]) / 2)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_clamp: Clamp bandwidth to [lo, hi] interval
//
// Inputs:
//   h   : real scalar  - candidate bandwidth
//   lo  : real scalar  - lower bound (missing = no lower bound)
//   hi  : real scalar  - upper bound (missing = no upper bound)
//
// Output:
//   returns : real scalar  - h clamped to [lo, hi]
//
// Notes:
//   - Bounds are skipped when set to missing (.)
//   - Enforces bwcheck constraints: ensures minimum observations
//     within the bandwidth window
// ---------------------------------------------------------------------------
real scalar _rdbw2d_loc_clamp(real scalar h, real scalar lo, real scalar hi)
{
    if (lo < . & h < lo) h = lo
    if (hi < . & h > hi) h = hi
    return(h)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_inv: Robust matrix inverse with fallback chain
//
// Inputs:
//   A : real matrix  - square matrix to invert
//
// Output:
//   returns : real matrix  - inverse (or pseudo-inverse) of A
//
// Notes:
//   - Attempts QR inverse first (qrinv)
//   - Falls back to invsym if QR produces missing values
//   - Falls back to Moore-Penrose pseudo-inverse as last resort
//   - Handles near-singular design matrices gracefully
// ---------------------------------------------------------------------------
real matrix _rdbw2d_loc_inv(real matrix A)
{
    real matrix B

    B = qrinv(A)
    if (hasmissing(B)) B = invsym(A)
    if (hasmissing(B)) B = pinv(A)
    return(B)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_fullrank: Test whether symmetric matrix has full rank
//
// Inputs:
//   A : real matrix  - square symmetric matrix (symmetrized internally)
//
// Output:
//   returns : real scalar  - 1 if full rank, 0 otherwise
//
// Notes:
//   - Symmetrizes A = (A+A')/2 before eigendecomposition
//   - Relative tolerance: tol = 1e-12 * max(|eigenvalues|)
//   - Full rank iff all eigenvalues exceed tol
//   - Used to verify design matrix rank before bandwidth estimation
// ---------------------------------------------------------------------------
real scalar _rdbw2d_loc_fullrank(real matrix A)
{
    real matrix eigvec
    real rowvector eigval, pos
    real scalar maxeig, tol

    A = (A + A') / 2
    symeigensystem(A, eigvec, eigval)
    maxeig = max(abs(eigval))
    if (!(maxeig > 0)) return(0)
    tol = 1e-12 * maxeig
    pos = select(eigval, eigval :> tol)
    return(length(pos) == cols(A))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_meat: Compute meat of the sandwich variance estimator
//
// Inputs:
//   X          : real matrix    - N x k design matrix (scaled basis)
//   w          : real colvector - kernel weights
//   resid      : real colvector - residuals (HC-adjusted)
//   C          : real colvector - cluster identifiers
//   hascluster : real scalar    - 1 if cluster VCE, 0 otherwise
//   h          : real scalar    - bandwidth (for normalization)
//   vce        : string scalar  - VCE type (hc0/hc1/hc2/hc3)
//
// Output:
//   returns : real matrix  - k x k meat matrix, scaled by h^2
//
// Notes:
//   - Non-clustered: sum of outer products of weighted scores
//   - Clustered: aggregates scores within clusters then forms
//     outer products; applies finite-sample correction for hc1
//   - Errors if fewer than 2 clusters found
// ---------------------------------------------------------------------------
real matrix _rdbw2d_loc_meat(real matrix X, real colvector w, real colvector resid,
    real colvector C, real scalar hascluster, real scalar h, string scalar vce)
{
    real scalar n, k, g, factor, j
    real matrix meat, score
    real colvector clusters, cidx

    n = rows(X)
    k = cols(X)
    if (!hascluster) {
        return(quadcross(X :* (w :* resid), X :* (w :* resid)) * h^2)
    }

    clusters = uniqrows(sort(C, 1))
    g = rows(clusters)
    if (g <= 1) {
        errprintf("cluster() must contain at least two clusters\n")
        _error(2001)
    }
    factor = 1
    if (vce == "hc1") factor = ((n - 1) / (n - k)) * (g / (g - 1))
    meat = J(k, k, 0)
    for (j = 1; j <= g; j++) {
        cidx = (C :== clusters[j])
        score = quadcross(select(X, cidx), select(w, cidx) :* select(resid, cidx))
        meat = meat + score * score'
    }
    return(meat * h^2 * factor)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_coeff: Estimate pilot bias coefficients for bandwidth
//                    constant computation
//
// Inputs:
//   x1     : real colvector  - first running var (centered)
//   x2     : real colvector  - second running var (centered)
//   dist   : real colvector  - distances (recomputed internally)
//   p      : real scalar     - polynomial order of main fit
//   target : real rowvector  - target functional vector e_nu
//   h      : real scalar     - pilot bandwidth
//   kernel : string scalar   - kernel name
//   ktype  : string scalar   - kernel type (prod/rad)
//
// Output:
//   returns : real rowvector  - 1 x kp1 bias direction vector;
//             all missing if insufficient observations
//
// Notes:
//   - Fits order-(p+1) polynomial to extract bias direction
//   - Returns vecq = (0, target * G^{-1} * X'WS) where S is the
//     higher-order basis block
// ---------------------------------------------------------------------------
real rowvector _rdbw2d_loc_coeff(real colvector x1, real colvector x2,
    real colvector dist, real scalar p, real rowvector target, real scalar h,
    string scalar kernel, string scalar ktype)
{
    real colvector w
    real matrix Raug, R, S, invG
    real colvector idx
    real scalar kp, kp1

    x1 = vec(x1)
    x2 = vec(x2)
    dist = _rdbw2d_loc_distance(x1, x2, ktype)
    w = vec(_rdbw2d_loc_weights(x1, x2, dist, h, kernel, ktype))
    idx = selectindex(w :> 0)
    if (length(idx) < (p + 1) * (p + 2) / 2) return(J(1, (p + 2) * (p + 3) / 2, .))

    x1 = x1[idx] :/ h
    x2 = x2[idx] :/ h
    w = w[idx]

    kp = (p + 1) * (p + 2) / 2
    kp1 = (p + 2) * (p + 3) / 2
    Raug = _rdbw2d_loc_basis(x1, x2, p + 1)
    R = Raug[, 1..kp]
    S = Raug[, (kp + 1)..kp1]
    invG = _rdbw2d_loc_inv(quadcross(R, w, R))

    return((J(1, kp, 0), target * invG * quadcross(R, w, S)))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_lm_beta_cov: Weighted local polynomial fit with optional
//                           sandwich covariance
//
// Inputs:
//   y          : real colvector  - outcome variable
//   x1         : real colvector  - first running variable (centered)
//   x2         : real colvector  - second running variable (centered)
//   dist       : real colvector  - distances (recomputed internally)
//   C          : real colvector  - cluster variable
//   hascluster : real scalar     - 1 if clustered, 0 otherwise
//   p          : real scalar     - polynomial order
//   h          : real scalar     - bandwidth
//   vce        : string scalar   - VCE type (hc0/hc1/hc2/hc3)
//   kernel     : string scalar   - kernel name
//   ktype      : string scalar   - kernel type (prod/rad)
//   wantcov    : real scalar     - 1: return beta+cov; 0: beta only
//
// Output:
//   returns : real matrix  - wantcov=0: k x 1 coefficient vector;
//             wantcov=1: (2k) x k matrix with beta and covariance
//
// Notes:
//   - Rescales coordinates by h before fitting (u = x/h)
//   - Applies HC0/HC1/HC2/HC3 residual adjustment
//   - Coefficients back-transformed to original scale via H^{-1}
// ---------------------------------------------------------------------------
real matrix _rdbw2d_loc_lm_beta_cov(real colvector y, real colvector x1,
    real colvector x2, real colvector dist, real colvector C, real scalar hascluster,
    real scalar p, real scalar h, string scalar vce, string scalar kernel,
    string scalar ktype, real scalar wantcov)
{
    real colvector yv, x1v, x2v, w, wv, idx, resid, hii, Cv
    real matrix R, invG, beta, H, sigma, cov, Xw, out
    real scalar n, k

    y = vec(y)
    x1 = vec(x1)
    x2 = vec(x2)
    C = vec(C)
    dist = _rdbw2d_loc_distance(x1, x2, ktype)
    w = vec(_rdbw2d_loc_weights(x1, x2, dist, h, kernel, ktype))
    idx = selectindex(w :> 0)
    n = length(idx)
    k = (p + 1) * (p + 2) / 2
    if (n < k) return(J(k + wantcov * k, k, .))

    yv = y[idx]
    x1v = x1[idx] :/ h
    x2v = x2[idx] :/ h
    wv = w[idx]
    if (hascluster) Cv = C[idx]
    else Cv = J(n, 1, .)

    R = _rdbw2d_loc_basis(x1v, x2v, p)
    invG = _rdbw2d_loc_inv(quadcross(R, wv, R))
    H = _rdbw2d_loc_H(h, p, 1)
    beta = H * invG * quadcross(R, wv, yv)
    if (wantcov == 0) return(beta)

    resid = yv - R * _rdbw2d_loc_H(h, p, 0) * beta
    if (vce == "hc1") {
        if (n <= k) return(J(k + k, k, .))
        if (!hascluster) resid = resid :* sqrt(n / (n - k))
    }
    else if (vce == "hc2" | vce == "hc3") {
        Xw = sqrt(wv) :* R
        hii = diagonal(Xw * invG * Xw')
        if (min(1 :- hii) <= 0) return(J(k + k, k, .))
        if (vce == "hc2") resid = resid :* sqrt(1 :/ (1 :- hii))
        else resid = resid :* (1 :/ (1 :- hii))
    }

    sigma = _rdbw2d_loc_meat(R, wv, resid, Cv, hascluster, h, vce)
    cov = invG' * sigma * invG
    out = J(2 * k, k, .)
    out[1..k, 1] = beta
    out[(k + 1)..(2 * k), 1..k] = cov
    return(out)
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_bwconst: Compute MSE bandwidth constants (bias+variance)
//
// Inputs:
//   y          : real colvector  - outcome variable
//   x1         : real colvector  - first running var (centered)
//   x2         : real colvector  - second running var (centered)
//   dist       : real colvector  - distances (recomputed internally)
//   C          : real colvector  - cluster variable
//   hascluster : real scalar     - 1 if clustered
//   p          : real scalar     - polynomial order for main fit
//   target     : real rowvector  - target functional vector
//   dn         : real scalar     - pilot bandwidth for variance
//   bn1        : real scalar     - pilot bw for bias (order p+1)
//   bn2        : real scalar     - pilot bw for regularization
//                                  (order p+2); missing = skip
//   vce        : string scalar   - VCE type
//   kernel     : string scalar   - kernel name
//   ktype      : string scalar   - kernel type (prod/rad)
//
// Output:
//   returns : real rowvector  - (B, V, Regb, Regv, n) where:
//     B    = leading bias constant
//     V    = variance constant (sandwich estimator)
//     Regb = regularization bias constant
//     Regv = regularization variance constant
//     n    = effective sample size within pilot bandwidth
//
// Notes:
//   - B = vecq * beta_{p+1}: bias times pilot coefficients
//   - V = target' * Sigma * target: pointwise variance
//   - Regb/Regv prevent near-zero-bias instability
//   - Returns missing components when pilot fits infeasible
// ---------------------------------------------------------------------------
real rowvector _rdbw2d_loc_bwconst(real colvector y, real colvector x1,
    real colvector x2, real colvector dist, real colvector C, real scalar hascluster,
    real scalar p, real rowvector target, real scalar dn, real scalar bn1,
    real scalar bn2, string scalar vce, string scalar kernel, string scalar ktype)
{
    real colvector yv, x1v, x2v, w, resid, hii, idx, Cv
    real matrix Raug, R, S, T, invG, H, beta, sigma, fit1, fit2
    real rowvector vecq, vect
    real scalar n, kp, kp1, kp2, B, V, Regb, Regv

    y = vec(y)
    x1 = vec(x1)
    x2 = vec(x2)
    C = vec(C)
    dist = _rdbw2d_loc_distance(x1, x2, ktype)

    w = vec(_rdbw2d_loc_weights(x1, x2, dist, dn, kernel, ktype))
    idx = selectindex(w :> 0)
    n = length(idx)
    kp = (p + 1) * (p + 2) / 2
    if (n < kp) return((., ., ., ., n))

    yv = y[idx]
    x1v = x1[idx] :/ dn
    x2v = x2[idx] :/ dn
    w = w[idx]
    if (hascluster) Cv = C[idx]
    else Cv = J(n, 1, .)

    kp1 = (p + 2) * (p + 3) / 2
    if (bn2 < .) {
        kp2 = (p + 3) * (p + 4) / 2
        Raug = _rdbw2d_loc_basis(x1v, x2v, p + 2)
    }
    else {
        kp2 = kp1
        Raug = _rdbw2d_loc_basis(x1v, x2v, p + 1)
    }
    R = Raug[, 1..kp]
    S = Raug[, (kp + 1)..kp1]
    invG = _rdbw2d_loc_inv(quadcross(R, w, R))

    vecq = (J(1, kp, 0), target * invG * quadcross(R, w, S))
    if (bn2 < .) {
        T = Raug[, (kp1 + 1)..kp2]
        vect = (J(1, kp1, 0), target * invG * quadcross(R, w, T))
    }

    H = _rdbw2d_loc_H(dn, p, 1)
    beta = H * invG * quadcross(R, w, yv)
    resid = yv - R * _rdbw2d_loc_H(dn, p, 0) * beta

    if (vce == "hc1") {
        if (n <= kp) return((., ., ., ., n))
        if (!hascluster) resid = resid :* sqrt(n / (n - kp))
    }
    else if (vce == "hc2" | vce == "hc3") {
        hii = diagonal((sqrt(w) :* R) * invG * (sqrt(w) :* R)')
        if (min(1 :- hii) <= 0) return((., ., ., ., n))
        if (vce == "hc2") resid = resid :* sqrt(1 :/ (1 :- hii))
        else resid = resid :* (1 :/ (1 :- hii))
    }

    sigma = _rdbw2d_loc_meat(R, w, resid, Cv, hascluster, dn, vce)
    V = (target * invG' * sigma * invG * target')[1, 1]

    fit1 = _rdbw2d_loc_lm_beta_cov(y, x1, x2, dist, C, hascluster, p + 1, bn1, vce, kernel, ktype, 1)
    if (missing(fit1[1, 1])) return((., V, ., ., n))
    beta = fit1[1..kp1, 1]
    B = (vecq * beta)[1, 1]
    Regv = (vecq * fit1[(kp1 + 1)..(2 * kp1), 1..kp1] * vecq')[1, 1] / (bn1^(2 + 2 * (p + 1)))

    Regb = 0
    if (bn2 < .) {
        fit2 = _rdbw2d_loc_lm_beta_cov(y, x1, x2, dist, C, hascluster, p + 2, bn2, vce, kernel, ktype, 0)
        if (missing(fit2[1, 1])) return((B, V, 0, Regv, n))
        Regb = (dn * vect * fit2[1..kp2, 1])[1, 1]
    }

    return((B, V, Regb, Regv, n))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_side: Local polynomial fit for one side at boundary point
//
// Inputs:
//   y       : real colvector  - outcome (full sample)
//   x1      : real colvector  - first running variable (full sample)
//   x2      : real colvector  - second running variable (full sample)
//   t       : real colvector  - treatment indicator (0/1)
//   side    : real scalar     - 0: control side; nonzero: treated
//   p       : real scalar     - polynomial order
//   kernel  : string scalar   - kernel name
//   ktype   : string scalar   - kernel type (prod/rad)
//   vce     : string scalar   - VCE type
//   h       : real scalar     - bandwidth
//   b1      : real scalar     - boundary point x1-coordinate
//   b2      : real scalar     - boundary point x2-coordinate
//   bwcheck : real scalar     - minimum obs constraint (unused)
//
// Output:
//   returns : real rowvector  - (n, beta0, var0, 0, 0) where:
//     n     = effective sample size
//     beta0 = intercept estimate (treatment effect component)
//     var0  = variance of intercept (sandwich [1,1])
//
// Notes:
//   - Centers coordinates at (b1, b2) before fitting
//   - Radial kernel uses sqrt(2)*h as effective bandwidth
//   - Returns all missing if insufficient observations
// ---------------------------------------------------------------------------
real rowvector _rdbw2d_loc_side(real colvector y, real colvector x1, real colvector x2,
    real colvector t, real scalar side, real scalar p, string scalar kernel,
    string scalar ktype, string scalar vce, real scalar h, real scalar b1, real scalar b2,
    real scalar bwcheck)
{
    real colvector xc1, xc2, w, fituse, resid, hii
    real matrix X, bread, ibread, beta, Xw, meat
    real scalar n, k, hr, out

    if (side == 0) fituse = (t :== 0)
    else fituse = (t :!= 0)

    xc1 = x1 :- b1
    xc2 = x2 :- b2
    if (ktype == "prod") {
        w = _rdbw2d_loc_kernel(xc1 :/ h, kernel) :* _rdbw2d_loc_kernel(xc2 :/ h, kernel) :/ (h^2)
    }
    else {
        hr = sqrt(2) * h
        w = _rdbw2d_loc_kernel(sqrt(xc1:^2 + xc2:^2) :/ hr, kernel) :/ (hr^2)
    }

    fituse = fituse :& (w :> 0)
    y = select(y, fituse)
    xc1 = select(xc1, fituse)
    xc2 = select(xc2, fituse)
    w = select(w, fituse)
    n = rows(y)
    k = (p + 1) * (p + 2) / 2
    if (n < k) return((., ., ., ., .))

    X = _rdbw2d_loc_basis(xc1, xc2, p)
    bread = quadcross(X, w, X)
    ibread = _rdbw2d_loc_inv(bread)
    beta = ibread * quadcross(X, w, y)
    resid = y - X * beta

    if (vce == "hc1") {
        resid = resid :* sqrt(n / (n - k))
    }
    else if (vce == "hc2" | vce == "hc3") {
        Xw = sqrt(w) :* X
        hii = diagonal(Xw * ibread * Xw')
        if (min(1 :- hii) <= 0) return((., ., ., ., .))
        if (vce == "hc2") resid = resid :* sqrt(1 :/ (1 :- hii))
        else resid = resid :* (1 :/ (1 :- hii))
    }

    meat = quadcross(X :* (w :* resid), X :* (w :* resid))
    out = (ibread * meat * ibread)[1,1]
    return((n, beta[1], out, 0, 0))
}

// ---------------------------------------------------------------------------
// _rdbw2d_loc_mata: Main dispatch for MSE-optimal bandwidth selection
//                   in 2D boundary RD
//
// Inputs:
//   yname        : string scalar  - Stata varname for outcome
//   x1name       : string scalar  - Stata varname for running var 1
//   x2name       : string scalar  - Stata varname for running var 2
//   tname        : string scalar  - Stata varname for treatment (0/1)
//   tousename    : string scalar  - Stata varname for sample marker
//   b1           : real scalar    - boundary point x1-coordinate
//   b2           : real scalar    - boundary point x2-coordinate
//   p            : real scalar    - polynomial order
//   kernel       : string scalar  - kernel name
//   ktype        : string scalar  - kernel type (prod/rad)
//   bwselect     : string scalar  - mserd/imserd/msetwo/imsetwo
//   method       : string scalar  - dpi or rot
//   vce          : string scalar  - VCE type (hc0/hc1/hc2/hc3)
//   bwcheck      : real scalar    - minimum observations constraint
//   masspoints   : string scalar  - check/adjust/off
//   scaleregul   : real scalar    - regularization scale (>= 0)
//   scalebiascrct: real scalar    - bias correction scale (>= 0)
//   stdflag      : real scalar    - 1: standardize vars; 0: no
//   cname        : string scalar  - cluster var name (empty=none)
//   targetname   : string scalar  - Stata matrix for target vector
//   derivsum     : real scalar    - sum of derivative orders
//   bwsname      : string scalar  - output matrix for bandwidths
//   constsname   : string scalar  - output matrix for MSE constants
//   massname     : string scalar  - output matrix for masspoints info
//
// Output:
//   void - stores results in Stata matrices: bwsname, constsname,
//          massname
//
// Notes:
//   - Full DPI/ROT algorithm: standardize, pilot bandwidth,
//     bias/variance constants, MSE-optimal formula, clamp, rescale
//   - mserd/imserd: common bw; msetwo/imsetwo: separate per side
//   - Final bandwidths: (b1, b2, h01, h02, h11, h12)
// ---------------------------------------------------------------------------
void _rdbw2d_loc_mata(
    string scalar yname,
    string scalar x1name,
    string scalar x2name,
    string scalar tname,
    string scalar tousename,
    real scalar b1,
    real scalar b2,
    real scalar p,
    string scalar kernel,
    string scalar ktype,
    string scalar bwselect,
    string scalar method,
    string scalar vce,
    real scalar bwcheck,
    string scalar masspoints,
    real scalar scaleregul,
    real scalar scalebiascrct,
    real scalar stdflag,
    string scalar cname,
    string scalar targetname,
    real scalar derivsum,
    string scalar bwsname,
    string scalar constsname,
    string scalar massname)
{
    real colvector y, x1, x2, t, C, idx0, idx1, x10, x20, x11, x21, d0, d1, d0check, d1check
    real colvector sort0, sort1
    real matrix unique0, unique1, uniqueall
    real scalar N, N0, N1, M, M0, M1, mass, sd1, sd2, dn, dn0, dn1, h0, h1, rawb1, rawb2
    real scalar minfit
    real scalar bwmin0, bwmin1, bwmax0, bwmax1, kp, g0, g1
    real rowvector e1, vecq0, vecq1, c0, c1, bnconst0, bnconst1, hnconst0, hnconst1
    real matrix bws, consts, massinfo, rankbasis0, rankbasis1

    st_view(y = ., ., yname, tousename)
    st_view(x1 = ., ., x1name, tousename)
    st_view(x2 = ., ., x2name, tousename)
    st_view(t = ., ., tname, tousename)
    if (cname != "") st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)

    rawb1 = b1
    rawb2 = b2
    N = rows(y)
    if (bwcheck > 0 & N < bwcheck) {
        errprintf("not enough observations to perform bandwidth calculations\n")
        _error(2001)
    }

    // stdvars: standardize running variables to unit variance before
    // bandwidth computation; final bandwidths are rescaled to original units
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
    N0 = sum(idx0)
    N1 = sum(idx1)
    minfit = (p + 2) * (p + 3) / 2
    if (N0 < minfit | N1 < minfit) {
        errprintf("each side must contain enough observations for bandwidth constant fits\n")
        _error(2001)
    }
    rankbasis0 = _rdbw2d_loc_basis(select(x1, idx0) :- b1, select(x2, idx0) :- b2, p + 1)
    rankbasis1 = _rdbw2d_loc_basis(select(x1, idx1) :- b1, select(x2, idx1) :- b2, p + 1)
    if (!_rdbw2d_loc_fullrank(quadcross(rankbasis0, rankbasis0)) |
        !_rdbw2d_loc_fullrank(quadcross(rankbasis1, rankbasis1))) {
        errprintf("bandwidth constant design is rank deficient on at least one side\n")
        _error(498)
    }

    if (masspoints != "off") {
        unique0 = uniqrows(select((x1, x2), idx0))
        unique1 = uniqrows(select((x1, x2), idx1))
        uniqueall = uniqrows((x1, x2))
        M0 = rows(unique0)
        M1 = rows(unique1)
        M = rows(uniqueall)
        mass = 1 - M / N
    }
    else {
        M = N
        M0 = N0
        M1 = N1
        mass = 0
    }
    massinfo = (M, M0, M1, mass)

    d0 = sqrt((select(x1, idx0) :- b1):^2 + (select(x2, idx0) :- b2):^2)
    d1 = sqrt((select(x1, idx1) :- b1):^2 + (select(x2, idx1) :- b2):^2)
    if (ktype == "prod") {
        d0 = abs(select(x1, idx0) :- b1)
        d1 = abs(select(x1, idx1) :- b1)
        d0 = (d0 :> abs(select(x2, idx0) :- b2)) :* d0 :+ (d0 :<= abs(select(x2, idx0) :- b2)) :* abs(select(x2, idx0) :- b2)
        d1 = (d1 :> abs(select(x2, idx1) :- b2)) :* d1 :+ (d1 :<= abs(select(x2, idx1) :- b2)) :* abs(select(x2, idx1) :- b2)
    }
    d0check = d0
    d1check = d1
    if (masspoints == "adjust") {
        d0check = sqrt((unique0[, 1] :- b1):^2 + (unique0[, 2] :- b2):^2)
        d1check = sqrt((unique1[, 1] :- b1):^2 + (unique1[, 2] :- b2):^2)
        if (ktype == "prod") {
            d0check = (abs(unique0[, 1] :- b1) :> abs(unique0[, 2] :- b2)) :*
                abs(unique0[, 1] :- b1) :+
                (abs(unique0[, 1] :- b1) :<= abs(unique0[, 2] :- b2)) :*
                abs(unique0[, 2] :- b2)
            d1check = (abs(unique1[, 1] :- b1) :> abs(unique1[, 2] :- b2)) :*
                abs(unique1[, 1] :- b1) :+
                (abs(unique1[, 1] :- b1) :<= abs(unique1[, 2] :- b2)) :*
                abs(unique1[, 2] :- b2)
        }
    }
    dn = _rdbw2d_loc_rot(x1, x2, kernel)
    if (!(dn > 0)) {
        errprintf("rule-of-thumb bandwidth could not be calculated\n")
        _error(498)
    }

    kp = (p + 1) * (p + 2) / 2
    e1 = st_matrix(targetname)
    if (cols(e1) != kp) {
        errprintf("deriv()/tangvec() target dimension mismatch\n")
        _error(498)
    }

    bwmin0 = .
    bwmin1 = .
    bwmax0 = .
    bwmax1 = .
    if (bwcheck > 0) {
        sort0 = sort(d0check, 1)
        sort1 = sort(d1check, 1)
        bwmin0 = sort0[min((bwcheck, rows(sort0)))]
        bwmin1 = sort1[min((bwcheck, rows(sort1)))]
        bwmax0 = sort0[rows(sort0)]
        bwmax1 = sort1[rows(sort1)]
    }
    dn0 = dn
    dn1 = dn
    if (bwcheck > 0) {
        dn0 = _rdbw2d_loc_clamp(dn0, bwmin0, bwmax0)
        dn1 = _rdbw2d_loc_clamp(dn1, bwmin1, bwmax1)
    }

    h0 = _rdbw2d_loc_median(d0)
    h1 = _rdbw2d_loc_median(d1)
    // DPI method: use pilot bandwidth to estimate bias and variance constants,
    // then compute h* = (C_V / C_B)^{1/(2p+4)} where C_V is the variance
    // constant and C_B is the squared bias constant (plus regularization).
    if (method == "dpi") {
        vecq0 = _rdbw2d_loc_coeff(select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0, p, e1, dn0, kernel, ktype)
        vecq1 = _rdbw2d_loc_coeff(select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1, p, e1, dn1, kernel, ktype)
        if (bwcheck <= 0 & missing(vecq0[1])) {
            dn0 = max((dn0, max(d0) * (1 + 1e-8)))
            vecq0 = _rdbw2d_loc_coeff(select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0, p, e1, dn0, kernel, ktype)
        }
        if (bwcheck <= 0 & missing(vecq1[1])) {
            dn1 = max((dn1, max(d1) * (1 + 1e-8)))
            vecq1 = _rdbw2d_loc_coeff(select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1, p, e1, dn1, kernel, ktype)
        }

        bnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
            select(C, idx0), cname != "", p + 1, vecq0, dn0, h0, ., vce, kernel, ktype)
        bnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
            select(C, idx1), cname != "", p + 1, vecq1, dn1, h1, ., vce, kernel, ktype)
        if (missing(bnconst0[1])) {
            h0 = max((h0, dn0, max(d0) * (1 + 1e-8)))
            bnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
                select(C, idx0), cname != "", p + 1, vecq0, dn0, h0, ., vce, kernel, ktype)
        }
        if (missing(bnconst1[1])) {
            h1 = max((h1, dn1, max(d1) * (1 + 1e-8)))
            bnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
                select(C, idx1), cname != "", p + 1, vecq1, dn1, h1, ., vce, kernel, ktype)
        }
        // MSE-optimal bandwidth formula: h* = ((d+2*s)*V / ((2*s+d-2*s) * B^2))^{1/(2p+6)}
        // where d=dimension, s=derivative order for the pilot step
        if (!(missing(bnconst0[1]) | missing(bnconst1[1]))) {
            h0 = ((2 + 2 * (p + 1)) * bnconst0[2] / ((2 * (p + 1) + 2 - 2 * (p + 1)) * (bnconst0[1]^2 + scaleregul * bnconst0[4])))^(1 / (2 * p + 6))
            h1 = ((2 + 2 * (p + 1)) * bnconst1[2] / ((2 * (p + 1) + 2 - 2 * (p + 1)) * (bnconst1[1]^2 + scaleregul * bnconst1[4])))^(1 / (2 * p + 6))
            if (!(h0 > 0) | missing(h0)) h0 = dn0
            if (!(h1 > 0) | missing(h1)) h1 = dn1
        }
        if (bwcheck > 0) {
            h0 = _rdbw2d_loc_clamp(h0, bwmin0, bwmax0)
            h1 = _rdbw2d_loc_clamp(h1, bwmin1, bwmax1)
        }
    }

    // Final MSE-optimal bandwidth: h* = ((d+2*v)*V / ((2p+d-2v) * B^2))^{1/(2p+4)}
    // where d=dimension=2, v=derivsum, V=variance constant, B=bias constant
    // scaleregul adds regularization to prevent division by near-zero bias
    if (bwselect == "mserd" | bwselect == "imserd") {
        hnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
            select(C, idx0), cname != "", p, e1, dn0, h0, _rdbw2d_loc_median(d0), vce, kernel, ktype)
        hnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
            select(C, idx1), cname != "", p, e1, dn1, h1, _rdbw2d_loc_median(d1), vce, kernel, ktype)
        if (bwcheck <= 0 & missing(hnconst0[1])) {
            h0 = max((h0, dn0, max(d0) * (1 + 1e-8)))
            hnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
                select(C, idx0), cname != "", p, e1, dn0, h0, _rdbw2d_loc_median(d0), vce, kernel, ktype)
        }
        if (bwcheck <= 0 & missing(hnconst1[1])) {
            h1 = max((h1, dn1, max(d1) * (1 + 1e-8)))
            hnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
                select(C, idx1), cname != "", p, e1, dn1, h1, _rdbw2d_loc_median(d1), vce, kernel, ktype)
        }
        if (missing(hnconst0[1]) | missing(hnconst1[1])) {
            errprintf("final bandwidth constants could not be calculated\n")
            _error(498)
        }
        h0 = ((2 + 2 * derivsum) * (hnconst0[2] + hnconst1[2]) / ((2 * p + 2 - 2 * derivsum) * ((hnconst0[1] + scalebiascrct * hnconst0[3] - hnconst1[1] - scalebiascrct * hnconst1[3])^2 + scaleregul * hnconst0[4] + scaleregul * hnconst1[4])))^(1 / (2 * p + 4))
        if (!(h0 > 0) | missing(h0)) h0 = max((dn0, dn1))
        h1 = h0
        if (bwcheck > 0) {
            h0 = _rdbw2d_loc_clamp(h0, max((bwmin0, bwmin1)), max((bwmax0, bwmax1)))
            h1 = h0
        }
    }
    else {
        hnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
            select(C, idx0), cname != "", p, e1, dn0, h0, _rdbw2d_loc_median(d0), vce, kernel, ktype)
        hnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
            select(C, idx1), cname != "", p, e1, dn1, h1, _rdbw2d_loc_median(d1), vce, kernel, ktype)
        if (bwcheck <= 0 & missing(hnconst0[1])) {
            h0 = max((h0, dn0, max(d0) * (1 + 1e-8)))
            hnconst0 = _rdbw2d_loc_bwconst(select(y, idx0), select(x1, idx0) :- b1, select(x2, idx0) :- b2, d0,
                select(C, idx0), cname != "", p, e1, dn0, h0, _rdbw2d_loc_median(d0), vce, kernel, ktype)
        }
        if (bwcheck <= 0 & missing(hnconst1[1])) {
            h1 = max((h1, dn1, max(d1) * (1 + 1e-8)))
            hnconst1 = _rdbw2d_loc_bwconst(select(y, idx1), select(x1, idx1) :- b1, select(x2, idx1) :- b2, d1,
                select(C, idx1), cname != "", p, e1, dn1, h1, _rdbw2d_loc_median(d1), vce, kernel, ktype)
        }
        if (missing(hnconst0[1]) | missing(hnconst1[1])) {
            errprintf("final bandwidth constants could not be calculated\n")
            _error(498)
        }
        h0 = ((2 + 2 * derivsum) * hnconst0[2] / ((2 * p + 2 - 2 * derivsum) * ((hnconst0[1] + scalebiascrct * hnconst0[3])^2 + scaleregul * hnconst0[4])))^(1 / (2 * p + 4))
        h1 = ((2 + 2 * derivsum) * hnconst1[2] / ((2 * p + 2 - 2 * derivsum) * ((hnconst1[1] + scalebiascrct * hnconst1[3])^2 + scaleregul * hnconst1[4])))^(1 / (2 * p + 4))
        if (!(h0 > 0) | missing(h0)) h0 = dn0
        if (!(h1 > 0) | missing(h1)) h1 = dn1
        if (bwcheck > 0) {
            h0 = _rdbw2d_loc_clamp(h0, bwmin0, bwmax0)
            h1 = _rdbw2d_loc_clamp(h1, bwmin1, bwmax1)
        }
    }

    bws = (rawb1, rawb2, h0 * sd1, h0 * sd2, h1 * sd1, h1 * sd2)
    consts = (hnconst0[5], hnconst1[5], hnconst0[1], hnconst1[1], hnconst0[2], hnconst1[2], hnconst0[3], hnconst1[3], hnconst0[4], hnconst1[4])

    st_matrix(bwsname, bws)
    st_matrix(constsname, consts)
    st_matrix(massname, massinfo)
}
end
