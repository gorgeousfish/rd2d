*! version 1.1.0 24jun2026
program define rdbw2d_dist, rclass
    version 16.0
    if `"`1'"' == "version" | `"`1'"' == ",version" {
        di as txt "rdbw2d_dist version 1.1.0 (24 June 2026)"
        exit
    }
    syntax varlist(min=2 numeric) [if] [in], ///
        [P(integer 1) KINK(string) KERnel(string) BWSELect(string) VCE(string) ///
         BWCHeck(string) MASSPoints(string) SCALEregul(real 1) CQT(real 0.5) ///
         CLuster(varname) ///
         MP(string) SR(real -1) BWS(string)]

    gettoken yvar dvars : varlist
    local neval : word count `dvars'
    if (`neval' == 0) {
        di as err "at least one distance variable is required after the outcome variable"
        exit 198
    }
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
        di as err "p() must be a nonnegative integer (0, 1, 2, ...); got `p'"
        exit 198
    }
    if (`scaleregul' >= . | `scaleregul' < 0) {
        di as err "scaleregul() must be finite and nonnegative"
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
        di as err "kernel() must be uniform, triangular, epanechnikov, or gaussian (aliases: uni, tri, epa, gau); got '`kernel''"
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

    local kink = lower("`kink'")
    if ("`kink'" == "") local kink "off"
    if !inlist("`kink'", "off", "on") {
        di as err "kink() must be 'on' or 'off'; kink(on) uses undersmoothing instead of RBC. Got '`kink''"
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
                exit 198
            }
        }
    }

    local firstd : word 1 of `dvars'
    quietly count if `touse' & `firstd' < 0
    local N0 = r(N)
    quietly count if `touse' & `firstd' >= 0
    local N1 = r(N)
    if (`N0' < `p' + 2 | `N1' < `p' + 2) {
        local failside = cond(`N0' < `p' + 2, "control", "treated")
        local failnh = cond(`N0' < `p' + 2, `N0', `N1')
        di as err "insufficient observations on `failside' side (Nh=`failnh', min=`=`p'+2'); consider reducing p()"
        exit 2001
    }

    tempname bws mseconsts massinfo
    mata: _rdbw2d_dist_mata("`yvar'", "`dvars'", "`touse'", "`clusterwork'", ///
        `p', "`kernel'", "`bwselect'", "`kink'", "`vce'", `bwcheck', ///
        `scaleregul', `cqt', "`masspoints'", "`bws'", "`mseconsts'", "`massinfo'")

    matrix colnames `bws' = b1 b2 h0 h1 Nh0 Nh1
    matrix colnames `mseconsts' = h0 h1 bias0 bias1 var0 var1 reg0 reg1 Nh0 Nh1 bwmin0 bwmin1 bwmax0 bwmax1
    // massinfo matrix columns: [M, M0, M1, mass_ratio]
    // M = total unique distance values, M0/M1 = unique per side
    // mass_ratio = 1 - M/N; values >= 0.2 trigger warnings
    matrix colnames `massinfo' = M M0 M1 mass
    matrix rownames `bws' = `dvars'
    matrix rownames `mseconsts' = `dvars'
    matrix rownames `massinfo' = `dvars'

    _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')

    // ===================================================================
    // TABLE DISPLAY SECTION
    // Layout modes: ultra (<50), compact (<63), narrow (<79), normal (>=79)
    // Columns: Point, h0, h1, Nh0, Nh1
    // ===================================================================
    * --- Table layout parameters ---
    local line_width = min(79, c(linesize))
    local hline_rule `"di as txt \"{hline `line_width'}\""'
    local ultra_table = c(linesize) < 50
    local compact_table = c(linesize) < 63
    local narrow_table = c(linesize) < 79
    di as txt _newline "Distance bandwidth selection"
    di as txt "{hline `line_width'}"
    if (`ultra_table' | `compact_table') {
        di as txt "  Eval points: " as res %9.0f `neval'
        di as txt "  Observations: " as res %9.0f `N'
        di as txt "  Selector: " as res "`bwselect'"
        di as txt "  Kernel: " as res "`kernel'" as txt "  VCE: " as res "`vce'"
        di as txt "  Kink: " as res "`kink'"
    }
    else {
        di as txt "  Evaluation points: " as res %9.0f `neval' ///
            as txt "    Observations: " as res %9.0f `N'
        di as txt "  Selector: " as res "`bwselect'" ///
            as txt "    Kernel: " as res "`kernel'" ///
            as txt "    Kink: " as res "`kink'" ///
            as txt "    VCE: " as res "`vce'"
    }
    di as txt "{hline `line_width'}"
    if (`ultra_table') {
        di as txt %4s "Pt" _col(6) %6s "h0" _col(13) %6s "h1" ///
            _col(20) %6s "Nh0" _col(27) %6s "Nh1"
    }
    else if (`compact_table') {
        di as txt %6s "Point" _col(8) %6s "h0" _col(15) %6s "h1" ///
            _col(22) %6s "Nh0" _col(29) %6s "Nh1"
    }
    else if (`narrow_table') {
        di as txt %8s "Point" _col(10) %8s "h0" _col(19) %8s "h1" ///
            _col(28) %8s "Nh0" _col(37) %8s "Nh1"
    }
    else {
        di as txt %10s "Point" _col(12) %10s "h0" _col(24) %10s "h1" ///
            _col(36) %10s "Nh0" _col(48) %10s "Nh1"
    }
    di as txt "{hline `line_width'}"
    forvalues j = 1/`neval' {
        local rname : word `j' of `dvars'
        local suffix : display %02.0f `j'
        local label_width = cond(`ultra_table', 4, cond(`compact_table', 6, cond(`narrow_table', 8, 10)))
        local prefix_len = max(0, `label_width' - strlen("`suffix'") - 1)
        local dname = cond(strlen("`rname'") > `label_width', substr("`rname'", 1, `prefix_len') + "~`suffix'", "`rname'")
        local dh0 = el(`bws', `j', 3)
        local dh1 = el(`bws', `j', 4)
        local dnh0 = el(`bws', `j', 5)
        local dnh1 = el(`bws', `j', 6)
        if (`ultra_table') {
            di as res %4s "`dname'" _col(6) %6.3g `dh0' ///
                _col(13) %6.3g `dh1' _col(20) %6.0g `dnh0' ///
                _col(27) %6.0g `dnh1'
        }
        else if (`compact_table') {
            di as res %6s "`dname'" _col(8) %6.4g `dh0' ///
                _col(15) %6.4g `dh1' _col(22) %6.0g `dnh0' ///
                _col(29) %6.0g `dnh1'
        }
        else if (`narrow_table') {
            di as res %8s "`dname'" _col(10) %8.4g `dh0' ///
                _col(19) %8.4g `dh1' _col(28) %8.0g `dnh0' ///
                _col(37) %8.0g `dnh1'
        }
        else {
            di as res %10s "`dname'" _col(12) %10.4g `dh0' ///
                _col(24) %10.4g `dh1' _col(36) %10.0g `dnh0' ///
                _col(48) %10.0g `dnh1'
        }
    }
    di as txt "{hline `line_width'}"

    return matrix bws = `bws'
    return matrix mseconsts = `mseconsts'
    return matrix masspoints = `massinfo'
    return scalar N = `N'
    return scalar N0 = `N0'
    return scalar N1 = `N1'
    return scalar neval = `neval'
    return scalar p = `p'
    return scalar bwcheck = `bwcheck'
    return scalar scaleregul = `scaleregul'
    return scalar cqt = `cqt'
    return local kernel "`kernel'"
    return local bwselect "`bwselect'"
    return local kink "`kink'"
    return local vce "`vce'"
    return local masspoints_opt "`masspoints'"
    return local version "1.1.0"
end

mata:
real colvector _rdbw2d_dist_kernel(real colvector u, string scalar kernel)
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

real matrix _rdbw2d_dist_basis(real colvector x, real scalar p)
{
    real scalar j
    real matrix X

    X = J(rows(x), p + 1, 1)
    for (j = 2; j <= p + 1; j++) {
        X[, j] = x:^(j - 1)
    }
    return(X)
}

real scalar _rdbw2d_dist_fullrank(real matrix A)
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

real scalar _rdbw2d_dist_quantile(real colvector x, real scalar q)
{
    real scalar n, h, lo, hi, g
    real colvector sx

    sx = sort(x, 1)
    n = rows(sx)
    if (n == 0) return(.)
    if (n == 1) return(sx[1])
    h = (n - 1) * q + 1
    lo = floor(h)
    hi = ceil(h)
    g = h - lo
    if (lo == hi) return(sx[lo])
    return((1 - g) * sx[lo] + g * sx[hi])
}

real scalar _rdbw2d_dist_rot(real colvector x, string scalar kernel)
{
    real scalar mu2K, l2K, N, varhat, traceconst

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

    N = rows(x)
    varhat = mean(x:^2) / 2
    if (varhat <= 0) return(.)
    traceconst = 1 / (2 * pi() * varhat^3)
    return((2 * l2K / (N * mu2K * traceconst))^(1 / 6))
}

real rowvector _rdbw2d_dist_ols_deriv(real colvector y, real colvector x, real scalar p)
{
    real matrix X, XXi, V
    real colvector beta, resid
    real scalar n, k, s2, last

    X = _rdbw2d_dist_basis(x, p + 1)
    n = rows(X)
    k = cols(X)
    if (n <= k) {
        errprintf("not enough observations for derivative pilot fit\n")
        _error(2001)
    }
    if (!_rdbw2d_dist_fullrank(quadcross(X, X))) {
        errprintf("derivative pilot design is rank deficient\n")
        _error(498)
    }
    XXi = invsym(quadcross(X, X))
    beta = XXi * quadcross(X, y)
    resid = y - X * beta
    s2 = quadcross(resid, resid) / (n - k)
    V = s2 * XXi
    last = rows(beta)
    return((beta[last], sqrt(V[last, last])))
}

real matrix _rdbw2d_dist_meat(real matrix WR, real colvector resid, real colvector C,
    real scalar hascluster, real scalar h, string scalar vce)
{
    real scalar n, k, i, g, factor
    real matrix M, scores
    real colvector clusters, idx

    n = rows(WR)
    k = cols(WR)
    if (!hascluster) {
        return(quadcross(resid :* WR, resid :* WR) * h^2)
    }

    clusters = uniqrows(sort(C, 1))
    g = rows(clusters)
    if (g <= 1) {
        errprintf("cluster() must contain at least two clusters\n")
        _error(2001)
    }
    M = J(k, k, 0)
    for (i = 1; i <= g; i++) {
        idx = (C :== clusters[i])
        scores = quadcross(select(WR, idx), select(resid, idx))
        M = M + scores * scores'
    }
    factor = 1
    if (vce == "hc1") factor = ((n - 1) / (n - k)) * (g / (g - 1))
    return(M * h^2 * factor)
}

real rowvector _rdbw2d_dist_side_constants(real colvector y, real colvector x,
    real colvector C, real scalar hascluster, real scalar p, string scalar kernel,
    string scalar vce, real scalar dn, real scalar dn_base, real scalar bwcheck,
    real scalar cqt, string scalar masspoints)
{
    real scalar qcut, bwmin, bwmax, eN, k, coeff, bias, reg, variance
    real rowvector deriv
    real colvector pilot_idx, xpilot, ypilot, u, w, inband, xs, ys, ws, resid, leverage, vcew, Cs
    real matrix R, gammai, sigmahalf, sigma, pmatrix, sqrtwR

    if (rows(x) < p + 2) {
        errprintf("each side must contain at least p()+2 observations\n")
        _error(2001)
    }
    // Pilot fit: use cqt-quantile of |D| to select observations for the
    // derivative estimate. This limits the pilot sample to the fraction closest
    // to the boundary, reducing contamination from far-away observations.
    qcut = _rdbw2d_dist_quantile(x, cqt)
    pilot_idx = (x :<= qcut)
    xpilot = select(x, pilot_idx)
    ypilot = select(y, pilot_idx)
    deriv = _rdbw2d_dist_ols_deriv(ypilot, xpilot, p)

    bwmin = .
    bwmax = .
    if (bwcheck > 0) {
        xs = sort(x, 1)
        if (masspoints == "adjust") xs = uniqrows(xs)
        bwmin = xs[min((bwcheck, rows(xs)))]
        bwmax = xs[rows(xs)]
        dn = max((dn, bwmin))
        dn = min((dn, bwmax))
    }

    u = x :/ dn
    w = _rdbw2d_dist_kernel(u, kernel) :/ (dn^2)
    inband = (w :> 0)
    eN = sum(inband)
    if (eN < p + 1) {
        errprintf("not enough effective observations inside the bandwidth\n")
        _error(2001)
    }

    xs = select(x, inband)
    ys = select(y, inband)
    ws = select(w, inband)
    Cs = select(C, inband)
    R = _rdbw2d_dist_basis(xs :/ dn, p)
    if (!_rdbw2d_dist_fullrank(quadcross(R, ws, R))) {
        errprintf("bandwidth constant design is rank deficient\n")
        _error(498)
    }
    gammai = invsym(quadcross(R, ws, R))
    resid = ys - _rdbw2d_dist_basis(xs, p + 1) * invsym(quadcross(_rdbw2d_dist_basis(xpilot, p + 1), _rdbw2d_dist_basis(xpilot, p + 1))) * quadcross(_rdbw2d_dist_basis(xpilot, p + 1), ypilot)

    k = p + 1
    sqrtwR = R :* sqrt(ws)
    if (vce == "hc0") {
        vcew = J(rows(resid), 1, 1)
    }
    else if (vce == "hc1") {
        if (eN <= k) {
            errprintf("not enough effective observations for hc1 degrees-of-freedom correction\n")
            _error(2001)
        }
        vcew = J(rows(resid), 1, (hascluster ? 1 : sqrt(eN / (eN - k))))
    }
    else {
        if (eN <= k) {
            errprintf("not enough effective observations for hc2/hc3 leverage correction\n")
            _error(2001)
        }
        leverage = rowsum((sqrtwR * gammai) :* sqrtwR)
        if (min(1 :- leverage) <= 0) {
            errprintf("leverage is too high for hc2/hc3 variance calculation\n")
            _error(2001)
        }
        if (vce == "hc2") vcew = sqrt(1 :/ (1 :- leverage))
        else vcew = 1 :/ (1 :- leverage)
    }
    resid = resid :* vcew
    sigmahalf = R :* ws
    sigma = _rdbw2d_dist_meat(sigmahalf, resid, Cs, hascluster, dn, vce)

    // MSE constant computation: bias = e1' * Gamma^{-1} * P * beta_{p+1}
    // where P projects the (p+1)-th order monomial onto the p-th order basis
    pmatrix = quadcross(R, ((xs :/ dn_base):^(p + 1)) :* ws)
    coeff = (gammai * pmatrix)[1]
    bias = coeff * deriv[1]
    reg = coeff^2 * deriv[2]^2
    variance = (gammai * sigma * gammai)[1, 1]

    return((bias, variance, reg, eN, bwmin, bwmax, dn))
}

void _rdbw2d_dist_mata(string scalar yname, string scalar dnames, string scalar tousename,
    string scalar cname, real scalar p, string scalar kernel, string scalar bwselect,
    string scalar kink, string scalar vce, real scalar bwcheck, real scalar scaleregul,
    real scalar cqt, string scalar masspoints, string scalar bwsname,
    string scalar constsname, string scalar massname)
{
    real scalar N, neval, j, dn, N0, N1, M, M0, M1, mass, h0, h1
    real colvector y, sdist, x, d1, C, idx0, idx1, uniquevals
    real matrix D, bws, consts, massinfo
    real rowvector c0, c1
    string rowvector dvarlist
    real scalar V, B, V0, V1, B0, B1, H

    dvarlist = tokens(dnames)
    st_view(y = ., ., yname, tousename)
    st_view(D = ., ., dvarlist, tousename)
    if (cname != "") st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)

    N = rows(y)
    neval = cols(D)
    d1 = D[, 1]
    N0 = sum(d1 :< 0)
    N1 = sum(d1 :>= 0)
    bws = J(neval, 6, .)
    consts = J(neval, 14, .)
    massinfo = J(neval, 4, .)

    for (j = 1; j <= neval; j++) {
        // signed distance: D<0 = control side, D>=0 = treated side
        sdist = D[, j]
        x = abs(sdist)
        idx0 = (sdist :< 0)
        idx1 = (sdist :>= 0)
        N0 = sum(idx0)
        N1 = sum(idx1)
        if (sum(idx0) < p + 2 | sum(idx1) < p + 2) {
            errprintf("each distance column must contain at least p()+2 observations on each side\n")
            _error(2001)
        }

        // massinfo: [M, M0, M1, mass_ratio] per distance column
        M = N
        M0 = N0
        M1 = N1
        mass = 0
        if (masspoints != "off") {
            uniquevals = uniqrows(sort(sdist, 1))
            M0 = sum(uniquevals :< 0)
            M1 = sum(uniquevals :>= 0)
            M = M0 + M1
            mass = 1 - M / N
        }
        massinfo[j, ] = (M, M0, M1, mass)

        dn = _rdbw2d_dist_rot(x, kernel)
        if (dn == . | dn <= 0) {
            errprintf("rule-of-thumb bandwidth could not be calculated\n")
            _error(498)
        }

        c0 = _rdbw2d_dist_side_constants(select(y, idx0), select(x, idx0), select(C, idx0),
            cname != "", p, kernel, vce, dn, dn, bwcheck, cqt, masspoints)
        c1 = _rdbw2d_dist_side_constants(select(y, idx1), select(x, idx1), select(C, idx1),
            cname != "", p, kernel, vce, dn, dn, bwcheck, cqt, masspoints)

        consts[j, 3] = c0[1]
        consts[j, 4] = c1[1]
        consts[j, 5] = c0[2]
        consts[j, 6] = c1[2]
        consts[j, 7] = c0[3]
        consts[j, 8] = c1[3]
        consts[j, 9] = c0[4]
        consts[j, 10] = c1[4]
        consts[j, 11] = c0[5]
        consts[j, 12] = c1[5]
        consts[j, 13] = c0[6]
        consts[j, 14] = c1[6]
    }

    // MSE-optimal bandwidth: h* = (2*V / ((2p+2) * B^2))^{1/(2p+4)}
    // mserd: common h for both sides, minimizing combined MSE
    if (bwselect == "mserd") {
        for (j = 1; j <= neval; j++) {
            H = (2 * (consts[j, 5] + consts[j, 6]) /
                ((2 * p + 2) * ((consts[j, 3] - consts[j, 4])^2 +
                scaleregul * consts[j, 7] + scaleregul * consts[j, 8])))^(1 / (2 * p + 4))
            if (bwcheck > 0) H = min((max((H, consts[j, 11], consts[j, 12])), consts[j, 13], consts[j, 14]))
            consts[j, 1] = H
            consts[j, 2] = H
        }
    }
    // msetwo: separate h0, h1 for control and treated sides
    else if (bwselect == "msetwo") {
        for (j = 1; j <= neval; j++) {
            h0 = (2 * consts[j, 5] / ((2 * p + 2) * (consts[j, 3]^2 + scaleregul * consts[j, 7])))^(1 / (2 * p + 4))
            h1 = (2 * consts[j, 6] / ((2 * p + 2) * (consts[j, 4]^2 + scaleregul * consts[j, 8])))^(1 / (2 * p + 4))
            if (bwcheck > 0) {
                h0 = min((max((h0, consts[j, 11])), consts[j, 13]))
                h1 = min((max((h1, consts[j, 12])), consts[j, 14]))
            }
            consts[j, 1] = h0
            consts[j, 2] = h1
        }
    }
    // imserd: integrated MSE-optimal, averages constants across eval points
    else if (bwselect == "imserd") {
        V = mean(consts[, 5]) + mean(consts[, 6])
        B = mean((consts[, 3] - consts[, 4]):^2 + scaleregul * consts[, 7] + scaleregul * consts[, 8])
        H = (2 * V / ((2 * p + 2) * B))^(1 / (2 * p + 4))
        for (j = 1; j <= neval; j++) {
            consts[j, 1] = H
            consts[j, 2] = H
            if (bwcheck > 0) {
                consts[j, 1] = min((max((consts[j, 1], consts[j, 11], consts[j, 12])), consts[j, 13], consts[j, 14]))
                consts[j, 2] = consts[j, 1]
            }
        }
    }
    // imsetwo: integrated MSE-optimal with separate bandwidths per side
    else {
        V0 = mean(consts[, 5])
        V1 = mean(consts[, 6])
        B0 = mean(consts[, 3]:^2 + scaleregul * consts[, 7])
        B1 = mean(consts[, 4]:^2 + scaleregul * consts[, 8])
        h0 = (2 * V0 / ((2 * p + 2) * B0))^(1 / (2 * p + 4))
        h1 = (2 * V1 / ((2 * p + 2) * B1))^(1 / (2 * p + 4))
        for (j = 1; j <= neval; j++) {
            consts[j, 1] = h0
            consts[j, 2] = h1
            if (bwcheck > 0) {
                consts[j, 1] = min((max((consts[j, 1], consts[j, 11])), consts[j, 13]))
                consts[j, 2] = min((max((consts[j, 2], consts[j, 12])), consts[j, 14]))
            }
        }
    }

    // kink(on): boundary kink causes distance-based conditional expectation
    // to lose differentiability, degrading bias from h^{p+1} to h.
    // Solution: use undersmoothing (q=p) instead of RBC (q=p+1).
    // Minimax convergence rate at kink: n^{-1/4}.
    // Bandwidth adjustment: rescale from MSE-optimal rate n^{-1/(2p+4)} to
    // kink-optimal rate n^{-1/4} using the ratio N^{-1/4} / N^{-1/(2p+4)}.
    if (kink == "on") {
        for (j = 1; j <= neval; j++) {
            if (bwselect == "mserd" | bwselect == "imserd") {
                consts[j, 1] = consts[j, 1] * massinfo[j, 1]^(-1/4) / massinfo[j, 1]^(-1/(2 * p + 4))
                consts[j, 2] = consts[j, 2] * massinfo[j, 1]^(-1/4) / massinfo[j, 1]^(-1/(2 * p + 4))
            }
            else {
                consts[j, 1] = consts[j, 1] * massinfo[j, 2]^(-1/4) / massinfo[j, 2]^(-1/(2 * p + 4))
                consts[j, 2] = consts[j, 2] * massinfo[j, 3]^(-1/4) / massinfo[j, 3]^(-1/(2 * p + 4))
            }
        }
    }

    for (j = 1; j <= neval; j++) {
        if (!(consts[j, 1] > 0) | !(consts[j, 2] > 0) | missing(consts[j, 1]) | missing(consts[j, 2])) {
            errprintf("final bandwidth could not be calculated\n")
            _error(498)
        }
    }

    bws[, 3] = consts[, 1]
    bws[, 4] = consts[, 2]
    bws[, 5] = consts[, 9]
    bws[, 6] = consts[, 10]
    st_matrix(bwsname, bws)
    st_matrix(constsname, consts)
    st_matrix(massname, massinfo)
}
end
