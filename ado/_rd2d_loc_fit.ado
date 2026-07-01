*! version 1.1.0 24jun2026
program define _rd2d_loc_fit, rclass
    version 16.0
    syntax varlist(min=4 max=4 numeric) [if] [in], AT(string asis) HX(real) ///
        [HY(string) P(integer 1) SIDE(string) KERnel(string) KTYPE(string) STDVars VCE(string) ///
         CLuster(varname)]

    gettoken yvar rest : varlist
    gettoken x1var rest : rest
    gettoken x2var dvar : rest
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

    if (`p' < 0) {
        di as err "p() must be a nonnegative integer"
        exit 198
    }
    if (`hx' >= . | `hx' <= 0) {
        di as err "hx() must be a finite positive number"
        exit 198
    }
    local has_hy = ("`hy'" != "")
    local hyval = .
    if (`has_hy') {
        capture confirm number `hy'
        if (_rc) {
            di as err "hy() must be numeric"
            exit 198
        }
        local hyval = real("`hy'")
        if (`hyval' >= . | `hyval' <= 0) {
            di as err "hy() must be a finite positive number"
            exit 198
        }
    }
    else {
        local hyval = `hx'
    }

    local side = lower("`side'")
    if ("`side'" == "") local side "control"
    if !inlist("`side'", "control", "treated") {
        di as err "side() must be control or treated"
        exit 198
    }

    local kernel = lower("`kernel'")
    if ("`kernel'" == "") local kernel "triangular"
    if inlist("`kernel'", "uni", "unif") local kernel "uniform"
    if inlist("`kernel'", "tri", "triag") local kernel "triangular"
    if inlist("`kernel'", "epa", "epan") local kernel "epanechnikov"
    if inlist("`kernel'", "gau") local kernel "gaussian"
    if !inlist("`kernel'", "uniform", "triangular", "epanechnikov", "gaussian") {
        di as err "kernel() must be uniform, triangular, epanechnikov, or gaussian"
        exit 198
    }

    local ktype = lower("`ktype'")
    if ("`ktype'" == "") local ktype "prod"
    if !inlist("`ktype'", "prod", "rad") {
        di as err "ktype() must be prod or rad"
        exit 198
    }

    local vce = lower("`vce'")
    if ("`vce'" == "") local vce "hc0"
    if !inlist("`vce'", "hc0", "hc1", "hc2", "hc3") {
        di as err "vce() must be hc0, hc1, hc2, or hc3"
        exit 198
    }
    if ("`clustername'" != "" & !inlist("`vce'", "hc0", "hc1")) {
        di as txt "note: vce(`vce') not available with cluster(); using vce(hc1) cluster-robust"
        local vce "hc1"
    }

    local stdflag = 0
    if ("`stdvars'" != "") local stdflag = 1

    local atlist "`at'"
    local atcount : word count `atlist'
    if (`atcount' != 2) {
        di as err "at() must contain exactly two numbers"
        exit 198
    }
    local at1 : word 1 of `atlist'
    local at2 : word 2 of `atlist'
    capture confirm number `at1'
    if (_rc) {
        di as err "at() must contain exactly two numbers"
        exit 198
    }
    capture confirm number `at2'
    if (_rc) {
        di as err "at() must contain exactly two numbers"
        exit 198
    }
    if (`at1' >= . | `at2' >= .) {
        di as err "at() must contain exactly two finite numbers"
        exit 198
    }

    quietly count if `touse'
    if (r(N) == 0) {
        di as err "no observations"
        exit 2000
    }
    quietly count if `touse' & !inlist(`dvar', 0, 1)
    if (r(N) > 0) {
        di as err "treatment indicator must contain only 0 and 1"
        exit 198
    }

    tempname b V bread meat residuals weights sd1 sd2 nh rank cond fallback
    mata: _rd2d_loc_fit_mata("`yvar'", "`x1var'", "`x2var'", "`dvar'", "`touse'", ///
        `at1', `at2', `hx', `hyval', `p', "`side'", "`kernel'", "`ktype'", "`vce'", ///
        `stdflag', "`clusterwork'", "`b'", "`V'", "`bread'", "`meat'", "`residuals'", "`weights'", ///
        "`sd1'", "`sd2'", "`nh'", "`rank'", "`cond'", "`fallback'")

    return scalar N_h = `nh'
    return scalar rank = `rank'
    return scalar condition = `cond'
    return scalar fallback = `fallback'
    return scalar p = `p'
    return scalar hx = `hx'
    return scalar hy = `hyval'
    return scalar sd_x1 = `sd1'
    return scalar sd_x2 = `sd2'
    return local side "`side'"
    return local kernel "`kernel'"
    return local ktype "`ktype'"
    return local vce "`vce'"
    return local cluster "`clustername'"
    return local stdvars "`stdflag'"
    return matrix b = `b'
    return matrix V = `V'
    return matrix bread = `bread'
    return matrix meat = `meat'
    return matrix residuals = `residuals'
    return matrix weights = `weights'
end

mata:
real colvector _rd2d_loc_kernel(real colvector u, real scalar kflag)
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

real matrix _rd2d_loc_basis(real colvector x1, real colvector x2, real scalar p)
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

real matrix _rd2d_loc_bread_inverse(real matrix bread, real rowvector diag)
{
    real matrix A

    A = (bread + bread') / 2
    if (diag[3] != 0) return(pinv(A))
    return(invsym(A))
}

void _rd2d_loc_fit_mata(
    string scalar yname,
    string scalar x1name,
    string scalar x2name,
    string scalar dname,
    string scalar tousename,
    real scalar at1,
    real scalar at2,
    real scalar h1,
    real scalar h2,
    real scalar p,
    string scalar side,
    string scalar kernel,
    string scalar ktype,
    string scalar vce,
    real scalar stdvars,
    string scalar cname,
    string scalar bname,
    string scalar vname,
    string scalar breadname,
    string scalar meatname,
    string scalar residualsname,
    string scalar weightsname,
    string scalar sd1name,
    string scalar sd2name,
    string scalar nhname,
    string scalar rankname,
    string scalar condname,
    string scalar fallbackname)
{
    real colvector y, x1, x2, d, C, xc1, xc2, fituse, w, resid, hii, adj
    real colvector clusters, cidx
    real matrix X, bread, ibread, beta, Xw, Xwr, V, meat, score, Xwr_b
    real rowvector diag
    real scalar n, k, hr, sd1, sd2, hascluster, g, j, factor, kflag
    real scalar blocksize, bstart, bend

    kflag = (kernel == "uniform") * 1 + (kernel == "triangular") * 2 + (kernel == "epanechnikov") * 3 + (kernel == "gaussian") * 4
    if (kflag == 0) kflag = 4

    st_view(y = ., ., yname, tousename)
    st_view(x1 = ., ., x1name, tousename)
    st_view(x2 = ., ., x2name, tousename)
    st_view(d = ., ., dname, tousename)
    hascluster = (cname != "")
    if (hascluster) st_view(C = ., ., cname, tousename)
    else C = J(rows(y), 1, .)

    sd1 = 1
    sd2 = 1
    if (stdvars != 0) {
        sd1 = sqrt(variance(x1))
        sd2 = sqrt(variance(x2))
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

    st_numscalar(sd1name, sd1)
    st_numscalar(sd2name, sd2)

    xc1 = x1 :- at1
    xc2 = x2 :- at2

    if (side == "control") {
        fituse = (d :== 0)
    }
    else {
        fituse = (d :== 1)
    }

    if (ktype == "prod") {
        w = _rd2d_loc_kernel(xc1 :/ h1, kflag) :* _rd2d_loc_kernel(xc2 :/ h2, kflag) :/ (h1 * h2)
    }
    else {
        hr = sqrt(h1^2 + h2^2)
        w = _rd2d_loc_kernel(sqrt(xc1:^2 + xc2:^2) :/ hr, kflag) :/ (hr^2)
    }

    fituse = fituse :& (w :> 0)
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
        errprintf("not enough observations for hc1 degrees-of-freedom correction\n")
        _error(2001)
    }

    X = _rd2d_loc_basis(xc1, xc2, p)
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
        if (vce == "hc2") {
            resid = resid :* sqrt(1 :/ (1 :- hii))
        }
        else {
            resid = resid :* (1 :/ (1 :- hii))
        }
    }

    if (!hascluster) {
        if (n > 200000) {
            // Blocked accumulation for large samples
            // Mathematically equivalent: A'A = sum_b (A_b'A_b) via associativity
            blocksize = 50000
            meat = J(k, k, 0)
            for (bstart = 1; bstart <= n; bstart = bstart + blocksize) {
                bend = min((bstart + blocksize - 1, n))
                Xwr_b = X[bstart..bend, .] :* (w[bstart..bend] :* resid[bstart..bend])
                meat = meat + quadcross(Xwr_b, Xwr_b)
            }
        }
        else {
            Xwr = X :* (w :* resid)
            meat = quadcross(Xwr, Xwr)
        }
    }
    else {
        clusters = uniqrows(sort(C, 1))
        g = rows(clusters)
        if (g <= 1) {
            errprintf("cluster() must contain at least two clusters inside the bandwidth on the %s side\n", side)
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
        meat = meat * factor
    }
    V = ibread * meat * ibread

    st_matrix(bname, beta')
    st_matrix(vname, V)
    st_matrix(breadname, bread)
    st_matrix(meatname, meat)
    st_matrix(residualsname, resid')
    st_matrix(weightsname, w')
    st_numscalar(nhname, n)
    st_numscalar(rankname, diag[1])
    st_numscalar(condname, diag[2])
    st_numscalar(fallbackname, diag[3])
}
end
