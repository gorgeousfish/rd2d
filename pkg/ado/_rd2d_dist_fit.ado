*! version 1.1.0 24jun2026
program define _rd2d_dist_fit, rclass
    version 16.0
    syntax varlist(min=2 max=2 numeric) [if] [in], H(real) P(integer) Side(string) ///
        [KERnel(string) VCE(string)]

    gettoken yvar dvar : varlist
    marksample touse
    markout `touse' `varlist'

    if (`h' >= . | `h' <= 0) {
        di as err "h() must be finite and positive"
        exit 198
    }
    if (`p' < 0) {
        di as err "p() must be a nonnegative integer"
        exit 198
    }

    local side = lower("`side'")
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

    local vce = lower("`vce'")
    if ("`vce'" == "") local vce "hc0"
    if !inlist("`vce'", "hc0", "hc1", "hc2", "hc3") {
        di as err "vce() must be hc0, hc1, hc2, or hc3"
        exit 198
    }

    tempvar fituse absd
    quietly gen double `absd' = abs(`dvar') if `touse'
    if ("`side'" == "control") {
        quietly gen byte `fituse' = `touse' & `dvar' < 0
    }
    else {
        quietly gen byte `fituse' = `touse' & `dvar' >= 0
    }
    if ("`kernel'" == "uniform") {
        quietly replace `fituse' = `fituse' & `absd' <= `h'
    }
    else if inlist("`kernel'", "triangular", "epanechnikov") {
        quietly replace `fituse' = `fituse' & `absd' < `h'
    }

    quietly count if `fituse'
    local nh = r(N)
    if (`nh' < `p' + 1) {
        di as err "insufficient observations on `side' side (Nh=`nh', min=`=`p'+1'); consider increasing h() or reducing p()"
        exit 2001
    }
    if ("`vce'" == "hc1" & `nh' <= `p' + 1) {
        di as err "insufficient observations on `side' side for hc1 degrees-of-freedom correction (Nh=`nh', min=`=`p'+2'); consider increasing h() or reducing p()"
        exit 2001
    }

    tempname b V bread meat residuals weights rank cond fallback
    mata: _rd2d_dist_fit_mata("`yvar'", "`dvar'", "`fituse'", `h', `p', "`kernel'", "`vce'", ///
        "`b'", "`V'", "`bread'", "`meat'", "`residuals'", "`weights'", ///
        "`rank'", "`cond'", "`fallback'")

    return scalar N_h = `nh'
    return scalar h = `h'
    return scalar p = `p'
    return scalar mu = `b'[1, 1]
    return scalar se = sqrt(`V'[1, 1])
    return scalar rank = `rank'
    return scalar condition = `cond'
    return scalar fallback = `fallback'
    return local side "`side'"
    return local kernel "`kernel'"
    return local vce "`vce'"
    return matrix b = `b'
    return matrix V = `V'
    return matrix bread = `bread'
    return matrix meat = `meat'
    return matrix residuals = `residuals'
    return matrix weights = `weights'
end

mata:
real colvector _rd2d_dist_kernel(real colvector u, string scalar kernel)
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

void _rd2d_dist_fit_mata(
    string scalar yname,
    string scalar dname,
    string scalar tousename,
    real scalar h,
    real scalar p,
    string scalar kernel,
    string scalar vce,
    string scalar bname,
    string scalar vname,
    string scalar breadname,
    string scalar meatname,
    string scalar residualsname,
    string scalar weightsname,
    string scalar rankname,
    string scalar condname,
    string scalar fallbackname)
{
    real colvector y, d, x, u, w, residuals, hii
    real matrix X, bread, ibread, beta, Xwr, meat, V, sqrtwX
    real rowvector diag
    real scalar n, j, k

    st_view(y = ., ., yname, tousename)
    st_view(d = ., ., dname, tousename)

    n = rows(y)
    k = p + 1
    x = abs(d)
    u = x :/ h
    w = _rd2d_dist_kernel(u, kernel) :/ (h^2)

    X = J(n, k, 1)
    for (j = 2; j <= k; j++) {
        X[, j] = x:^(j - 1)
    }

    bread = quadcross(X, w, X)
    diag = _rd2d_dist_bread_diag(bread)
    ibread = _rd2d_dist_bread_inverse(bread, diag)
    beta = ibread * quadcross(X, w, y)
    residuals = y - X * beta

    if (vce == "hc2" | vce == "hc3") {
        sqrtwX = sqrt(w) :* X
        hii = rowsum((sqrtwX * ibread) :* sqrtwX)
        if (min(1 :- hii) <= 0) {
            errprintf("leverage is too high for hc2/hc3 variance calculation\n")
            _error(2001)
        }
        if (vce == "hc2") {
            residuals = residuals :* sqrt(1 :/ (1 :- hii))
        }
        else {
            residuals = residuals :* (1 :/ (1 :- hii))
        }
    }

    Xwr = X :* (w :* residuals)
    meat = quadcross(Xwr, Xwr)
    if (vce == "hc1") {
        meat = meat * n / (n - k)
    }
    V = ibread * meat * ibread

    st_matrix(bname, beta')
    st_matrix(vname, V)
    st_matrix(breadname, bread)
    st_matrix(meatname, meat)
    st_matrix(residualsname, residuals')
    st_matrix(weightsname, w')
    st_numscalar(rankname, diag[1])
    st_numscalar(condname, diag[2])
    st_numscalar(fallbackname, diag[3])
}
end
