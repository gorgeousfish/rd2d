version 16.0
clear all
set more off

capture which rdbw2d
if _rc {
    foreach pkgpath in "rd2d-stata/ado" "rd2d-stata/pkg/ado" "ado" "pkg/ado" "../ado" "../pkg/ado" {
        capture confirm file "`pkgpath'/rdbw2d.ado"
        if !_rc {
            adopath ++ "`pkgpath'"
        }
    }
    foreach helppath in "rd2d-stata/help" "rd2d-stata/pkg/help" "help" "pkg/help" "../help" "../pkg/help" {
        capture confirm file "`helppath'/rdbw2d.sthlp"
        if !_rc {
            adopath ++ "`helppath'"
        }
    }
}

set obs 400
generate double row = mod(_n - 1, 20) - 9.5
generate double col = floor((_n - 1) / 20) - 9.5
generate double x1 = row / 4
generate double x2 = col / 4
generate byte treat = x1 >= 0
generate double y = 1 + .5 * x1 - .25 * x2 + 1.2 * treat + .1 * x1 * x2 + .05 * mod(_n, 7)

rdbw2d y x1 x2 treat, at(0 0) method(rot) kernel(triangular) ktype(prod) bwcheck(20)
matrix location_bws = r(bws)
matrix location_mse = r(mseconsts)
assert rowsof(location_bws) == 1
assert colsof(location_bws) == 6
assert rowsof(location_mse) == 1
assert colsof(location_mse) == 10
assert r(N) == 400
assert r(N0) == 200
assert r(N1) == 200
assert r(p) == 1
assert "`r(kernel)'" == "triangular"
assert "`r(ktype)'" == "prod"
assert "`r(method)'" == "rot"
assert location_bws[1, 3] > 0
assert location_bws[1, 4] > 0
assert location_bws[1, 5] > 0
assert location_bws[1, 6] > 0

rd2d y x1 x2 treat, at(0 0) h(1.25) p(1) q(2) kernel(triangular) ///
    ktype(prod) bwcheck(20) cbands repp(200)
matrix location_results = e(results)
matrix location_fit_bws = e(bws)
matrix location_diag = e(diagnostics)
matrix location_cov = e(cov_q)
matrix location_corr = e(corr_q)
local location_target "at(0,0)"
local location_fallback "`e(fallback)'"
local location_h01 = location_fit_bws[1, 3]
local location_h02 = location_fit_bws[1, 4]
local location_h11 = location_fit_bws[1, 5]
local location_h12 = location_fit_bws[1, 6]
local location_nh0 = location_fit_bws[1, 7]
local location_nh1 = location_fit_bws[1, 8]
assert rowsof(location_results) == 1
assert colsof(location_results) == 18
assert rowsof(location_fit_bws) == 1
assert colsof(location_fit_bws) == 8
assert rowsof(location_diag) == 1
assert colsof(location_diag) == 14
assert rowsof(location_cov) == 1
assert rowsof(location_corr) == 1
assert e(N) == 400
assert e(N0) == 200
assert e(N1) == 200
assert e(neval) == 1
assert e(p) == 1
assert e(q) == 2
assert e(repp) == 200
assert "`e(cmd)'" == "rd2d"
assert "`e(cbands)'" == "on"
assert "`e(bwsource)'" == "user"
assert "`e(kernel)'" == "triangular"
assert "`e(ktype)'" == "prod"
assert location_results[1, 5] < .
assert location_results[1, 6] > 0
assert location_results[1, 11] < .
assert location_results[1, 12] < .
assert location_fit_bws[1, 3] == 1.25
assert location_fit_bws[1, 4] == 1.25
assert location_fit_bws[1, 5] == 1.25
assert location_fit_bws[1, 6] == 1.25
assert location_diag[1, 5] == 0
assert location_diag[1, 8] == 0
assert location_diag[1, 11] == 0
assert location_diag[1, 14] == 0

twoway ///
    (scatter y x1 if treat == 0, msymbol(oh) msize(vsmall) mcolor(navy)) ///
    (scatter y x1 if treat == 1, msymbol(oh) msize(vsmall) mcolor(maroon)) ///
    (lfit y x1 if treat == 0, lcolor(navy) lwidth(medthick)) ///
    (lfit y x1 if treat == 1, lcolor(maroon) lwidth(medthick)), ///
    scheme(s1color) ///
    xline(0, lcolor(gs8) lpattern(dash) lwidth(medium)) ///
    title("Location RD example", size(medsmall)) ///
    subtitle("Outcome by running coordinate x1", size(small)) ///
    xtitle("Running coordinate x1", size(small)) ///
    ytitle("Outcome", size(small)) ///
    xlabel(, labsize(small) grid glcolor(gs14)) ///
    ylabel(, labsize(small) angle(horizontal) grid glcolor(gs14)) ///
    legend(order(1 "Control observations" 2 "Treated observations" ///
        3 "Control fit" 4 "Treated fit") rows(2) position(6) size(vsmall) ///
        region(lstyle(none))) ///
    note("Synthetic N=400 example (200 control / 200 treated); side-specific linear fits; dashed line marks the boundary.", ///
        size(vsmall)) ///
    xsize(8) ysize(4.8) ///
    graphregion(color(white)) plotregion(color(white) lcolor(gs12)) ///
    name(rd2d_location_overview, replace)
capture graph display rd2d_location_overview
assert _rc == 0

display as txt "rd2d location example reporting row: target=" as res "`location_target'" ///
    as txt ", Est.q=" as res %9.6f location_results[1, 5] ///
    as txt ", Se.q=" as res %9.6f location_results[1, 6] ///
    as txt ", h=(" as res %5.2f `location_h01' as txt "," ///
    as res %5.2f `location_h02' as txt "," ///
    as res %5.2f `location_h11' as txt "," ///
    as res %5.2f `location_h12' as txt ")" ///
    as txt ", Nh=(" as res %4.0f `location_nh0' as txt "," ///
    as res %4.0f `location_nh1' as txt ")" ///
    as txt ", fallback=" as res "`location_fallback'" ///
    as txt ", graph=rd2d_location_overview."
