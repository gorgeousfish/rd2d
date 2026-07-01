version 16.0
clear all
set more off

capture which rdbw2d_dist
if _rc {
    foreach pkgpath in "rd2d-stata/ado" "rd2d-stata/pkg/ado" "ado" "pkg/ado" "../ado" "../pkg/ado" {
        capture confirm file "`pkgpath'/rdbw2d_dist.ado"
        if !_rc {
            adopath ++ "`pkgpath'"
        }
    }
    foreach helppath in "rd2d-stata/help" "rd2d-stata/pkg/help" "help" "pkg/help" "../help" "../pkg/help" {
        capture confirm file "`helppath'/rdbw2d_dist.sthlp"
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
generate double dist0 = x1
generate byte treat = dist0 >= 0
generate double y = 1 + .5 * x1 - .25 * x2 + 1.2 * treat + .1 * x1 * x2 + .05 * mod(_n, 7)

rdbw2d_dist y dist0, kernel(triangular) bwcheck(20)
matrix distance_bws = r(bws)
matrix distance_mse = r(mseconsts)
matrix distance_mass = r(masspoints)
assert rowsof(distance_bws) == 1
assert colsof(distance_bws) == 6
assert rowsof(distance_mse) == 1
assert colsof(distance_mse) == 14
assert rowsof(distance_mass) == 1
assert colsof(distance_mass) == 4
assert r(N) == 400
assert r(N0) == 200
assert r(N1) == 200
assert r(neval) == 1
assert r(p) == 1
assert "`r(kernel)'" == "triangular"
assert "`r(kink)'" == "off"
assert distance_bws[1, 3] > 0
assert distance_bws[1, 4] > 0
assert distance_bws[1, 5] > 0
assert distance_bws[1, 6] > 0

rd2d_dist y dist0, h(1.25) p(1) q(2) kernel(triangular) bwcheck(20) ///
    cbands repp(200)
matrix distance_results = e(results)
matrix distance_fit_bws = e(bws)
matrix distance_diag = e(diagnostics)
matrix distance_cov = e(cov_q)
matrix distance_corr = e(corr_q)
local distance_target "dist0"
local distance_kink "`e(kink)'"
local distance_fallback "`e(fallback)'"
local distance_h0 = distance_fit_bws[1, 3]
local distance_h1 = distance_fit_bws[1, 4]
local distance_nh0 = distance_fit_bws[1, 5]
local distance_nh1 = distance_fit_bws[1, 6]
assert rowsof(distance_results) == 1
assert colsof(distance_results) == 18
assert rowsof(distance_fit_bws) == 1
assert colsof(distance_fit_bws) == 6
assert rowsof(distance_diag) == 1
assert colsof(distance_diag) == 14
assert rowsof(distance_cov) == 1
assert rowsof(distance_corr) == 1
assert e(N) == 400
assert e(N0) == 200
assert e(N1) == 200
assert e(neval) == 1
assert e(p) == 1
assert e(q) == 2
assert e(repp) == 200
assert "`e(cmd)'" == "rd2d_dist"
assert "`e(cbands)'" == "on"
assert "`e(bwsource)'" == "user"
assert "`e(kernel)'" == "triangular"
assert "`e(kink)'" == "off"
assert distance_results[1, 5] < .
assert distance_results[1, 6] > 0
assert distance_results[1, 11] < .
assert distance_results[1, 12] < .
assert distance_fit_bws[1, 3] == 1.25
assert distance_fit_bws[1, 4] == 1.25
assert distance_diag[1, 5] == 0
assert distance_diag[1, 8] == 0
assert distance_diag[1, 11] == 0
assert distance_diag[1, 14] == 0

rdbw2d_dist y dist0, p(1) kink(on) kernel(triangular) bwcheck(20)
matrix distance_kink_selector_bws = r(bws)
matrix distance_kink_selector_mass = r(masspoints)
assert rowsof(distance_kink_selector_bws) == 1
assert colsof(distance_kink_selector_bws) == 6
assert rowsof(distance_kink_selector_mass) == 1
assert colsof(distance_kink_selector_mass) == 4
assert "`r(kink)'" == "on"
assert distance_kink_selector_bws[1, 3] > 0
assert distance_kink_selector_bws[1, 4] > 0

rd2d_dist y dist0, h(1.25) p(1) kink(on) kernel(triangular) bwcheck(20)
matrix distance_kink_results = e(results)
matrix distance_kink_fit_bws = e(bws)
matrix distance_kink_diag = e(diagnostics)
local distance_kink_target "dist0"
local distance_kink_status "`e(kink)'"
local distance_kink_fallback "`e(fallback)'"
local distance_kink_q = e(q)
local distance_kink_h0 = distance_kink_fit_bws[1, 3]
local distance_kink_h1 = distance_kink_fit_bws[1, 4]
local distance_kink_nh0 = distance_kink_fit_bws[1, 5]
local distance_kink_nh1 = distance_kink_fit_bws[1, 6]
assert rowsof(distance_kink_results) == 1
assert colsof(distance_kink_results) == 18
assert rowsof(distance_kink_fit_bws) == 1
assert colsof(distance_kink_fit_bws) == 6
assert rowsof(distance_kink_diag) == 1
assert colsof(distance_kink_diag) == 14
assert e(N) == 400
assert e(p) == 1
assert e(q) == e(p)
assert "`e(cmd)'" == "rd2d_dist"
assert "`e(kink)'" == "on"
assert "`e(bwsource)'" == "user"
assert "`e(kernel)'" == "triangular"
assert distance_kink_results[1, 5] < .
assert distance_kink_results[1, 6] > 0
assert distance_kink_fit_bws[1, 3] == 1.25
assert distance_kink_fit_bws[1, 4] == 1.25
assert distance_kink_diag[1, 5] == 0
assert distance_kink_diag[1, 8] == 0
assert distance_kink_diag[1, 11] == 0
assert distance_kink_diag[1, 14] == 0

twoway ///
    (scatter y dist0 if dist0 < 0, msymbol(oh) msize(vsmall) mcolor(navy)) ///
    (scatter y dist0 if dist0 >= 0, msymbol(oh) msize(vsmall) mcolor(maroon)) ///
    (lfit y dist0 if dist0 < 0, lcolor(navy) lwidth(medthick)) ///
    (lfit y dist0 if dist0 >= 0, lcolor(maroon) lwidth(medthick)), ///
    scheme(s1color) ///
    xline(0, lcolor(gs8) lpattern(dash) lwidth(medium)) ///
    title("Distance RD example", size(medsmall)) ///
    subtitle("Outcome by signed distance", size(small)) ///
    xtitle("Signed distance", size(small)) ///
    ytitle("Outcome", size(small)) ///
    xlabel(, labsize(small) grid glcolor(gs14)) ///
    ylabel(, labsize(small) angle(horizontal) grid glcolor(gs14)) ///
    legend(order(1 "Control observations" 2 "Treated observations" ///
        3 "Control fit" 4 "Treated fit") rows(2) position(6) size(vsmall) ///
        region(lstyle(none))) ///
    note("Synthetic N=400 example (200 control / 200 treated); side-specific linear fits; dashed line marks zero distance.", ///
        size(vsmall)) ///
    xsize(8) ysize(4.8) ///
    graphregion(color(white)) plotregion(color(white) lcolor(gs12)) ///
    name(rd2d_distance_overview, replace)
capture graph display rd2d_distance_overview
assert _rc == 0

display as txt "rd2d distance example reporting row: target=" as res "`distance_target'" ///
    as txt ", kink=" as res "`distance_kink'" ///
    as txt ", Est.q=" as res %9.6f distance_results[1, 5] ///
    as txt ", Se.q=" as res %9.6f distance_results[1, 6] ///
    as txt ", h=(" as res %5.2f `distance_h0' as txt "," ///
    as res %5.2f `distance_h1' as txt ")" ///
    as txt ", Nh=(" as res %4.0f `distance_nh0' as txt "," ///
    as res %4.0f `distance_nh1' as txt ")" ///
    as txt ", fallback=" as res "`distance_fallback'" ///
    as txt ", graph=rd2d_distance_overview."

display as txt "rd2d distance kink reporting row: target=" as res "`distance_kink_target'" ///
    as txt ", kink=" as res "`distance_kink_status'" ///
    as txt ", q=" as res %2.0f `distance_kink_q' ///
    as txt ", Est.q=" as res %9.6f distance_kink_results[1, 5] ///
    as txt ", Se.q=" as res %9.6f distance_kink_results[1, 6] ///
    as txt ", h=(" as res %5.2f `distance_kink_h0' as txt "," ///
    as res %5.2f `distance_kink_h1' as txt ")" ///
    as txt ", Nh=(" as res %4.0f `distance_kink_nh0' as txt "," ///
    as res %4.0f `distance_kink_nh1' as txt ")" ///
    as txt ", fallback=" as res "`distance_kink_fallback'" ///
    as txt ", graph=rd2d_distance_overview."
