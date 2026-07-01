version 16.0
clear all
set more off

// =====================================================================
// rd2d_imse_example.do
// 展示 IMSE 最优带宽选择 / Integrated MSE bandwidth selection demo
// 对比 mserd (per-point) vs imserd (integrated) vs imsetwo 带宽结构
// =====================================================================

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

// --- 合成数据生成 / Synthetic data ---
set seed 20260624
set obs 500
gen double row = mod(_n - 1, 25) - 12
gen double col = floor((_n - 1) / 25) - 9
gen double x1 = row / 6
gen double x2 = col / 5
gen byte d = (x1 + 0.3 * x2 >= 0)
gen double y = 1 + 0.4 * x1 - 0.2 * x2 + 1.3 * d + 0.08 * x1 * x2 + rnormal() * 0.1

// --- (1) rdbw2d: mserd 每点选不同带宽 / per-point MSE bandwidth ---
rdbw2d y x1 x2 d, at(0 0) p(1) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(rot) vce(hc1) bwcheck(20) masspoints(off)
matrix bw_mserd_1 = r(bws)

rdbw2d y x1 x2 d, at(0.1 0.3) p(1) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(rot) vce(hc1) bwcheck(20) masspoints(off)
matrix bw_mserd_2 = r(bws)

rdbw2d y x1 x2 d, at(-0.2 0.1) p(1) kernel(triangular) ktype(prod) ///
    bwselect(mserd) method(rot) vce(hc1) bwcheck(20) masspoints(off)
matrix bw_mserd_3 = r(bws)

// mserd: 每个点有自己的 4 个带宽 / Each point gets its own 4 bandwidths
display as txt "mserd per-point bandwidths:"
display as txt "  at(0,0):    h01=" as res %6.3f bw_mserd_1[1, 3] ///
    as txt " h02=" as res %6.3f bw_mserd_1[1, 4] ///
    as txt " h11=" as res %6.3f bw_mserd_1[1, 5] ///
    as txt " h12=" as res %6.3f bw_mserd_1[1, 6]
display as txt "  at(0.1,0.3): h01=" as res %6.3f bw_mserd_2[1, 3] ///
    as txt " h02=" as res %6.3f bw_mserd_2[1, 4] ///
    as txt " h11=" as res %6.3f bw_mserd_2[1, 5] ///
    as txt " h12=" as res %6.3f bw_mserd_2[1, 6]
display as txt "  at(-0.2,0.1): h01=" as res %6.3f bw_mserd_3[1, 3] ///
    as txt " h02=" as res %6.3f bw_mserd_3[1, 4] ///
    as txt " h11=" as res %6.3f bw_mserd_3[1, 5] ///
    as txt " h12=" as res %6.3f bw_mserd_3[1, 6]

// --- (2) rdbw2d: imserd 积分MSE统一带宽 / integrated MSE single bandwidth ---
rdbw2d y x1 x2 d, at(0 0) p(1) kernel(triangular) ktype(prod) ///
    bwselect(imserd) method(rot) vce(hc1) bwcheck(20) masspoints(off)
matrix bw_imserd_1 = r(bws)

// 验证 imserd 结构: 6 列矩阵 / Verify structure: 6-column matrix
assert colsof(bw_imserd_1) == 6
assert bw_imserd_1[1, 3] > 0

// --- (3) rd2d: imserd 多点自动带宽 / automatic IMSE for multiple points ---
rd2d y x1 x2 d, at(0 0 0.1 0.3 -0.2 0.1) p(1) q(2) ///
    kernel(triangular) ktype(prod) vce(hc1) bwcheck(20) masspoints(off) ///
    bwselect(imserd) method(rot) rbc(on)
matrix res_imserd = e(results)
matrix bws_imserd = e(bws)

assert "`e(bwselect)'" == "imserd"
assert "`e(bwsource)'" == "automatic"
assert e(neval) == 3
assert rowsof(bws_imserd) == 3
assert colsof(bws_imserd) == 8

// imserd 关键性质: 所有点使用相同的带宽 / Key property: same h across all points
local imse_h01 = bws_imserd[1, 3]
local imse_h02 = bws_imserd[1, 4]
local imse_h11 = bws_imserd[1, 5]
local imse_h12 = bws_imserd[1, 6]
forvalues j = 1/3 {
    assert abs(bws_imserd[`j', 3] - `imse_h01') < 1e-10
    assert abs(bws_imserd[`j', 4] - `imse_h02') < 1e-10
    assert abs(bws_imserd[`j', 5] - `imse_h11') < 1e-10
    assert abs(bws_imserd[`j', 6] - `imse_h12') < 1e-10
}

display as txt "imserd uniform bandwidths (same for all 3 points):"
display as txt "  h01=" as res %6.3f `imse_h01' ///
    as txt " h02=" as res %6.3f `imse_h02' ///
    as txt " h11=" as res %6.3f `imse_h11' ///
    as txt " h12=" as res %6.3f `imse_h12'

// --- (4) rd2d: mserd 多点各自带宽 / per-point MSE for comparison ---
rd2d y x1 x2 d, at(0 0 0.1 0.3 -0.2 0.1) p(1) q(2) ///
    kernel(triangular) ktype(prod) vce(hc1) bwcheck(20) masspoints(off) ///
    bwselect(mserd) method(rot) rbc(on)
matrix res_mserd = e(results)
matrix bws_mserd = e(bws)

assert "`e(bwselect)'" == "mserd"
assert e(neval) == 3

// mserd: 各点带宽不同 / mserd: each point has different bandwidths
local mserd_h01_1 = bws_mserd[1, 3]
local mserd_h01_2 = bws_mserd[2, 3]
local mserd_h01_3 = bws_mserd[3, 3]
display as txt "mserd per-point h01: at1=" as res %6.3f `mserd_h01_1' ///
    as txt " at2=" as res %6.3f `mserd_h01_2' ///
    as txt " at3=" as res %6.3f `mserd_h01_3'

// --- (5) imsetwo: 侧别不同的统一带宽 / side-specific integrated bandwidth ---
rd2d y x1 x2 d, at(0 0 0.1 0.3 -0.2 0.1) p(1) q(2) ///
    kernel(triangular) ktype(prod) vce(hc1) bwcheck(20) masspoints(off) ///
    bwselect(imsetwo) method(rot) rbc(on)
matrix res_imsetwo = e(results)
matrix bws_imsetwo = e(bws)

assert "`e(bwselect)'" == "imsetwo"
assert e(neval) == 3

// imsetwo: 所有点共享同一对 (h0, h1)，但 h0 ≠ h1 / All points share (h0, h1), h0 != h1
local imsetwo_h01 = bws_imsetwo[1, 3]
local imsetwo_h11 = bws_imsetwo[1, 5]
forvalues j = 1/3 {
    assert abs(bws_imsetwo[`j', 3] - `imsetwo_h01') < 1e-10
    assert abs(bws_imsetwo[`j', 5] - `imsetwo_h11') < 1e-10
}
display as txt "imsetwo side-specific bandwidths (same across points):"
display as txt "  h0=" as res %6.3f `imsetwo_h01' ///
    as txt " h1=" as res %6.3f `imsetwo_h11'

// --- (6) r(bws) 结构差异总结 / Structure difference summary ---
display _n as txt "=== IMSE bandwidth example report ==="
display as txt "r(bws) structure for 3 evaluation points:"
display as txt "  mserd:   3 rows x 6 cols, each row has point-specific bandwidths"
display as txt "  imserd:  3 rows x 6 cols, all rows share identical bandwidths"
display as txt "  imsetwo: 3 rows x 6 cols, all rows share h0/h1 but h0 != h1"
display as txt ""
display as txt "imserd: Est.q at1=" as res %9.6f res_imserd[1, 5] ///
    as txt " at2=" as res %9.6f res_imserd[2, 5] ///
    as txt " at3=" as res %9.6f res_imserd[3, 5]
display as txt "mserd:  Est.q at1=" as res %9.6f res_mserd[1, 5] ///
    as txt " at2=" as res %9.6f res_mserd[2, 5] ///
    as txt " at3=" as res %9.6f res_mserd[3, 5]

di as txt "EXIT:0"
