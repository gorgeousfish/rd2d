version 16.0
clear all
set more off

// =====================================================================
// rd2d_cbands_example.do
// 展示置信带(confidence bands)的使用 / Uniform confidence bands demo
// 对比 pointwise CI 与 uniform CB 的宽度差异
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

// --- (1) rd2d 多评估点，仅 pointwise CI / Multiple at() points, pointwise CI only ---
rd2d y x1 x2 d, at(0 0 0.1 0.3 -0.2 0.1) h(1.5) p(1) q(2) ///
    kernel(triangular) ktype(prod) vce(hc1) bwcheck(20) masspoints(off)
matrix res_ci = e(results)
assert "`e(cbands)'" == "off"
assert e(neval) == 3
assert rowsof(res_ci) == 3

// 提取 pointwise CI 宽度 / Extract pointwise CI widths
local ci_w1 = res_ci[1, 10] - res_ci[1, 9]
local ci_w2 = res_ci[2, 10] - res_ci[2, 9]
local ci_w3 = res_ci[3, 10] - res_ci[3, 9]

// --- (2) rd2d 带 cbands / with uniform confidence bands ---
set seed 20260624
rd2d y x1 x2 d, at(0 0 0.1 0.3 -0.2 0.1) h(1.5) p(1) q(2) ///
    kernel(triangular) ktype(prod) vce(hc1) bwcheck(20) masspoints(off) ///
    cbands repp(500)
matrix res_cb = e(results)
matrix cb_cov = e(cov_q)
matrix cb_corr = e(corr_q)

// 验证 cbands 返回值 / Verify cbands returns
assert "`e(cbands)'" == "on"
assert e(repp) == 500
assert e(cb_crit) < .
assert e(neval) == 3
assert rowsof(res_cb) == 3
assert colsof(res_cb) == 18
assert rowsof(cb_cov) == 3
assert colsof(cb_cov) == 3
assert rowsof(cb_corr) == 3
assert colsof(cb_corr) == 3

// 提取 uniform CB 宽度 / Extract uniform CB widths
local cb_w1 = res_cb[1, 12] - res_cb[1, 11]
local cb_w2 = res_cb[2, 12] - res_cb[2, 11]
local cb_w3 = res_cb[3, 12] - res_cb[3, 11]

// --- (3) 验证 CB 代数：CB = Est.q +/- cb_crit * Se.q ---
local cb_crit_val = e(cb_crit)
forvalues j = 1/3 {
    assert abs(res_cb[`j', 11] - (res_cb[`j', 5] - `cb_crit_val' * res_cb[`j', 6])) < 1e-8
    assert abs(res_cb[`j', 12] - (res_cb[`j', 5] + `cb_crit_val' * res_cb[`j', 6])) < 1e-8
}

// --- (4) 验证 cov_q 对角线 = Se.q^2 / Diagonal of cov_q equals Se.q^2 ---
forvalues j = 1/3 {
    assert abs(cb_cov[`j', `j'] - res_cb[`j', 6]^2) < 1e-10
}

// --- (5) 对比 pointwise CI 与 uniform CB 宽度 / Compare widths ---
// uniform CB 始终宽于 pointwise CI（cb_crit > z_{0.975}）
local z_975 = invnormal(0.975)
display as txt "cb_crit=" as res %9.6f `cb_crit_val' ///
    as txt " vs z_0.975=" as res %9.6f `z_975'
assert `cb_crit_val' > `z_975'

// 宽度比 / Width ratio
local ratio1 = `cb_w1' / `ci_w1'
display as txt "CB/CI width ratio at at1: " as res %6.3f `ratio1'

// --- (6) rd2d_dist 带 cbands / distance path with cbands ---
gen double dist = cond(d == 1, sqrt(x1^2 + x2^2), -sqrt(x1^2 + x2^2))
gen double dist2 = sign(dist) * (abs(dist) + 0.1)

set seed 20260624
rd2d_dist y dist dist2, h(1.5) p(1) q(2) kernel(triangular) ///
    vce(hc1) bwcheck(20) masspoints(off) cbands repp(500)
matrix dres_cb = e(results)
assert "`e(cbands)'" == "on"
assert e(neval) == 2
assert e(cb_crit) < .

// distance path CB 代数同样成立 / CB algebra holds for distance path too
local dcb_crit = e(cb_crit)
forvalues j = 1/2 {
    assert abs(dres_cb[`j', 11] - (dres_cb[`j', 5] - `dcb_crit' * dres_cb[`j', 6])) < 1e-8
    assert abs(dres_cb[`j', 12] - (dres_cb[`j', 5] + `dcb_crit' * dres_cb[`j', 6])) < 1e-8
}

// --- 格式化报告 / Formatted report ---
display _n as txt "=== Confidence bands example report ==="
display as txt "rd2d 3-point cbands: cb_crit=" as res %9.6f `cb_crit_val' ///
    as txt ", z_0.975=" as res %9.6f `z_975'
forvalues j = 1/3 {
    display as txt "  at`j': Est.q=" as res %9.6f res_cb[`j', 5] ///
        as txt " CI=[" as res %9.6f res_ci[`j', 9] as txt "," ///
        as res %9.6f res_ci[`j', 10] as txt "]" ///
        as txt " CB=[" as res %9.6f res_cb[`j', 11] as txt "," ///
        as res %9.6f res_cb[`j', 12] as txt "]"
}
display as txt "rd2d_dist 2-col cbands: cb_crit=" as res %9.6f `dcb_crit'
forvalues j = 1/2 {
    display as txt "  dist`j': Est.q=" as res %9.6f dres_cb[`j', 5] ///
        as txt " CB=[" as res %9.6f dres_cb[`j', 11] as txt "," ///
        as res %9.6f dres_cb[`j', 12] as txt "]"
}

di as txt "EXIT:0"
