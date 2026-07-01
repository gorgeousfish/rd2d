version 16.0
clear all
set more off

// =====================================================================
// rd2d_cluster_example.do
// 展示 cluster-robust VCE 的使用 / Cluster-robust inference demo
// 生成含cluster分组的合成数据，对比有无cluster的标准误
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

// --- 合成数据生成 / Synthetic data with cluster structure ---
set seed 20260624
set obs 600
gen int region = ceil(_n / 25)
gen double x1 = runiform() * 2 - 1
gen double x2 = runiform() * 2 - 1
gen byte d = (x1 >= 0 & x2 >= 0)

// cluster-level random effect (region-level intercept shift)
gen double re = .
forvalues r = 1/24 {
    replace re = rnormal() * 0.4 in 1/600 if region == `r'
}

gen double y = 1 + 0.5 * x1 - 0.3 * x2 + 1.5 * d + 0.1 * x1 * x2 + re + rnormal() * 0.15

// --- (1) rd2d 无 cluster / without cluster ---
rd2d y x1 x2 d, at(0 0) h(1) p(1) q(2) kernel(triangular) ///
    ktype(prod) vce(hc1) bwcheck(20) masspoints(off)
matrix res_noclust = e(results)
local se_noclust = res_noclust[1, 6]
assert "`e(clustered)'" == "off"
assert "`e(cmd)'" == "rd2d"

// --- (2) rd2d 带 cluster / with cluster ---
rd2d y x1 x2 d, at(0 0) h(1) p(1) q(2) kernel(triangular) ///
    ktype(prod) vce(hc1) cluster(region) bwcheck(20) masspoints(off)
matrix res_clust = e(results)
local se_clust = res_clust[1, 6]

// 验证 cluster 返回值 / Verify cluster returns
assert "`e(clustered)'" == "on"
assert "`e(cluster)'" == "region"
assert "`e(cmd)'" == "rd2d"
assert e(N) == 600
assert e(neval) == 1
assert rowsof(res_clust) == 1
assert colsof(res_clust) == 18

// cluster SE 通常更大（组内相关导致） / Cluster SE is typically larger
display as txt "rd2d SE comparison: no-cluster=" as res %9.6f `se_noclust' ///
    as txt ", cluster=" as res %9.6f `se_clust'

// --- (3) rd2d 带 cluster 和 cbands / with cluster and confidence bands ---
set seed 20260624
rd2d y x1 x2 d, at(0 0) h(1) p(1) q(2) kernel(triangular) ///
    ktype(prod) vce(hc1) cluster(region) bwcheck(20) masspoints(off) ///
    cbands repp(500)
matrix res_cb = e(results)
assert "`e(clustered)'" == "on"
assert "`e(cbands)'" == "on"
assert e(cb_crit) < .
assert res_cb[1, 11] < .
assert res_cb[1, 12] > res_cb[1, 11]

// --- (4) rd2d_dist 带 cluster / signed-distance path with cluster ---
gen double dist = cond(d == 1, sqrt(x1^2 + x2^2), -sqrt(x1^2 + x2^2))

rd2d_dist y dist, h(1) p(1) q(2) kernel(triangular) ///
    vce(hc1) bwcheck(20) masspoints(off)
matrix dres_noclust = e(results)
local dse_noclust = dres_noclust[1, 6]

rd2d_dist y dist, h(1) p(1) q(2) kernel(triangular) ///
    vce(hc1) cluster(region) bwcheck(20) masspoints(off)
matrix dres_clust = e(results)
local dse_clust = dres_clust[1, 6]

assert "`e(clustered)'" == "on"
assert "`e(cluster)'" == "region"
assert "`e(cmd)'" == "rd2d_dist"

display as txt "rd2d_dist SE comparison: no-cluster=" as res %9.6f `dse_noclust' ///
    as txt ", cluster=" as res %9.6f `dse_clust'

// --- (5) vce(hc2) + cluster 自动降级为 hc1 / auto-downgrade ---
rd2d y x1 x2 d, at(0 0) h(1) p(1) q(2) kernel(triangular) ///
    ktype(prod) vce(hc2) cluster(region) bwcheck(20) masspoints(off)
assert "`e(vce)'" == "hc1"
assert "`e(clustered)'" == "on"

// --- 格式化报告 / Formatted report ---
display _n as txt "=== Cluster-robust VCE example report ==="
display as txt "rd2d at(0,0): Est.q=" as res %9.6f res_clust[1, 5] ///
    as txt ", Se.q(noclust)=" as res %9.6f `se_noclust' ///
    as txt ", Se.q(cluster)=" as res %9.6f `se_clust' ///
    as txt ", clustered=" as res "`e(clustered)'" ///
    as txt ", cluster=" as res "`e(cluster)'"
display as txt "rd2d_dist: Est.q=" as res %9.6f dres_clust[1, 5] ///
    as txt ", Se.q(noclust)=" as res %9.6f `dse_noclust' ///
    as txt ", Se.q(cluster)=" as res %9.6f `dse_clust'
display as txt "rd2d cbands: cb_crit=" as res %9.6f e(cb_crit) ///
    as txt ", CB=[" as res %9.6f res_cb[1, 11] as txt ", " ///
    as res %9.6f res_cb[1, 12] as txt "]"

di as txt "EXIT:0"
