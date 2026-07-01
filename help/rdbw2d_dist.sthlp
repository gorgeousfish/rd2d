{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rdbw2d_dist}

{p 4 4 2}
{cmd:rdbw2d_dist} selects bandwidths for signed-distance boundary regression
discontinuity designs.

{title:Syntax}

{p 8 8 2}
{cmd:rdbw2d_dist} {it:yvar} {it:distvarlist} {ifin}
[{cmd:,}
{cmd:p(}{it:#}{cmd:)}
{cmd:kink(}{it:on|off}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:vce(}{it:hc0|hc1|hc2|hc3}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:check|adjust|off}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:cqt(}{it:#}{cmd:)}
{cmd:cluster(}{it:varname}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rdbw2d_dist} selects bandwidths for one or more signed-distance running
variables. Negative distances denote control observations and nonnegative
distances denote treated observations.

{p 4 4 2}
The target-defining inputs are the signed-distance variables. Each row in
{cmd:r(bws)} is indexed by the distance-variable name and the sign convention,
with {cmd:kink()} recording whether the nonsmooth-boundary path was used. The
selector has no boundary-coordinate input, so any coordinate used later in a
table is an application label supplied by the do-file.

{p 4 4 2}
The implementation follows the distance-based rate rules in Cattaneo,
Titiunik, and Yu: smooth-boundary selectors use the usual local-polynomial
rate, while {cmd:kink(on)} shrinks selected bandwidths to the nonsmooth-boundary
rate and should not be interpreted as robust bias correction.

{p 4 4 2}
After selection, {cmd:rdbw2d_dist} prints a compact bandwidth table with one
row per distance column, including side-specific bandwidths and effective
sample sizes. The complete numeric contract remains in {cmd:r(bws)},
{cmd:r(mseconsts)}, and {cmd:r(masspoints)}. Long distance-column names are
shortened in the printed table only, using a row suffix to keep labels
distinct; stored matrix row names keep the original variable names. Printed
tables adapt to narrower Stata {cmd:linesize} settings without changing stored
matrix names or columns. At compact {cmd:linesize(50)} widths, the printed
distance labels stay readable in the same compact table layout.

{p 4 4 2}
Kernel aliases follow the R interface where practical: {cmd:uni}, {cmd:unif},
{cmd:tri}, {cmd:triag}, {cmd:epa}, {cmd:epan}, and {cmd:gau} normalize to the
corresponding full kernel names.

{p 4 4 2}
The Stata default {cmd:bwcheck()} is {cmd:50 + p + 1}, matching the paper and
the documented R interface. For cross-language numerical comparisons, set
{cmd:bwcheck()} explicitly so the support rule is identical across runs and
software versions.

{p 4 4 2}
{cmd:cluster()} supplies a cluster identifier to the bandwidth selector. If
{cmd:vce(hc2)} or {cmd:vce(hc3)} is requested with {cmd:cluster()},
{cmd:vce()} is reset to {cmd:hc1}.
Cluster {cmd:hc1} uses the single finite-sample multiplier
{cmd:((N_h - 1)/(N_h - k))*(G/(G - 1))}. It does not also apply the
non-cluster {cmd:hc1} residual multiplier.

{p 4 4 2}
With {cmd:masspoints(check)} or {cmd:masspoints(adjust)}, the command stores
unique-support diagnostics and warns when mass points are detected.
{cmd:masspoints(check)} also suggests {cmd:masspoints(adjust)}.

{p 4 4 2}
With {cmd:masspoints(off)}, {cmd:r(masspoints)} stores the full sample count
and side counts with {cmd:mass = 0}; it does not compute unique-support counts
or print mass-point warnings.

{p 4 4 2}
With {cmd:kink(on)}, those full counts are also the counts used for kink
bandwidth scaling. Use {cmd:masspoints(check)} or {cmd:masspoints(adjust)}
when repeated support should affect kink scaling through unique-support
diagnostics.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:p()}}local polynomial order; default {cmd:p(1)}{p_end}
{synopt:{cmd:kink()}}{cmd:off} for smooth-boundary distance bandwidths or {cmd:on} for nonsmooth-boundary rates; default {cmd:off}{p_end}
{synopt:{cmd:kernel()}}{cmd:triangular}, {cmd:uniform}, {cmd:epanechnikov}, or {cmd:gaussian}; default {cmd:triangular}{p_end}
{synopt:{cmd:bwselect()}}{cmd:mserd}, {cmd:imserd}, {cmd:msetwo}, or {cmd:imsetwo}; default {cmd:mserd}{p_end}
{synopt:{cmd:vce()}}{cmd:hc0}, {cmd:hc1}, {cmd:hc2}, or {cmd:hc3}; default {cmd:hc1}{p_end}
{synopt:{cmd:bwcheck()}}nonnegative integer minimum preliminary support size; default {cmd:50 + p + 1}{p_end}
{synopt:{cmd:masspoints()}}{cmd:check}, {cmd:adjust}, or {cmd:off}; default {cmd:check}{p_end}
{synopt:{cmd:scaleregul()}}nonnegative regularization scale; default {cmd:1}{p_end}
{synopt:{cmd:cqt()}}quantile fraction for preliminary bias estimation, between 0 and 1; default {cmd:.5}{p_end}
{synopt:{cmd:cluster()}}cluster identifier; HC2/HC3 auto-downgraded to HC1 under clustering{p_end}
{synoptline}

{phang}
{opt cluster(varname)} supplies a cluster identifier to the bandwidth selector
variance constants.  Only {cmd:vce(hc0)} and {cmd:vce(hc1)} are supported
under clustering.  If {cmd:vce(hc2)} or {cmd:vce(hc3)} is specified together
with {cmd:cluster()}, the variance estimator is automatically downgraded to
{cmd:hc1} and a warning is displayed.  The reason is that HC2 and HC3
leverage-based adjustments have prohibitive computational cost and unclear
statistical properties under cluster dependence.
{p_end}

{phang}
{opt scaleregul(#)} controls the scale of the regularization term in the
bandwidth selector.  The Stata default for distance commands is
{cmd:scaleregul(1)}, matching the R package default.  The location commands
({cmd:rdbw2d}, {cmd:rd2d}) default to {cmd:scaleregul(3)}, providing stronger
regularization for the bivariate bandwidth problem.  Users comparing results
across the location and distance commands should note this distinction.
{p_end}

{title:Stored Results}

{p 4 4 2}
{cmd:rdbw2d_dist} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(bws)}}matrix with columns {cmd:b1 b2 h0 h1 Nh0 Nh1}{p_end}
{synopt:{cmd:r(mseconsts)}}matrix with columns {cmd:h0 h1 bias0 bias1 var0 var1 reg0 reg1 Nh0 Nh1 bwmin0 bwmin1 bwmax0 bwmax1}{p_end}
{p 8 8 2}
In {cmd:r(mseconsts)}, {cmd:Nh0} and {cmd:Nh1} are the in-band effective sample sizes used for the selector bandwidth-constant fits. The {cmd:Nh0} and {cmd:Nh1} columns in {cmd:r(bws)} carry those selector counts, while estimation commands report final p-fit counts in {cmd:e(bws)}.{p_end}
{p 8 8 2}
The {cmd:bwmin0}, {cmd:bwmin1}, {cmd:bwmax0}, and {cmd:bwmax1} columns record the side-specific {cmd:bwcheck()} clamp bounds. With {cmd:masspoints(check)} or {cmd:masspoints(off)}, these bounds use the raw repeated support; with {cmd:masspoints(adjust)}, they use unique support. When {cmd:bwcheck(0)} is used, these bound columns are missing.{p_end}
{synopt:{cmd:r(masspoints)}}matrix with columns {cmd:M M0 M1 mass}{p_end}
{synopt:{cmd:r(N)}, {cmd:r(N0)}, {cmd:r(N1)}}sample counts{p_end}
{synopt:{cmd:r(neval)}}number of distance columns / evaluation points{p_end}
{synopt:{cmd:r(p)}, {cmd:r(bwcheck)}}polynomial order and support check used{p_end}
{synopt:{cmd:r(scaleregul)}}scalar: bandwidth regularization scale used{p_end}
{synopt:{cmd:r(cqt)}}scalar: preliminary quantile fraction for kink bias estimation, between 0 and 1{p_end}
{synopt:{cmd:r(kernel)}}kernel function used{p_end}
{synopt:{cmd:r(bwselect)}}bandwidth selector branch used{p_end}
{synopt:{cmd:r(kink)}}string: {cmd:on} or {cmd:off}, recording whether the nonsmooth-boundary undersmoothing path was used{p_end}
{synopt:{cmd:r(vce)}}variance estimator used{p_end}
{synopt:{cmd:r(masspoints_opt)}}mass-points option in effect{p_end}

{p 4 4 2}
The Stata distance selector currently has no boundary-coordinate input, so
{cmd:b1} and {cmd:b2} are missing in {cmd:r(bws)}.

{p 4 4 2}
With multiple evaluation points, {cmd:r(N0)} and {cmd:r(N1)} report the first
point's side counts, while {cmd:r(masspoints)} reports the per-row counts.

{title:Design notes on stored results}

{p 4 4 2}
The following differences between {cmd:rdbw2d_dist} and {cmd:rdbw2d} stored
results reflect intentional design distinctions grounded in the mathematical
structure of the two bandwidth-selection problems.

{p 4 4 2}
{cmd:r(cqt)} and {cmd:r(kink)} are returned by {cmd:rdbw2d_dist} but not by
{cmd:rdbw2d}. The distance selectors support nonsmooth-boundary (kink) inference,
where {cmd:cqt} is the preliminary quantile fraction for bias estimation and
{cmd:kink} records whether the nonsmooth-boundary undersmoothing path was used.
Location selectors do not have a kink path; derivative targeting is achieved
through {cmd:r(deriv)} and {cmd:r(tangvec)}.

{p 4 4 2}
{cmd:r(neval)} is returned by {cmd:rdbw2d_dist} but not by {cmd:rdbw2d}.
{cmd:rdbw2d} accepts exactly one {cmd:at()} coordinate pair, so the number of
evaluation points is implicitly 1. {cmd:rdbw2d_dist} accepts multiple
signed-distance columns, so {cmd:r(neval)} records the number of distance
variables supplied.

{p 4 4 2}
{cmd:r(scalebiascrct)}, {cmd:r(derivsum)}, and {cmd:r(ktype)} are returned by
{cmd:rdbw2d} but not by {cmd:rdbw2d_dist}. Location selectors use a 2D local
polynomial fit where bias correction requires a separate scale parameter to
handle asymmetric bias across the two running dimensions. Distance selectors
use a 1D scalar local polynomial fit, where the bias correction is embedded
directly in the q-order polynomial and does not need a separate scale. Location
selectors also support both product ({cmd:prod}) and radial ({cmd:rad}) kernel
types, while distance selectors use a scalar univariate kernel with no
kernel-type distinction.

{title:Reporting note}

{p 4 4 2}
Use {cmd:rdbw2d_dist} returns to document how automatic signed-distance
bandwidths were chosen. Keep the distance-variable label, sign convention,
{cmd:kink()} choice, selector branch, support rule, and mass-point diagnostics
with the bandwidth record. For a publication row, run {cmd:rd2d_dist} with the
same distance variables and report the final estimate, interval, bandwidths,
and local samples from the estimator returns. Selector bandwidths in
{cmd:r(bws)} explain the choice, but they are not a substitute for the final
estimation row in {cmd:e(results)} and {cmd:e(bws)}. Because the distance
selector has no boundary-coordinate input, {cmd:b1} and {cmd:b2} in
{cmd:r(bws)} are missing by design.

{title:Scope and limitations}

{p 4 4 2}
{cmd:rdbw2d_dist} selects bandwidths for the supplied signed-distance variables
and the chosen {cmd:kink()} convention. It does not choose the scientific
estimand, certify the identifying assumptions of a boundary RD design, or turn
a distance-score row into a boundary-coordinate effect.

{p 4 4 2}
Selector support records and mass-point diagnostics in {cmd:r()} are reporting
aids for documenting the signed-distance bandwidth path before estimation.
They are not identification proofs and should travel with, rather than
replace, the final {cmd:rd2d_dist} estimation row.

{title:Examples}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 400}{p_end}
{phang2}{cmd:. generate double row = mod(_n - 1, 20) - 9.5}{p_end}
{phang2}{cmd:. generate double col = floor((_n - 1) / 20) - 9.5}{p_end}
{phang2}{cmd:. generate double x1 = row / 4}{p_end}
{phang2}{cmd:. generate double x2 = col / 4}{p_end}
{phang2}{cmd:. generate double dist0 = x1}{p_end}
{phang2}{cmd:. generate byte treat = dist0 >= 0}{p_end}
{phang2}{cmd:. generate double y = 1 + .5*x1 - .25*x2 + 1.2*treat + .1*x1*x2 + .05*mod(_n,7)}{p_end}
{phang2}{cmd:. rdbw2d_dist y dist0, bwcheck(20)}{p_end}
{phang2}{cmd:. matrix B = r(bws)}{p_end}
{phang2}{cmd:. assert rowsof(B) == 1}{p_end}
{phang2}{cmd:. assert B[1,3] > 0}{p_end}

{p 4 4 2}
For a generated-example signed-distance selector tutorial in a downloaded
repository copy, see {cmd:rd2d-stata/README.md}. It pairs
{cmd:rd2d-stata/data/data_rd2d.csv} with the generated-example signed-distance
matrix {cmd:rd2d-stata/data/D.csv}. Those CSV files are repository tutorial
inputs, not external empirical data and not files retrieved by {cmd:net get}.
The manuscript provenance file records their checksums and redistribution
limits. The matching selector call is
{cmd:rdbw2d_dist y v1 v5, p(1) kernel(triangular) bwselect(mserd) kink(on)}
with {cmd:bwcheck(22)}, {cmd:masspoints(check)}, {cmd:scaleregul(1)}, and
{cmd:cqt(0.5)}. The returned {cmd:r(bws)} matrix preserves the distance-column
row names {cmd:v1 v5} and reports positive side-specific bandwidths for the
full 20,000-observation generated-example sample.

{p 4 4 2}
Packaged example do-files also create named overview graphs with compact
legends, clear boundary reference lines, and side-specific fitted curves. Use
the location graph to inspect the boundary in the original two-dimensional
running-score space and the distance graph to inspect the zero cutoff in the
supplied signed-distance score.

{p 4 4 2}
The packaged distance example also runs a compact {cmd:kink(on)} selector and
estimator check. It prints a second reporting row with the effective
{cmd:q=p} nonsmooth-boundary path, so users can see how the selector option,
estimator option, returned objects, and reporting row travel together.

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu's distance-based boundary RD source motivates the
signed-distance target, bandwidth-selection path, and nonsmooth-boundary
{cmd:kink(on)} convention used by {cmd:rdbw2d_dist}. Their {cmd:rd2d} software
source describes the broader R, Python, and Stata software contract. Use the
Stata help files and stored results as the command-level reference for syntax,
selector returns, support diagnostics, and reporting.

{p 4 8 2}
For neighboring Stata RD workflows, see the {cmd:rdrobust} and {cmd:rdmulti}
Stata Journal articles. Those sources are scalar-cutoff and multi-cutoff or
multi-score references; they do not replace the signed-distance bandwidth
record documented here.

{title:See Also}

{p 4 4 2}
{help rd2d_dist}, {help rd2d}, {help rdbw2d}
{p_end}
