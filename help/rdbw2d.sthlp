{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rdbw2d}

{p 4 4 2}
{cmd:rdbw2d} selects bandwidths for location-based bivariate local polynomial
boundary regression discontinuity designs.

{title:Syntax}

{p 8 8 2}
{cmd:rdbw2d} {it:yvar} {it:x1var} {it:x2var} {it:tvar} {ifin}
{cmd:,} {cmd:at(}{it:b1 b2}{cmd:)}
[{cmd:p(}{it:#}{cmd:)}
{cmd:deriv(}{it:d1 d2}{cmd:)} {cmd:tangvec(}{it:v1 v2}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:ktype(}{it:prod|rad}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:method(}{it:dpi|rot}{cmd:)}
{cmd:vce(}{it:hc0|hc1|hc2|hc3}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:check|adjust|off}{cmd:)}
{cmd:c(}{it:varname}{cmd:)} {cmd:cluster(}{it:varname}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:scalebiascrct(}{it:#}{cmd:)}
{cmd:stdvars}]

{title:Description}

{p 4 4 2}
{cmd:rdbw2d} selects bandwidths for bivariate local polynomial boundary RD
designs with a binary treatment indicator. Values of {it:tvar} must be 0 for
control observations and 1 for treated observations. The command returns
boundary coordinates, control and treatment bandwidths, effective sample sizes,
selector constants, and mass-point diagnostics in {cmd:r()}.

{p 4 4 2}
The target-defining input is {cmd:at()}. The selector row in {cmd:r(bws)}
records bandwidths for that single boundary coordinate, optionally qualified
by {cmd:deriv()} or {cmd:tangvec()}. Store the coordinate, derivative or
tangent-vector target, selector branch, support rule, and mass-point
diagnostics with the selector record before using the bandwidths in
{cmd:rd2d}.

{p 4 4 2}
{cmd:at()} accepts exactly one boundary coordinate pair. Estimation at multiple
points is handled by {help rd2d}, which calls {cmd:rdbw2d} once per point when
automatic bandwidths are requested.

{p 4 4 2}
When {cmd:stdvars} is specified, bandwidths are computed on standardized
coordinates and returned on the original coordinate scale; {cmd:b1} and
{cmd:b2} remain the user-provided {cmd:at()} coordinates.

{p 4 4 2}
{cmd:deriv()} changes the local-polynomial target used by the selector.
{cmd:tangvec()} requests a first-order directional derivative target and
overrides {cmd:deriv()}; the direction vector must be finite and nonzero.

{p 4 4 2}
After selection, {cmd:rdbw2d} prints a compact bandwidth table with the
evaluation coordinate and side-specific coordinate bandwidths. The complete
numeric contract remains in {cmd:r(bws)}, {cmd:r(mseconsts)}, and
{cmd:r(masspoints)}; each returned matrix has one row named {cmd:1}. Printed
tables adapt to narrower Stata {cmd:linesize} settings without changing stored
matrix names or columns. At compact {cmd:linesize(50)} widths, the printed
bandwidth labels stay readable in the same compact table layout.

{p 4 4 2}
Kernel aliases follow the R interface where practical: {cmd:uni}, {cmd:unif},
{cmd:tri}, {cmd:triag}, {cmd:epa}, {cmd:epan}, and {cmd:gau} normalize to the
corresponding full kernel names.

{p 4 4 2}
{cmd:c()} and {cmd:cluster()} supply the cluster identifier to the bandwidth
selector. If {cmd:vce(hc2)} or {cmd:vce(hc3)} is requested with either alias,
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

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:at()}}one boundary coordinate pair; required{p_end}
{synopt:{cmd:p()}}local polynomial order; default {cmd:p(1)}{p_end}
{synopt:{cmd:deriv()}}derivative order target; default {cmd:deriv(0 0)}; requires {cmd:d1+d2 <= p}{p_end}
{synopt:{cmd:tangvec()}}direction vector for directional derivative bandwidths; overrides {cmd:deriv()}{p_end}
{synopt:{cmd:kernel()}}{cmd:triangular}, {cmd:uniform}, {cmd:epanechnikov}, or {cmd:gaussian}; default {cmd:triangular}{p_end}
{synopt:{cmd:ktype()}}{cmd:prod} for product kernels or {cmd:rad} for radial kernels; default {cmd:prod}{p_end}
{synopt:{cmd:bwselect()}}{cmd:mserd}, {cmd:imserd}, {cmd:msetwo}, or {cmd:imsetwo}; default {cmd:mserd}{p_end}
{synopt:{cmd:method()}}{cmd:dpi} or {cmd:rot}; default {cmd:dpi}{p_end}
{synopt:{cmd:vce()}}{cmd:hc0}, {cmd:hc1}, {cmd:hc2}, or {cmd:hc3}; default {cmd:hc1}{p_end}
{synopt:{cmd:bwcheck()}}nonnegative integer minimum preliminary support size; default {cmd:50 + p + 1}{p_end}
{synopt:{cmd:masspoints()}}{cmd:check}, {cmd:adjust}, or {cmd:off}; default {cmd:check}{p_end}
{synopt:{cmd:c()}, {cmd:cluster()}}cluster identifier; HC2/HC3 auto-downgraded to HC1 under clustering{p_end}
{synopt:{cmd:scaleregul()}}nonnegative regularization scale; default {cmd:3}{p_end}
{synopt:{cmd:scalebiascrct()}}nonnegative bias-correction scale; default {cmd:1}{p_end}
{synopt:{cmd:stdvars}}standardize running variables before selecting bandwidths{p_end}
{synoptline}

{phang}
{opt deriv(d1 d2)} specifies the derivative orders with respect to {it:x1} and
{it:x2} used by the bandwidth selector.
Default is {cmd:deriv(0 0)}, which selects bandwidths for the level difference
(treatment effect) at the boundary point.  The constraint {cmd:d1 + d2 <= p}
must hold.
{p_end}

{phang}
{opt tangvec(v1 v2)} requests bandwidth selection for a first-order directional
derivative in the direction ({it:v1}, {it:v2}).  The direction vector must be
finite and nonzero.  When {cmd:tangvec()} is specified, it overrides
{cmd:deriv()}: the selector target is set internally by the tangent-vector
path regardless of any user-supplied {cmd:deriv()} values.
{p_end}

{phang}
{opt c(varname)} or equivalently {opt cluster(varname)} supplies a cluster
identifier to the bandwidth selector variance constants.  Only {cmd:vce(hc0)}
and {cmd:vce(hc1)} are supported under clustering.  If {cmd:vce(hc2)} or
{cmd:vce(hc3)} is specified together with {cmd:c()} or {cmd:cluster()}, the
variance estimator is automatically downgraded to {cmd:hc1} and a warning is
displayed.  The reason is that HC2 and HC3 leverage-based adjustments have
prohibitive computational cost and unclear statistical properties under cluster
dependence.
{p_end}

{phang}
{opt stdvars} standardizes ({it:x1}, {it:x2}) by their sample standard
deviations before bandwidth computation.  This ensures that the regularization
term {cmd:scaleregul()} operates on comparable scales across dimensions.
Without standardization, dimensions with larger variance dominate the bandwidth
selection.  Standardization is recommended when {it:x1} and {it:x2} have
different units or substantially different variances.
{p_end}

{phang}
{opt scaleregul(#)} controls the scale of the regularization term in the
bandwidth selector.  The Stata default {cmd:scaleregul(3)} follows the
recommendation in Cattaneo, Titiunik, and Yu (2025).  The R package currently
defaults to {cmd:scaleregul=1} (a known discrepancy).  To reproduce R package
results in Stata, specify {cmd:scaleregul(1)} explicitly.  The
{cmd:scaleregul(3)} default provides stronger regularization and is recommended
for general use.
{p_end}

{title:Stored Results}

{p 4 4 2}
{cmd:rdbw2d} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(bws)}}matrix with columns {cmd:b1 b2 h01 h02 h11 h12}{p_end}
{synopt:{cmd:r(mseconsts)}}matrix with columns {cmd:Nh.0 Nh.1 bias.0 bias.1 var.0 var.1 reg.bias.0 reg.bias.1 reg.var.0 reg.var.1}{p_end}
{p 8 8 2}
In {cmd:r(mseconsts)}, {cmd:Nh.0} and {cmd:Nh.1} are the in-band effective sample sizes used for the selector bandwidth-constant fits, not estimation-command final p-fit counts.{p_end}
{p 8 8 2}
The {cmd:bwcheck()} clamp affects the selected bandwidths in {cmd:r(bws)}, using raw support with {cmd:masspoints(check)} or {cmd:masspoints(off)} and unique support with {cmd:masspoints(adjust)}. Location selectors do not return separate {cmd:bwmin} or {cmd:bwmax} stored columns.{p_end}
{synopt:{cmd:r(masspoints)}}matrix with columns {cmd:M M0 M1 mass}{p_end}
{synopt:{cmd:r(N)}, {cmd:r(N0)}, {cmd:r(N1)}}sample counts{p_end}
{synopt:{cmd:r(p)}, {cmd:r(bwcheck)}}polynomial order and support check used{p_end}
{synopt:{cmd:r(scaleregul)}}scalar: bandwidth regularization scale used{p_end}
{synopt:{cmd:r(scalebiascrct)}}scalar: bias-correction scale parameter for the 2D polynomial fit{p_end}
{synopt:{cmd:r(kernel)}, {cmd:r(ktype)}}kernel options{p_end}
{synopt:{cmd:r(bwselect)}, {cmd:r(method)}}selector branch and method{p_end}
{synopt:{cmd:r(vce)}}variance estimator used{p_end}
{synopt:{cmd:r(masspoints_opt)}}mass-points option in effect{p_end}
{synopt:{cmd:r(stdvars)}}{cmd:on} if standardized coordinates were used{p_end}
{synopt:{cmd:r(deriv)}}string: derivative target orders {it:d1 d2} used by the selector{p_end}
{synopt:{cmd:r(tangvec)}}string: direction vector {it:v1 v2} for directional-derivative bandwidths, or empty{p_end}
{synopt:{cmd:r(derivsum)}}total derivative order used by the selector target{p_end}

{title:Design notes on stored results}

{p 4 4 2}
The following differences between {cmd:rdbw2d} and {cmd:rdbw2d_dist} stored
results reflect intentional design distinctions grounded in the mathematical
structure of the two bandwidth-selection problems.

{p 4 4 2}
{cmd:r(scalebiascrct)} and {cmd:r(derivsum)} are returned by {cmd:rdbw2d} but
not by {cmd:rdbw2d_dist}. Location selectors use a 2D local polynomial fit
where bias correction requires a separate scale parameter to handle asymmetric
bias across the two running dimensions. Distance selectors use a 1D scalar local
polynomial fit, where the bias correction is embedded directly in the q-order
polynomial and does not need a separate scale.

{p 4 4 2}
{cmd:r(ktype)} is returned by {cmd:rdbw2d} but not by {cmd:rdbw2d_dist}. Location
selectors support both product ({cmd:prod}) and radial ({cmd:rad}) kernel types
for the bivariate kernel. Distance selectors use a scalar univariate kernel and
have no kernel-type distinction.

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
evaluation points is implicitly 1. {cmd:rdbw2d_dist} accepts multiple signed-distance
columns, so {cmd:r(neval)} records the number of distance variables supplied.

{title:Reporting note}

{p 4 4 2}
Use {cmd:rdbw2d} returns to document how automatic location bandwidths were
chosen. Keep the {cmd:at()} coordinate, derivative or directional target,
selector branch, support rule, and mass-point diagnostics with the bandwidth
record. For a publication row, run {cmd:rd2d} with the same target and report
the final estimate, interval, bandwidths, and local samples from the estimator
returns. Selector bandwidths in {cmd:r(bws)} explain the choice, but they are
not a substitute for the final estimation row in {cmd:e(results)} and
{cmd:e(bws)}.

{title:Scope and limitations}

{p 4 4 2}
{cmd:rdbw2d} selects bandwidths for the single boundary-coordinate target
supplied through {cmd:at()} and any derivative or tangent-vector options. It
does not choose the scientific estimand, certify the identifying assumptions of
a boundary RD design, or define aggregation weights along the boundary.

{p 4 4 2}
Selector support records and mass-point diagnostics in {cmd:r()} are reporting
aids for documenting the bandwidth path before estimation. They are not
identification proofs and should travel with, rather than replace, the final
{cmd:rd2d} estimation row.

{title:Examples}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 400}{p_end}
{phang2}{cmd:. generate double row = mod(_n - 1, 20) - 9.5}{p_end}
{phang2}{cmd:. generate double col = floor((_n - 1) / 20) - 9.5}{p_end}
{phang2}{cmd:. generate double x1 = row / 4}{p_end}
{phang2}{cmd:. generate double x2 = col / 4}{p_end}
{phang2}{cmd:. generate byte treat = x1 >= 0}{p_end}
{phang2}{cmd:. generate double y = 1 + .5*x1 - .25*x2 + 1.2*treat + .1*x1*x2 + .05*mod(_n,7)}{p_end}
{phang2}{cmd:. rdbw2d y x1 x2 treat, at(0 0) method(rot) bwcheck(20)}{p_end}
{phang2}{cmd:. matrix B = r(bws)}{p_end}
{phang2}{cmd:. assert rowsof(B) == 1}{p_end}
{phang2}{cmd:. assert B[1,3] > 0}{p_end}

{p 4 4 2}
For the generated-example location tutorial in {cmd:rd2d-stata/README.md},
the repository CSV files are tutorial inputs, not external empirical data and
not files retrieved by {cmd:net get}. The manuscript provenance file records
their checksums and redistribution limits. For that tutorial,
the matching selector call is
{cmd:rdbw2d y x1 x2 treat, at(0 50) p(2) tangvec(1 -0.5) kernel(epanechnikov) ktype(prod) bwselect(msetwo) method(dpi) vce(hc3) bwcheck(52) masspoints(check) scaleregul(3) scalebiascrct(1) stdvars}
after loading {cmd:rd2d-stata/data/data_rd2d.csv}. The returned
{cmd:r(bws)} matrix keeps {cmd:at(0 50)} on the user scale, reports
{cmd:r(stdvars) = on}, and returns positive raw-scale coordinate bandwidths
on the full generated-example {cmd:N=20000}, {cmd:N0=6191},
{cmd:N1=13809} sample.

{p 4 4 2}
Packaged example do-files also create named overview graphs with compact
legends, clear boundary reference lines, and side-specific fitted curves. Use
the location graph to inspect the boundary in the original two-dimensional
running-score space and the distance graph to inspect the zero cutoff in the
supplied signed-distance score.

{marker bwselect_guide}{...}
{title:Bandwidth selection strategies}

{pstd}
{cmd:mserd} minimizes the pointwise mean squared error (MSE) at each
evaluation point independently.  The resulting bandwidths vary across
points, adapting to local data density and curvature.{p_end}

{pstd}
{cmd:imserd} minimizes the integrated MSE across all evaluation points,
producing a single common bandwidth.  This is appropriate when the goal is
aggregation (WBATE, AATE) or uniform inference via confidence bands.{p_end}

{pstd}
{cmd:msetwo} and {cmd:imsetwo} allow separate bandwidths for the control
and treated sides, accommodating asymmetric data distributions around the
boundary.{p_end}

{pstd}
For single-point reporting, {cmd:mserd} is appropriate.  For multi-point
aggregation (WBATE, AATE, LBATE) or uniform confidence bands ({cmd:cbands}),
{cmd:imserd} or {cmd:imsetwo} controls the global bias by selecting a
common bandwidth optimized for the integrated criterion.{p_end}

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu's boundary RD location-method source defines the
boundary-coordinate target and bandwidth-selection problem used by
{cmd:rdbw2d}. Their {cmd:rd2d} software source describes the broader R, Python,
and Stata software contract. Use the Stata help files and stored results as the
command-level reference for syntax, selector returns, support diagnostics, and
reporting.

{p 4 8 2}
For neighboring Stata RD workflows, see the {cmd:rdrobust} and {cmd:rdmulti}
Stata Journal articles. Those sources are scalar-cutoff and multi-cutoff or
multi-score references; they do not replace the boundary-coordinate bandwidth
record documented here.

{title:See Also}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rdbw2d_dist}
