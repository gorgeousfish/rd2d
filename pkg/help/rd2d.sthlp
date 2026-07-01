{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d}

{p 4 4 2}
{cmd:rd2d} implements location-based bivariate local polynomial boundary
regression discontinuity estimation and inference.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d} {it:yvar} {it:x1var} {it:x2var} {it:tvar} {ifin}
{cmd:,}
{cmd:at(}{it:b1 b2} [{it:b1 b2 ...}]{cmd:)}
[{cmd:h(}{it:hlist}{cmd:)}
{cmd:p(}{it:#}{cmd:)} {cmd:q(}{it:#}{cmd:)}
{cmd:deriv(}{it:d1 d2}{cmd:)} {cmd:tangvec(}{it:v1 v2}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)} {cmd:ktype(}{it:prod|rad}{cmd:)}
{cmd:vce(}{it:hc0|hc1|hc2|hc3}{cmd:)}
{cmd:cluster(}{it:varname}{cmd:)}
{cmd:level(}{it:#}{cmd:)} {cmd:repp(}{it:#}{cmd:)}
{cmd:side(}{it:two|left|right}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)} {cmd:method(}{it:dpi|rot}{cmd:)}
{cmd:rbc(}{it:on|off}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:check|adjust|off}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:scalebiascrct(}{it:#}{cmd:)}
{cmd:stdvars} {cmd:cbands}]

{title:Description}

{p 4 4 2}
{cmd:rd2d} estimates boundary treatment effects at one or more boundary
locations. The four required variables are the outcome, the two running
coordinates, and a binary treatment indicator. Values of {it:tvar} must be
0 for control observations and 1 for treated observations.

{p 4 4 2}
The target label for a table row comes from {cmd:at()}. Each coordinate pair
names a boundary point in the original two-score space, and {cmd:deriv()} or
{cmd:tangvec()} further qualifies the reported target when used. Printed row
names such as {cmd:at1} are compact labels; publication tables should keep the
coordinate, derivative or tangent-vector target, final bandwidths, local
samples, intervals, and diagnostics with the same row.

{p 4 4 2}
If {cmd:h()} is omitted, {cmd:rd2d} uses {help rdbw2d} selector constants for
the requested evaluation points. With {cmd:bwselect(imserd)} or
{cmd:bwselect(imsetwo)}, the command averages the point-specific constants
before forming the integrated bandwidths. If {cmd:h()} is supplied, one
positive bandwidth is reused for all coordinates and sides, four positive
bandwidths are reused for each evaluation point, or four bandwidths are
supplied per evaluation point.

{p 4 4 2}
Kernel aliases follow the R interface where practical: {cmd:uni}, {cmd:unif},
{cmd:tri}, {cmd:triag}, {cmd:epa}, {cmd:epan}, and {cmd:gau} normalize to the
corresponding full kernel names.

{p 4 4 2}
By default the command reports robust-bias-corrected inference using
{cmd:q(p+1)}. With {cmd:rbc(off)}, {cmd:q()} is reset to {cmd:p()} and the
reported q-order estimate equals the p-order estimand path.

{p 4 4 2}
{cmd:side()} chooses two-sided or one-sided pointwise confidence intervals.
When {cmd:cbands} is also specified, the same side convention is used for the
uniform confidence bands.

{p 4 4 2}
After estimation, {cmd:rd2d} prints a compact table with {cmd:Est.p},
{cmd:Se.p}, {cmd:Est.q}, {cmd:Se.q}, and pointwise confidence intervals for
each evaluation point. When {cmd:cbands} is requested, the printed table
focuses on {cmd:Est.q}, {cmd:Se.q}, pointwise confidence intervals, and
{cmd:CB.low}/{cmd:CB.high} uniform-band endpoints. Evaluation points are labeled
{cmd:at1}, {cmd:at2}, and so on in printed output and in the row and column
names of the stored matrices listed below. Printed tables adapt to narrower
Stata {cmd:linesize} settings without changing stored matrix names or columns.
At compact {cmd:linesize(50)} widths, estimation tables shorten endpoint
labels to {cmd:CI.lo}/{cmd:CI.hi} and {cmd:CB.lo}/{cmd:CB.hi}.
Open one-sided endpoints are printed as {cmd:-inf} or {cmd:inf}; the stored
matrix sentinels remain {cmd:-c(maxdouble)} and {cmd:c(maxdouble)}.

{p 4 4 2}
{cmd:deriv()} changes the derivative target in the bivariate local polynomial
basis, including the factorial scaling implied by raw monomial coefficients.
{cmd:tangvec()} estimates a first-order directional derivative and
overrides {cmd:deriv()}; the current implementation supports {cmd:tangvec()}
only with one {cmd:at()} point, and the direction vector must be finite and
nonzero.

{p 4 4 2}
With {cmd:cbands}, the command fills {cmd:CB_lower} and {cmd:CB_upper} in
{cmd:e(results)} using Gaussian-simulation critical values. {cmd:repp()} sets
the number of simulation draws; the default is 1000.

{p 4 4 2}
For one-sided bands, {cmd:side(left)} reports a finite upper band with
{cmd:CB_lower = -c(maxdouble)}, and {cmd:side(right)} reports a finite lower
band with {cmd:CB_upper = c(maxdouble)}. The simulated critical value remains
positive in both cases.

{p 4 4 2}
With {cmd:cluster()}, {cmd:rd2d} uses cluster-summed sandwich scores for
side-specific local polynomial fits. If {cmd:vce(hc2)} or {cmd:vce(hc3)} is
requested with clustering, {cmd:vce()} is reset to {cmd:hc1}.
Cluster {cmd:hc1} uses the single finite-sample multiplier
{cmd:((N_h - 1)/(N_h - k))*(G/(G - 1))}. It does not also apply the
non-cluster {cmd:hc1} residual multiplier.

{p 4 4 2}
When {cmd:h()} is omitted, {cmd:rd2d} delegates bandwidth selection to
{cmd:rdbw2d}. Mass-point diagnostics from that selector are still reported by
{cmd:rd2d}: {cmd:masspoints(check)} prints the mass-point warning and
{cmd:masspoints(adjust)} applies the adjusted support while suppressing only
the adjustment suggestion.

{p 4 4 2}
When {cmd:h()} is supplied, {cmd:masspoints(adjust)} does not alter the user
bandwidths. It stores the same unique-support diagnostics as
{cmd:masspoints(check)} and suppresses only the adjustment suggestion.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{dlgtab:Evaluation points}

{synopt:{cmd:at()}}one or more paired boundary coordinates; required{p_end}

{dlgtab:Core estimation}

{synopt:{cmd:p()}, {cmd:q()}}polynomial orders; default {cmd:p(1)} and {cmd:q(p+1)} with {cmd:rbc(on)}{p_end}
{synopt:{cmd:kernel()}}{cmd:triangular}, {cmd:uniform}, {cmd:epanechnikov}, or {cmd:gaussian}; default {cmd:triangular}{p_end}
{synopt:{cmd:ktype()}}{cmd:prod} for product kernels or {cmd:rad} for radial kernels; default {cmd:prod}{p_end}
{synopt:{cmd:deriv()}}derivative order target; default {cmd:deriv(0 0)}; requires {cmd:d1+d2 <= p}{p_end}
{synopt:{cmd:tangvec()}}direction vector for directional derivatives; overrides {cmd:deriv()}; one {cmd:at()} point only{p_end}

{dlgtab:Bandwidth selection}

{synopt:{cmd:h()}}manual bandwidths; omit for automatic {cmd:rdbw2d} bandwidths{p_end}
{synopt:{cmd:bwselect()}}{cmd:mserd}, {cmd:imserd}, {cmd:msetwo}, or {cmd:imsetwo}; default {cmd:mserd}{p_end}
{synopt:{cmd:method()}}{cmd:dpi} or {cmd:rot} bandwidth method; default {cmd:dpi}{p_end}
{synopt:{cmd:scaleregul()}}nonnegative regularization scale; default {cmd:3}{p_end}
{synopt:{cmd:scalebiascrct()}}nonnegative bias-correction scale; default {cmd:1}{p_end}
{synopt:{cmd:bwcheck()}}nonnegative integer minimum preliminary support size; default {cmd:50 + p + 1}{p_end}

{dlgtab:Inference}

{synopt:{cmd:vce()}}{cmd:hc0}, {cmd:hc1}, {cmd:hc2}, or {cmd:hc3}; default {cmd:hc1}{p_end}
{synopt:{cmd:cluster()}}cluster identifier; HC2/HC3 auto-downgraded to HC1 under clustering{p_end}
{synopt:{cmd:level()}}confidence level; default {cmd:95}{p_end}
{synopt:{cmd:side()}}{cmd:two}, {cmd:left}, or {cmd:right} confidence intervals and bands; default {cmd:two}{p_end}
{synopt:{cmd:cbands}}compute uniform confidence bands{p_end}
{synopt:{cmd:repp()}}Gaussian simulation draws for uniform bands; default {cmd:1000}{p_end}
{synopt:{cmd:rbc()}}{cmd:on} or {cmd:off}; default {cmd:on}{p_end}

{dlgtab:Variable handling}

{synopt:{cmd:stdvars}}standardize coordinates before automatic bandwidth selection{p_end}
{synopt:{cmd:masspoints()}}{cmd:check}, {cmd:adjust}, or {cmd:off}; default {cmd:check}{p_end}

{synoptline}

{phang}
{opt deriv(d1 d2)} specifies the derivative orders with respect to {it:x1} and
{it:x2} in the bivariate local polynomial basis.
Default is {cmd:deriv(0 0)}, which estimates the level difference (treatment
effect) at the boundary point.  The constraint {cmd:d1 + d2 <= p} must hold.
Higher-order derivatives are identified from the local polynomial fit of order
{cmd:p()}.
{p_end}

{phang}
{opt tangvec(v1 v2)} requests estimation of a first-order directional
derivative in the direction ({it:v1}, {it:v2}).  The direction vector must be
finite and nonzero.  When {cmd:tangvec()} is specified, it overrides
{cmd:deriv()}: the derivative target is set internally by the tangent-vector
path regardless of any user-supplied {cmd:deriv()} values.  Only a single
evaluation point (one {cmd:at()} pair) is supported with {cmd:tangvec()}.
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

{phang}
{opt cluster(varname)} requests cluster-robust inference using the groups
defined by {it:varname}.  Only {cmd:vce(hc0)} and {cmd:vce(hc1)} are supported
under clustering.  If {cmd:vce(hc2)} or {cmd:vce(hc3)} is specified together
with {cmd:cluster()}, the variance estimator is automatically downgraded to
{cmd:hc1} and a warning is displayed.  The reason is that HC2 and HC3
leverage-based adjustments have prohibitive computational cost and unclear
statistical properties under cluster dependence.
{p_end}

{phang}
{opt stdvars} standardizes ({it:x1}, {it:x2}) by their sample standard
deviations before bandwidth computation.  This ensures that the regularization
term {cmd:scaleregul()} operates on comparable scales across dimensions.
Without standardization, dimensions with larger variance dominate the bandwidth
selection.  Standardization is recommended when {it:x1} and {it:x2} have
different units or substantially different variances.
{p_end}

{title:Stored Results}

{p 4 4 2}
{cmd:rd2d} is {cmd:eclass}. It stores the following:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(results)}}matrix with columns {cmd:b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h01 h02 h11 h12 Nh0 Nh1}{p_end}
{synopt:{cmd:e(results_A0)}}control-side local-fit estimates, side-specific standard errors, bandwidths, and effective sample size{p_end}
{synopt:{cmd:e(results_A1)}}treated-side local-fit estimates, side-specific standard errors, bandwidths, and effective sample size{p_end}
{synopt:{cmd:e(bws)}}estimation bandwidths and final p-fit in-band effective sample sizes; when bandwidths are automatic, {cmd:h01}, {cmd:h02}, {cmd:h11}, and {cmd:h12} come from {cmd:rdbw2d}, while {cmd:Nh0} and {cmd:Nh1} describe the estimation fit at those bandwidths{p_end}
{synopt:{cmd:e(diagnostics)}}rank, condition number, and generalized-inverse fallback diagnostics{p_end}
{synopt:{cmd:e(masspoints)}}matrix with columns {cmd:M M0 M1 mass}{p_end}
{synopt:{cmd:e(cov_q)}, {cmd:e(corr_q)}}q-order covariance and correlation matrices{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}q-order coefficient vector and covariance matrix posted for Stata estimation commands{p_end}
{synopt:{cmd:e(N)}, {cmd:e(N0)}, {cmd:e(N1)}}sample counts{p_end}
{synopt:{cmd:e(neval)}}number of evaluation points{p_end}
{synopt:{cmd:e(p)}, {cmd:e(q)}}polynomial orders used{p_end}
{synopt:{cmd:e(level)}, {cmd:e(repp)}}confidence level and simulation draws used{p_end}
{synopt:{cmd:e(cb_crit)}}uniform-band critical value when {cmd:cbands} is used{p_end}
{synopt:{cmd:e(cb_psd_adjusted)}}scalar: 0 or 1, indicating whether the confidence-band covariance was repaired to positive semidefinite form{p_end}
{synopt:{cmd:e(cb_min_eigen)}}scalar: minimum eigenvalue of the q-order covariance before PSD repair{p_end}
{synopt:{cmd:e(bwcheck)}}minimum preliminary support size used{p_end}
{synopt:{cmd:e(scaleregul)}}scalar: bandwidth regularization scale used{p_end}
{synopt:{cmd:e(scalebiascrct)}}scalar: bias-correction scale parameter for the 2D polynomial fit{p_end}
{synopt:{cmd:e(derivsum)}}total derivative order used by the estimation target{p_end}
{synopt:{cmd:e(kernel)}, {cmd:e(ktype)}}kernel options used{p_end}
{synopt:{cmd:e(bwselect)}, {cmd:e(method)}}bandwidth selector branch and method used{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(rbc)}, {cmd:e(side)}}inference options used{p_end}
{synopt:{cmd:e(masspoints_opt)}}mass-point option in effect{p_end}
{synopt:{cmd:e(bwsource)}}string: {cmd:automatic} if bandwidths were selected by {cmd:rdbw2d} or {cmd:user} if specified via {cmd:h()}{p_end}
{synopt:{cmd:e(cbands)}}{cmd:on} if uniform confidence bands were requested{p_end}
{synopt:{cmd:e(clustered)}, {cmd:e(cluster)}}cluster path flag and cluster variable name{p_end}
{synopt:{cmd:e(stdvars)}}{cmd:on} if automatic bandwidth selection used standardized coordinates{p_end}
{synopt:{cmd:e(deriv)}, {cmd:e(tangvec)}}derivative target and directional vector used{p_end}
{synopt:{cmd:e(fallback)}}string: generalized-inverse method used, {cmd:invsym} or {cmd:pinv}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rd2d}{p_end}

{p 4 4 2}
The standard errors in {cmd:e(results_A0)} and {cmd:e(results_A1)} describe the
separate side-specific local fits. Treatment-effect standard errors are the
{cmd:Se_p} and {cmd:Se_q} columns in {cmd:e(results)}; q-order inference uses
the treatment-effect covariance in {cmd:e(cov_q)} and {cmd:e(V)}.

{p 4 4 2}
With {cmd:masspoints(off)}, {cmd:e(masspoints)} stores the full sample count
and side counts with {cmd:mass = 0}; it does not compute unique-support counts
or print mass-point warnings.

{title:Design notes on stored results}

{p 4 4 2}
The following differences between {cmd:rd2d} and {cmd:rd2d_dist} stored
results reflect intentional design distinctions grounded in the mathematical
structure of the two estimation problems.

{p 4 4 2}
{cmd:e(scalebiascrct)} and {cmd:e(derivsum)} are returned by {cmd:rd2d} but
not by {cmd:rd2d_dist}. Location commands use a 2D local polynomial fit where
bias correction requires a separate scale parameter to handle asymmetric bias
across the two running dimensions. Distance commands use a 1D scalar local
polynomial fit, where the bias correction is embedded directly in the q-order
polynomial and does not need a separate scale.

{p 4 4 2}
{cmd:e(ktype)} is returned by {cmd:rd2d} but not by {cmd:rd2d_dist}. Location
commands support both product ({cmd:prod}) and radial ({cmd:rad}) kernel types
for the bivariate kernel. Distance commands use a scalar univariate kernel and
have no kernel-type distinction.

{p 4 4 2}
{cmd:e(cqt)} and {cmd:e(kink)} are returned by {cmd:rd2d_dist} but not by
{cmd:rd2d}. The distance commands support nonsmooth-boundary (kink) inference,
where {cmd:cqt} is the preliminary quantile fraction for bias estimation and
{cmd:kink} records whether the nonsmooth-boundary undersmoothing path was used.
Location commands achieve derivative targeting through {cmd:e(deriv)} and
{cmd:e(tangvec)} instead.

{title:Reporting note}

{p 4 4 2}
For publication tables, keep the {cmd:at()} coordinate or derivative target
with the estimate. Read estimates, intervals, local samples, and final
estimation bandwidths from the same command run, using {cmd:e(results)} and
{cmd:e(bws)}. Selector returns such as {cmd:r(bws)} explain the automatic
bandwidth choice, but they are not a substitute for the final estimation row.
Inspect {cmd:e(diagnostics)}, {cmd:e(masspoints)}, {cmd:e(fallback)}, and any
confidence-band diagnostics before rounding the row for a table.

{p 4 4 2}
For table notes, use {cmd:e(fallback)} for generalized-inverse notes,
{cmd:e(masspoints)} and {cmd:e(masspoints_opt)} for support-rule notes, and
{cmd:e(cb_psd_adjusted)}, {cmd:e(cb_min_eigen)}, and {cmd:e(cb_crit)} for
confidence-band covariance notes. Attach the note to the same row whose
estimate, bandwidths, and local samples were read from {cmd:e(results)} and
{cmd:e(bws)}.

{p 4 4 2}
When reporting boundary-point or segment heterogeneity, use one location row
per {cmd:at()} target and keep the target label with that row. When reporting
an aggregate along the boundary, first save the source {cmd:rd2d} rows and then
define the aggregation weights or rule in the do-file. {cmd:rd2d} stores the
location rows used to construct such a table; it does not define
application-specific aggregation weights.

{title:Scope and limitations}

{p 4 4 2}
{cmd:rd2d} estimates the boundary-coordinate representation supplied through
{cmd:at()} and any derivative or tangent-vector options. The command does not
choose the scientific estimand, certify the identifying assumptions of a
boundary RD design, or define aggregation weights along the boundary.

{p 4 4 2}
Diagnostics in {cmd:e(diagnostics)}, {cmd:e(masspoints)}, {cmd:e(fallback)},
and confidence-band records are reporting aids. They help a do-file disclose
conditioning, support, and covariance-repair information for the row being
reported, but they are not identification proofs or substitutes for design
justification.

{title:Examples}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 400}{p_end}
{phang2}{cmd:. generate double row = mod(_n - 1, 20) - 9.5}{p_end}
{phang2}{cmd:. generate double col = floor((_n - 1) / 20) - 9.5}{p_end}
{phang2}{cmd:. generate double x1 = row / 4}{p_end}
{phang2}{cmd:. generate double x2 = col / 4}{p_end}
{phang2}{cmd:. generate byte treat = x1 >= 0}{p_end}
{phang2}{cmd:. generate double y = 1 + .5*x1 - .25*x2 + 1.2*treat + .1*x1*x2 + .05*mod(_n,7)}{p_end}
{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0) h(1.25) p(1) q(2) cbands repp(200) bwcheck(20)}{p_end}
{phang2}{cmd:. matrix R = e(results)}{p_end}
{phang2}{cmd:. assert rowsof(R) == 1}{p_end}
{phang2}{cmd:. assert e(N) == 400}{p_end}

{p 4 4 2}
Packaged example do-files also create named overview graphs with compact
legends, clear boundary reference lines, and side-specific fitted curves. Use
the location graph to inspect the boundary in the original two-dimensional
running-score space and the distance graph to inspect the zero cutoff in the
supplied signed-distance score.

{p 4 4 2}
A downloaded repository copy of {cmd:rd2d-stata/README.md} includes a
generated-example tutorial using {cmd:rd2d-stata/data/data_rd2d.csv}. Those
CSV files are repository tutorial inputs, not external empirical data and not
files retrieved by {cmd:net get}. The manuscript provenance file records their
checksums and redistribution limits. The tutorial
demonstrates {cmd:p(2) q(3) tangvec(1 -0.5)} with {cmd:stdvars},
{cmd:side(right)}, and {cmd:cbands} {cmd:repp(790)}.

{phang}Cluster-robust VCE example:{p_end}
{phang2}{stata "do rd2d_cluster_example.do":rd2d_cluster_example.do}{p_end}

{phang}Confidence bands example:{p_end}
{phang2}{stata "do rd2d_cbands_example.do":rd2d_cbands_example.do}{p_end}

{phang}IMSE bandwidth selection example:{p_end}
{phang2}{stata "do rd2d_imse_example.do":rd2d_imse_example.do}{p_end}

{marker relationships}{...}
{title:Result relationships}

{pstd}
The following relationships hold among stored results:{p_end}

{phang2}
{cmd:e(b)}[j] equals {cmd:e(results)}[j,5], the bias-corrected point estimate
at evaluation point j.{p_end}

{phang2}
{cmd:e(V)} equals {cmd:e(cov_q)} when {cmd:cbands} is specified.
When {cmd:cbands} is not specified, {cmd:e(V)} contains the diagonal of the
pointwise variance matrix.{p_end}

{phang2}
{cmd:e(corr_q)}[j,k] = {cmd:e(cov_q)}[j,k] / (Se_q[j] * Se_q[k]), where
Se_q[j] = {cmd:e(results)}[j,6].{p_end}

{phang2}
{cmd:e(diagnostics)}[j,5] > 0 or {cmd:e(diagnostics)}[j,8] > 0 indicates that
evaluation point j used a generalized inverse (pinv) due to rank deficiency.
Estimates at such points may have reduced numerical precision.{p_end}

{phang2}
{cmd:e(bws)}[j,3:6] equals {cmd:e(results)}[j,13:16], the bandwidths used
for point j.{p_end}

{marker troubleshooting}{...}
{title:Common issues}

{pstd}
{cmd:r(198)}: Invalid syntax or option value.  Check option spelling and
parameter types.  Common causes: non-numeric {cmd:at()} coordinates,
{cmd:deriv()} values exceeding p, zero {cmd:tangvec()} direction.{p_end}

{pstd}
{cmd:r(2000)}: No observations satisfy the sample conditions.  Verify that
the data contains observations on both sides of the boundary within the
specified {cmd:[if]} and {cmd:[in]} range.{p_end}

{pstd}
{cmd:r(2001)}: Insufficient observations for estimation.  The local sample
within the bandwidth window is smaller than the minimum required
((q+1)(q+2)/2 for location).  Consider increasing {cmd:h()}, reducing
{cmd:p()}, or using evaluation points in denser data regions.{p_end}

{pstd}
{cmd:r(498)}: Numerical computation failure.  Typically caused by a singular
or near-singular design matrix within the local window.  Check for
collinear masspoints ({cmd:masspoints(check)}), reduce {cmd:p()}, or
increase the bandwidth.{p_end}

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu's boundary RD location-method source defines the
boundary-coordinate target used by {cmd:rd2d}. Their {cmd:rd2d} software source
describes the broader R, Python, and Stata software contract. Use the Stata
help files and stored results as the command-level reference for syntax,
options, returned objects, and reporting.

{p 4 8 2}
For neighboring Stata RD workflows, see the {cmd:rdrobust} and {cmd:rdmulti}
Stata Journal articles. Those sources are scalar-cutoff and multi-cutoff or
multi-score references; they do not replace the boundary-coordinate reporting
contract documented here.

{title:See Also}

{p 4 4 2}
{help rdbw2d}, {help rd2d_dist}, {help rdbw2d_dist}
