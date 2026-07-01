{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d_dist}

{p 4 4 2}
{cmd:rd2d_dist} implements signed-distance local polynomial boundary
regression discontinuity estimation and inference.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_dist} {it:yvar} {it:distvarlist} {ifin}
[{cmd:,}
{cmd:h(}{it:hlist}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:kink(}{it:on|off}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:side(}{it:two|left|right}{cmd:)}
{cmd:nocbands}
{cmd:repp(}{it:#}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:vce(}{it:hc0|hc1|hc2|hc3}{cmd:)}
{cmd:cluster(}{it:varname}{cmd:)}
{cmd:rbc(}{it:on|off}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:check|adjust|off}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:cqt(}{it:#}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_dist} estimates treatment effects using one or more signed-distance
running variables. Negative distances denote control observations and
nonnegative distances denote treated observations.

{p 4 4 2}
The target label for a table row comes from the signed-distance variable name.
The sign convention and {cmd:kink()} choice are part of that target. Because
the command has no boundary-coordinate input, any coordinate printed next to a
distance row is a do-file label, not a command-returned target. Keep that label
separate from the distance-row name, final bandwidths, local samples,
intervals, and diagnostics.

{p 4 4 2}
If {cmd:h()} is omitted, bandwidths are selected by {help rdbw2d_dist}. If
{cmd:h()} is supplied, one positive bandwidth is reused for both sides, two
positive side-specific bandwidths are reused for every distance column, or two
bandwidths are supplied per distance column.

{p 4 4 2}
Kernel aliases follow the R interface where practical: {cmd:uni}, {cmd:unif},
{cmd:tri}, {cmd:triag}, {cmd:epa}, {cmd:epan}, and {cmd:gau} normalize to the
corresponding full kernel names.

{p 4 4 2}
With {cmd:rbc(on)} and {cmd:kink(off)}, the default q-order path uses
{cmd:q(p+1)}. With {cmd:kink(on)}, the command uses the nonsmooth-boundary
undersmoothing convention and resets q-order inference to {cmd:p()} when
needed; it should not be described as robust bias correction. With
{cmd:rbc(off)}, {cmd:q()} is reset to {cmd:p()} and the reported q-order
estimate equals the p-order estimand path.

{p 4 4 2}
{cmd:side()} chooses two-sided or one-sided pointwise confidence intervals.
The same side convention is used for the uniform confidence bands (computed
by default).

{p 4 4 2}
With uniform confidence bands (the default), {cmd:rd2d_dist} fills {cmd:CB_lower} and {cmd:CB_upper} in
{cmd:e(results)} using Gaussian-simulation critical values. The covariance uses
the q-order signed residual influence representation for the distance fits.

{p 4 4 2}
For one-sided bands, {cmd:side(left)} reports a finite upper band with
{cmd:CB_lower = -c(maxdouble)}, and {cmd:side(right)} reports a finite lower
band with {cmd:CB_upper = c(maxdouble)}. The simulated critical value remains
positive in both cases.

{p 4 4 2}
With {cmd:cluster()}, {cmd:rd2d_dist} uses cluster-summed sandwich scores for
pointwise standard errors and the q-order covariance matrix. If {cmd:vce(hc2)}
or {cmd:vce(hc3)} is requested with clustering, {cmd:vce()} is reset to
{cmd:hc1}.
Cluster {cmd:hc1} uses the single finite-sample multiplier
{cmd:((N_h - 1)/(N_h - k))*(G/(G - 1))}. It does not also apply the
non-cluster {cmd:hc1} residual multiplier.

{p 4 4 2}
When {cmd:h()} is omitted, {cmd:rd2d_dist} delegates bandwidth selection to
{cmd:rdbw2d_dist}. Mass-point diagnostics from that selector are still reported
by {cmd:rd2d_dist}: {cmd:masspoints(check)} prints the mass-point warning and
{cmd:masspoints(adjust)} applies the adjusted support while suppressing only
the adjustment suggestion.

{p 4 4 2}
When {cmd:h()} is supplied, {cmd:masspoints(adjust)} does not alter the user
bandwidths. It stores the same unique-support diagnostics as
{cmd:masspoints(check)} and suppresses only the adjustment suggestion.

{p 4 4 2}
After estimation, {cmd:rd2d_dist} prints a compact table with {cmd:Est.q},
{cmd:Se.q}, pointwise confidence intervals, and
{cmd:CB.low}/{cmd:CB.high} uniform-band endpoints for
each distance column. When {cmd:nocbands} is specified, the printed table
shows {cmd:Est.p}, {cmd:Se.p}, {cmd:Est.q}, {cmd:Se.q}, and pointwise confidence
intervals instead. Long distance-column labels
are shortened only in the printed table, using a row suffix to keep labels
distinct. The complete numeric contract and full row names remain in
{cmd:e(results)} and the related matrices listed below. Printed tables adapt to
narrower Stata {cmd:linesize} settings without changing stored matrix names or
columns. At compact {cmd:linesize(50)} widths, estimation tables shorten
endpoint labels to {cmd:CI.lo}/{cmd:CI.hi} and {cmd:CB.lo}/{cmd:CB.hi}.
Open one-sided endpoints are printed as {cmd:-inf} or {cmd:inf}; the stored
matrix sentinels remain {cmd:-c(maxdouble)} and {cmd:c(maxdouble)}.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}

{dlgtab:Core estimation}

{synopt:{cmd:p()}, {cmd:q()}}polynomial orders; default {cmd:p(1)}; {cmd:q(p+1)} with {cmd:rbc(on)} and {cmd:kink(off)}, {cmd:q(p)} with {cmd:kink(on)} or {cmd:rbc(off)}{p_end}
{synopt:{cmd:kink()}}{cmd:off} for smooth-boundary inference or {cmd:on} for kink-adaptive undersmoothing (sets {cmd:q=p}); default {cmd:off}{p_end}
{synopt:{cmd:kernel()}}{cmd:triangular}, {cmd:uniform}, {cmd:epanechnikov}, or {cmd:gaussian}; default {cmd:triangular}{p_end}

{dlgtab:Bandwidth selection}

{synopt:{cmd:h()}}manual scalar, side-specific, or per-distance-column bandwidths; omit for automatic {cmd:rdbw2d_dist} bandwidths{p_end}
{synopt:{cmd:bwselect()}}{cmd:mserd}, {cmd:imserd}, {cmd:msetwo}, or {cmd:imsetwo}; default {cmd:mserd}{p_end}
{synopt:{cmd:scaleregul()}}nonnegative regularization scale; default {cmd:1}{p_end}
{synopt:{cmd:bwcheck()}}nonnegative integer minimum preliminary support size; default {cmd:50 + p + 1}{p_end}
{synopt:{cmd:cqt()}}quantile fraction for preliminary bias estimation, between 0 and 1; default {cmd:.5}{p_end}

{dlgtab:Inference}

{synopt:{cmd:vce()}}{cmd:hc0}, {cmd:hc1}, {cmd:hc2}, or {cmd:hc3}; default {cmd:hc1}{p_end}
{synopt:{cmd:cluster()}}cluster identifier; HC2/HC3 auto-downgraded to HC1 under clustering{p_end}
{synopt:{cmd:level()}}confidence level; default {cmd:95}{p_end}
{synopt:{cmd:side()}}{cmd:two}, {cmd:left}, or {cmd:right} confidence intervals and bands; default {cmd:two}{p_end}
{synopt:{cmd:cbands}}compute uniform confidence bands (the default; retained for backward compatibility){p_end}
{synopt:{cmd:nocbands}}suppress uniform confidence bands{p_end}
{synopt:{cmd:repp()}}Gaussian simulation draws for uniform bands; default {cmd:1000}{p_end}
{synopt:{cmd:rbc()}}{cmd:on} or {cmd:off}; default {cmd:on}{p_end}

{dlgtab:Variable handling}

{synopt:{cmd:masspoints()}}{cmd:check}, {cmd:adjust}, or {cmd:off}; default {cmd:check}{p_end}

{synoptline}

{phang}
{opt kink(on)} enables kink-adaptive inference for boundaries with geometric
corner points (kinks).  At a kink, the conditional expectation of the signed
distance variable loses differentiability, making the standard robust bias
correction path invalid.  When {cmd:kink(on)} is specified, the command
automatically sets {cmd:q = p}, which is equivalent to turning off robust bias
correction and relying on undersmoothing instead.  Even if the user explicitly
specifies {cmd:rbc(on)}, the effective inference path remains {cmd:q = p}
(no higher-order polynomial bias correction) when {cmd:kink(on)} is active.
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

{marker kink_technical}{...}
{title:Technical note: boundary kinks}

{pstd}
A boundary kink is a point where the boundary curve is not differentiable
(the tangent direction changes abruptly).  At such points, the signed-distance
conditional expectation function loses smoothness, and the approximation bias
of local polynomial estimation degrades from the standard O(h^{c -(p+1)}) rate
to O(h), regardless of the polynomial order p (Cattaneo, Titiunik, and Yu,
2026, Journal of Econometrics).{p_end}

{pstd}
The minimax convergence rate at a kink is n^{c -1/4}, slower than the standard
rate n^{c -2(p+1)/(2p+3)} achievable away from kinks.{p_end}

{pstd}
When {cmd:kink(on)} is specified, the command sets q = p (effectively
disabling robust bias correction) and uses undersmoothing to achieve the
minimax rate.  This is the Stata implementation of the kink-adaptive strategy.
The R package achieves the same rate through bandwidth shrinkage to the
N^{c -1/3} scale.{p_end}

{pstd}
The location-based command {cmd:rd2d} is inherently adaptive to boundary
kinks because it estimates directly in the bivariate coordinate space without
relying on the smoothness of the distance transformation.  Users with data
near a known boundary kink may prefer {cmd:rd2d} over {cmd:rd2d_dist}.{p_end}

{title:Stored Results}

{p 4 4 2}
{cmd:rd2d_dist} is {cmd:eclass}. It stores the following:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(results)}}matrix with columns {cmd:b1 b2 Est_p Se_p Est_q Se_q z pvalue CI_lower CI_upper CB_lower CB_upper h0 h1 h0_rbc h1_rbc Nh0 Nh1}{p_end}
{synopt:{cmd:e(results_A0)}}control-side local-fit estimates, side-specific standard errors, bandwidths, and effective sample size{p_end}
{synopt:{cmd:e(results_A1)}}treated-side local-fit estimates, side-specific standard errors, bandwidths, and effective sample size{p_end}
{synopt:{cmd:e(bws)}}estimation bandwidths and final p-fit in-band effective sample sizes; when bandwidths are automatic, {cmd:h0} and {cmd:h1} come from {cmd:rdbw2d_dist}, while {cmd:Nh0} and {cmd:Nh1} describe the estimation fit at those bandwidths{p_end}
{synopt:{cmd:e(diagnostics)}}rank, condition number, and generalized-inverse fallback diagnostics{p_end}
{synopt:{cmd:e(masspoints)}}matrix with columns {cmd:M M0 M1 mass}{p_end}
{synopt:{cmd:e(cov_q)}, {cmd:e(corr_q)}}q-order covariance and correlation matrices{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}q-order coefficient vector and covariance matrix posted for Stata estimation commands{p_end}
{synopt:{cmd:e(N)}, {cmd:e(N0)}, {cmd:e(N1)}}sample counts{p_end}
{synopt:{cmd:e(neval)}}number of distance columns / evaluation points{p_end}
{synopt:{cmd:e(p)}, {cmd:e(q)}}polynomial orders used{p_end}
{synopt:{cmd:e(level)}, {cmd:e(repp)}}confidence level and simulation draws used{p_end}
{synopt:{cmd:e(cb_crit)}}uniform-band critical value (computed by default){p_end}
{synopt:{cmd:e(cb_psd_adjusted)}}scalar: 0 or 1, indicating whether the confidence-band covariance was repaired to positive semidefinite form{p_end}
{synopt:{cmd:e(cb_min_eigen)}}scalar: minimum eigenvalue of the q-order covariance before PSD repair{p_end}
{synopt:{cmd:e(bwcheck)}}minimum preliminary support size used{p_end}
{synopt:{cmd:e(scaleregul)}}scalar: bandwidth regularization scale used{p_end}
{synopt:{cmd:e(cqt)}}scalar: preliminary quantile fraction for kink bias estimation, between 0 and 1{p_end}
{synopt:{cmd:e(kernel)}}kernel function used{p_end}
{synopt:{cmd:e(bwselect)}}bandwidth selector branch used{p_end}
{synopt:{cmd:e(kink)}}string: {cmd:on} or {cmd:off}, recording whether the nonsmooth-boundary undersmoothing path was used{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(rbc)}, {cmd:e(side)}}inference options used{p_end}
{synopt:{cmd:e(masspoints_opt)}}mass-point option in effect{p_end}
{synopt:{cmd:e(bwsource)}}string: {cmd:automatic} if bandwidths were selected by {cmd:rdbw2d_dist} or {cmd:user} if specified via {cmd:h()}{p_end}
{synopt:{cmd:e(cbands)}}{cmd:on} if uniform confidence bands were computed (default); {cmd:off} if {cmd:nocbands} was specified{p_end}
{synopt:{cmd:e(clustered)}, {cmd:e(cluster)}}cluster path flag and cluster variable name{p_end}
{synopt:{cmd:e(fallback)}}string: generalized-inverse method used, {cmd:invsym} or {cmd:pinv}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rd2d_dist}{p_end}

{p 4 4 2}
The standard errors in {cmd:e(results_A0)} and {cmd:e(results_A1)} describe the
separate side-specific local fits. Treatment-effect standard errors are the
{cmd:Se_p} and {cmd:Se_q} columns in {cmd:e(results)}; q-order inference uses
the treatment-effect covariance in {cmd:e(cov_q)} and {cmd:e(V)}.

{p 4 4 2}
The Stata command does not accept boundary coordinate metadata for signed
distances. Therefore the distance-command result matrices include {cmd:b1}
and {cmd:b2} columns, but those entries are missing by design.

{p 4 4 2}
With multiple evaluation points, {cmd:e(N0)} and {cmd:e(N1)} report the first
point's side counts, while {cmd:e(masspoints)} reports the per-row counts.

{p 4 4 2}
With {cmd:masspoints(off)}, {cmd:e(masspoints)} stores the full sample count
and side counts with {cmd:mass = 0}; it does not compute unique-support counts
or print mass-point warnings.

{p 4 4 2}
With automatic bandwidths and {cmd:kink(on)}, those full counts are also the
counts used for kink bandwidth scaling, q-order covariance, and confidence
bands. Use {cmd:masspoints(check)} or {cmd:masspoints(adjust)} when repeated
support should affect the kink scaling through unique-support diagnostics.

{title:Design notes on stored results}

{p 4 4 2}
The following differences between {cmd:rd2d_dist} and {cmd:rd2d} stored
results reflect intentional design distinctions grounded in the mathematical
structure of the two estimation problems.

{p 4 4 2}
{cmd:e(cqt)} and {cmd:e(kink)} are returned by {cmd:rd2d_dist} but not by
{cmd:rd2d}. The distance commands support nonsmooth-boundary (kink) inference,
where {cmd:cqt} is the preliminary quantile fraction for bias estimation and
{cmd:kink} records whether the nonsmooth-boundary undersmoothing path was used.
Location commands achieve derivative targeting through {cmd:e(deriv)} and
{cmd:e(tangvec)} instead.

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

{title:Reporting note}

{p 4 4 2}
For publication tables, keep the distance-variable label, sign convention, and
{cmd:kink()} choice with the estimate. Read estimates, intervals, local
samples, and final estimation bandwidths from the same command run, using
{cmd:e(results)} and {cmd:e(bws)}. Selector returns such as {cmd:r(bws)}
explain the automatic bandwidth choice, but they are not a substitute for the
final estimation row. Because the distance command does not accept boundary
coordinates, any coordinate shown next to a distance row is an application
label added by the do-file, not a command-returned target. Inspect
{cmd:e(diagnostics)}, {cmd:e(masspoints)}, {cmd:e(fallback)}, and any
confidence-band diagnostics before rounding the row for a table.

{p 4 4 2}
For table notes, use {cmd:e(fallback)} for generalized-inverse notes,
{cmd:e(masspoints)} and {cmd:e(masspoints_opt)} for support-rule notes,
{cmd:e(kink)} for the signed-distance kink convention, and
{cmd:e(cb_psd_adjusted)}, {cmd:e(cb_min_eigen)}, and {cmd:e(cb_crit)} for
confidence-band covariance notes. Attach the note to the same row whose
estimate, bandwidths, and local samples were read from {cmd:e(results)} and
{cmd:e(bws)}.

{p 4 4 2}
Closest-distance pooling and boundary-point-specific distance grids should be
reported as distance-score targets. If a do-file supplies one distance column
per grid point, the command-returned row name is still the distance-variable
name. Any boundary coordinate printed next to that row is a script-supplied
label, and it must stay separate from the {cmd:rd2d_dist} returns.

{title:Scope and limitations}

{p 4 4 2}
{cmd:rd2d_dist} estimates the supplied signed-distance representation named by
the distance variables and the chosen {cmd:kink()} convention. The command
does not choose the scientific estimand, certify the identifying assumptions
of a boundary RD design, or turn a distance-score row into a
boundary-coordinate effect.

{p 4 4 2}
Diagnostics in {cmd:e(diagnostics)}, {cmd:e(masspoints)}, {cmd:e(fallback)},
{cmd:e(kink)}, and confidence-band records are reporting aids. They help a
do-file disclose conditioning, support, nonsmooth-boundary convention, and
covariance-repair information for the row being reported, but they are not
identification proofs or substitutes for design justification.

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
{phang2}{cmd:. rd2d_dist y dist0, h(1.25) p(1) q(2) cbands repp(200) bwcheck(20)}{p_end}
{phang2}{cmd:. matrix R = e(results)}{p_end}
{phang2}{cmd:. assert rowsof(R) == 1}{p_end}
{phang2}{cmd:. assert e(N) == 400}{p_end}

{p 4 4 2}
For a generated-example signed-distance tutorial in a downloaded repository
copy, see {cmd:rd2d-stata/README.md}. It pairs {cmd:rd2d-stata/data/data_rd2d.csv}
with the generated-example signed-distance matrix
{cmd:rd2d-stata/data/D.csv}. Those CSV files are repository tutorial inputs,
not external empirical data and not files retrieved by {cmd:net get}. The
manuscript provenance file records their checksums and redistribution limits.
The tutorial then runs
{cmd:rd2d_dist y v1 v5} with {cmd:kink(on)}, {cmd:side(left)}, and
{cmd:cbands} {cmd:repp(760)}. The tutorial verifies the full
20,000-observation generated-example sample, preserves the distance-column row names
{cmd:v1 v5}, and records the effective {cmd:q=p} nonsmooth-boundary path.

{p 4 4 2}
Packaged example do-files also create named overview graphs with compact
legends, clear boundary reference lines, and side-specific fitted curves. Use
the location graph to inspect the boundary in the original two-dimensional
running-score space and the distance graph to inspect the zero cutoff in the
supplied signed-distance score.

{p 4 4 2}
The packaged distance example also prints compact reporting rows after checking
the returned objects. One row uses the default smooth-distance path, and a
second row uses {cmd:kink(on)} and records the effective {cmd:q=p}
nonsmooth-boundary path. These rows are quick checks of reporting, not
substantive application estimates.

{phang}Cluster-robust VCE example:{p_end}
{phang2}{stata "do rd2d_cluster_example.do":rd2d_cluster_example.do}{p_end}

{phang}Confidence bands example:{p_end}
{phang2}{stata "do rd2d_cbands_example.do":rd2d_cbands_example.do}{p_end}

{marker relationships}{...}
{title:Result relationships}

{pstd}
The following relationships hold among stored results:{p_end}

{phang2}
{cmd:e(b)}[j] equals {cmd:e(results)}[j,5], the bias-corrected point estimate
at evaluation point j.{p_end}

{phang2}
{cmd:e(V)} equals {cmd:e(cov_q)} when confidence bands are computed (the default).
When {cmd:nocbands} is specified, {cmd:e(V)} contains the diagonal of the
pointwise variance matrix.{p_end}

{phang2}
{cmd:e(corr_q)}[j,k] = {cmd:e(cov_q)}[j,k] / (Se_q[j] * Se_q[k]), where
Se_q[j] = {cmd:e(results)}[j,6].{p_end}

{phang2}
{cmd:e(diagnostics)}[j,5] > 0 or {cmd:e(diagnostics)}[j,8] > 0 indicates that
evaluation point j used a generalized inverse (pinv) due to rank deficiency.
Estimates at such points may have reduced numerical precision.{p_end}

{phang2}
{cmd:e(bws)}[j,3:4] equals {cmd:e(results)}[j,13:14], the estimation
bandwidths used for point j.{p_end}

{marker troubleshooting}{...}
{title:Common issues}

{pstd}
{cmd:r(198)}: Invalid syntax or option value.  Check option spelling and
parameter types.  Common causes: non-numeric {cmd:h()} values,
{cmd:kink()} combined with unsupported options.{p_end}

{pstd}
{cmd:r(2000)}: No observations satisfy the sample conditions.  Verify that
the data contains observations on both sides of the cutoff (negative and
nonnegative distance values) within the specified {cmd:[if]} and {cmd:[in]}
range.{p_end}

{pstd}
{cmd:r(2001)}: Insufficient observations for estimation.  The local sample
within the bandwidth window is smaller than the minimum required
(q+1 for distance).  Consider increasing {cmd:h()}, reducing {cmd:p()}, or
using distance variables with denser support near zero.{p_end}

{pstd}
{cmd:r(498)}: Numerical computation failure.  Typically caused by a singular
or near-singular design matrix within the local window.  Check for
collinear masspoints ({cmd:masspoints(check)}), reduce {cmd:p()}, or
increase the bandwidth.{p_end}

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu's distance-based boundary RD source motivates the
signed-distance target and the nonsmooth-boundary {cmd:kink(on)} path used by
{cmd:rd2d_dist}. Their {cmd:rd2d} software source describes the broader R,
Python, and Stata software contract. Use the Stata help files and stored
results as the command-level reference for syntax, options, returned objects,
and reporting.

{p 4 8 2}
For neighboring Stata RD workflows, see the {cmd:rdrobust} and {cmd:rdmulti}
Stata Journal articles. Those sources are scalar-cutoff and multi-cutoff or
multi-score references; they do not make a signed-distance row a
boundary-coordinate effect.

{title:See Also}

{p 4 4 2}
{help rdbw2d_dist}, {help rd2d}, {help rdbw2d}
