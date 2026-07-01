{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d_aggregate}

{p 4 4 2}
{cmd:rd2d_aggregate} computes aggregated treatment effect summaries from
boundary-point estimates stored by {cmd:rd2d} or {cmd:rd2d_dist}.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_aggregate}
{cmd:,}
{cmd:method(}{it:wbate|aate|lbate}{cmd:)}
[{cmd:weights(}{it:numlist}{cmd:)}
{cmd:subset(}{it:numlist}{cmd:)}
{cmd:estimate(}{it:q|p}{cmd:)}
{cmd:level(}{it:#}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_aggregate} is a postestimation command that aggregates boundary
treatment effect estimates across evaluation points. It operates on the
estimation results left in memory by {cmd:rd2d} or {cmd:rd2d_dist} and returns
a scalar summary with inference.

{p 4 4 2}
Three aggregation methods are available. WBATE and AATE produce weighted
averages with normal-approximation inference derived from the stored
covariance matrix {cmd:e(cov_q)}. LBATE identifies the largest pointwise
estimate and constructs a confidence interval from the uniform band critical
value {cmd:e(cb_crit)}.

{p 4 4 2}
The command requires that the preceding {cmd:rd2d} or {cmd:rd2d_dist} call
used the {cmd:cbands} option so that the cross-point covariance or uniform
band critical value is available.

{p 4 4 2}
When {cmd:subset()} is specified, only the designated rows of
{cmd:e(results)} participate in the aggregation. The covariance submatrix is
extracted accordingly. This permits aggregation over scientifically meaningful
subsets of the boundary without re-estimation.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:method()}}aggregation method; required. One of {cmd:wbate},
{cmd:aate}, or {cmd:lbate}{p_end}
{synopt:{cmd:weights()}}nonneg. weights for WBATE, one per selected evaluation
point; required when {cmd:method(wbate)}{p_end}
{synopt:{cmd:subset()}}positive integer indices selecting a subset of
evaluation points for aggregation{p_end}
{synopt:{cmd:estimate()}}use {cmd:q}-order (RBC, default) or {cmd:p}-order
estimates; with {cmd:estimate(p)} only the point estimate is reported.
Not available with {cmd:method(lbate)} (Theorem 6 CI requires q-order
uniform band critical value){p_end}
{synopt:{cmd:level()}}confidence level in percent; default {cmd:95}{p_end}
{synoptline}

{phang}
{opt method(wbate|aate|lbate)} selects the aggregation rule. {cmd:wbate}
computes a user-weighted boundary average. {cmd:aate} uses equal weights
1/J. {cmd:lbate} reports the supremum over evaluation points.
{p_end}

{phang}
{opt weights(numlist)} supplies one nonneg. value per evaluation point (or per
selected subset point). The command normalizes to sum one internally. Required
for {cmd:method(wbate)}; ignored otherwise.
{p_end}

{phang}
{opt subset(numlist)} restricts aggregation to specified rows of
{cmd:e(results)}. Indices must be positive integers between 1 and
{cmd:e(neval)}.
{p_end}

{phang}
{opt estimate(q|p)} controls which column of {cmd:e(results)} supplies the
point estimates. Default is {cmd:q} (robust bias-corrected). With
{cmd:estimate(p)}, inference is suppressed because {cmd:e(cov_q)} covers only
q-order estimates.
{p_end}

{phang}
{opt level(#)} sets the confidence level for intervals. Default is 95. Must
be strictly between 0 and 100.
{p_end}

{title:Methods}

{p 4 4 2}
Let J denote the number of selected evaluation points and let tau_q(b_j)
denote the q-order boundary treatment effect estimate at evaluation point j.
Let SE_q(b_j) denote the corresponding standard error from the q-order
local polynomial fit.

{p 4 4 2}
{ul:WBATE: Weighted Boundary Average Treatment Effect}

{p 8 4 2}
Point estimate:

{p 12 4 2}
tau_WBATE = sum_{j=1}^{J} w_j * tau_q(b_j)

{p 8 4 2}
where w_j >= 0 are user-supplied weights normalized so that
sum_{j=1}^{J} w_j = 1.

{p 8 4 2}
Variance:

{p 12 4 2}
Var(tau_WBATE) = w' * Cov_q * w

{p 8 4 2}
where Cov_q is the q-order covariance matrix stored in {cmd:e(cov_q)} (or
its subset rows and columns when {cmd:subset()} is used).

{p 8 4 2}
Inference:

{p 12 4 2}
SE = sqrt(Var),  z = tau_WBATE / SE,  p = 2*(1 - Phi(|z|))

{p 12 4 2}
CI = [tau_WBATE - z_{alpha/2} * SE,  tau_WBATE + z_{alpha/2} * SE]

{p 8 4 2}
Theoretical basis: Cattaneo, Titiunik, and Yu (2025), Theorems 4 and 5.

{p 4 4 2}
{ul:AATE: Aggregated Average Treatment Effect}

{p 8 4 2}
Point estimate:

{p 12 4 2}
tau_AATE = (1/J) * sum_{j=1}^{J} tau_q(b_j)

{p 8 4 2}
AATE is the equal-weight special case of WBATE with w_j = 1/J for all j.
Variance, standard error, z-statistic, p-value, and confidence interval
follow identically from the WBATE formulae with uniform weights.

{p 4 4 2}
{ul:LBATE: Largest Boundary Average Treatment Effect}

{p 8 4 2}
Point estimate:

{p 12 4 2}
tau_LBATE = max_{j=1,...,J} tau_q(b_j)

{p 8 4 2}
Confidence interval:

{p 12 4 2}
CI_lower = max_{j} [ tau_q(b_j) - q_alpha * SE_q(b_j) ]

{p 12 4 2}
CI_upper = max_{j} [ tau_q(b_j) + q_alpha * SE_q(b_j) ]

{p 8 4 2}
where q_alpha is the uniform confidence band critical value stored in
{cmd:e(cb_crit)}. This construction inverts the uniform band to obtain valid
simultaneous coverage for the supremum functional.

{p 8 4 2}
No z-statistic or p-value is reported. The supremum of a Gaussian process has
no simple normal approximation, so pointwise z-inference is not meaningful
for the largest-effect functional.

{p 8 4 2}
The CI construction exploits the fact that the uniform band covers all J
pointwise estimates simultaneously, so the band inversion for the supremum
functional inherits the same nominal coverage probability.

{p 8 4 2}
LBATE requires the preceding estimation to have used {cmd:cbands}.

{p 8 4 2}
Theoretical basis: Cattaneo, Titiunik, and Yu (2025), Theorem 6.

{title:Stored results}

{p 4 4 2}
{cmd:rd2d_aggregate} is {cmd:rclass}. Stored scalars and locals depend on
the method.

{p 4 4 2}
{ul:WBATE and AATE (with estimate(q)):}

{p 4 4 2}
When {cmd:estimate(q)} is in effect (the default), all variance-based returns
are available. The covariance is extracted from {cmd:e(cov_q)}, which the
preceding {cmd:cbands} call populated.

{synoptset 18 tabbed}{...}
{synopt:{cmd:r(estimate)}}aggregated point estimate{p_end}
{synopt:{cmd:r(se)}}standard error{p_end}
{synopt:{cmd:r(z)}}z-statistic{p_end}
{synopt:{cmd:r(p)}}two-sided p-value{p_end}
{synopt:{cmd:r(ci_lower)}}lower confidence limit{p_end}
{synopt:{cmd:r(ci_upper)}}upper confidence limit{p_end}
{synopt:{cmd:r(variance)}}estimated variance of the aggregate{p_end}
{synopt:{cmd:r(J)}}number of evaluation points used{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(method)}}string: {cmd:wbate} or {cmd:aate}{p_end}
{synopt:{cmd:r(estimate_order)}}string: {cmd:q} or {cmd:p}{p_end}

{p 4 4 2}
With {cmd:estimate(p)}, only {cmd:r(estimate)}, {cmd:r(J)}, {cmd:r(level)},
{cmd:r(method)}, and {cmd:r(estimate_order)} are stored. Variance-based
returns are absent because {cmd:e(cov_q)} does not cover p-order estimates.

{p 4 4 2}
{ul:LBATE:}

{synoptset 18 tabbed}{...}
{synopt:{cmd:r(estimate)}}largest point estimate{p_end}
{synopt:{cmd:r(ci_lower)}}lower confidence limit (uniform band inversion){p_end}
{synopt:{cmd:r(ci_upper)}}upper confidence limit (uniform band inversion){p_end}
{synopt:{cmd:r(cb_crit)}}uniform band critical value used{p_end}
{synopt:{cmd:r(jmax)}}row index of the maximizing evaluation point{p_end}
{synopt:{cmd:r(J)}}number of evaluation points used{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(method)}}string: {cmd:lbate}{p_end}
{synopt:{cmd:r(estimate_order)}}string: {cmd:q} or {cmd:p}{p_end}

{title:Examples}

{pstd}Estimate with three boundary points and compute AATE:{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0  0 1  0 -1) h(1.5) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_aggregate, method(aate)}{p_end}
{phang2}{cmd:. return list}{p_end}

{pstd}User-weighted WBATE with a subset of points:{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0  0 1  0 -1  1 0) h(1.5) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_aggregate, method(wbate) weights(2 1 1) subset(1 2 3)}{p_end}

{pstd}LBATE from distance estimation:{p_end}

{phang2}{cmd:. rd2d_dist y dist treat, at(0 0  0 1  0 -1) h(1.2) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_aggregate, method(lbate)}{p_end}

{pstd}P-order point estimate only (no inference):{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0  0 1) h(1.5) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_aggregate, method(aate) estimate(p)}{p_end}

{pstd}Custom confidence level:{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0  0 1  0 -1) h(1.5) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_aggregate, method(wbate) weights(1 1 1) level(90)}{p_end}

{title:References}

{p 4 8 2}
Cattaneo, M. D., R. Titiunik, and R. Yu. 2025. Boundary discontinuity
regression design. Working paper.
{p_end}

{p 4 8 2}
Theorems 4 and 5 establish the joint asymptotic normality of pointwise
boundary treatment effects and the validity of weighted linear combinations
under the covariance structure provided by the bivariate local polynomial
fit. Theorem 6 provides uniform inference for the supremum functional using
Gaussian simulation critical values, enabling the LBATE confidence set.
{p_end}

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rdbw2d}, {help rdbw2d_dist}

{p 4 4 2}
For boundary treatment effect estimation, see {help rd2d} (location method)
and {help rd2d_dist} (distance method). For automatic bandwidth selection,
see {help rdbw2d} and {help rdbw2d_dist}.
