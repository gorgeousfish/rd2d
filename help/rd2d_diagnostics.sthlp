{smcl}
{* *! version 1.1.0 30jun2026}{...}
{title:rd2d_diagnostics}

{p 4 4 2}
{cmd:rd2d_diagnostics} formats and displays estimation quality diagnostics
stored by {cmd:rd2d} or {cmd:rd2d_dist}. It reads from {cmd:e()} without
modifying any stored estimation object.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_diagnostics}
[{cmd:,}
{cmd:output(}{it:summary}|{it:full}|{it:warnings}{cmd:)}
{cmd:subset(}{it:numlist}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_diagnostics} examines the numerical health of local polynomial fits
at each evaluation point. It reports matrix rank, condition numbers, inversion
method (standard vs. pseudo-inverse fallback), effective sample sizes, and
mass-point ratios, then flags problematic evaluation points with actionable
suggestions.

{p 4 4 2}
Three display modes are available via {cmd:output()}:

{p 8 8 2}
{cmd:summary} (default): a compact table with one row per evaluation point
showing key diagnostics and a status flag ({cmd:ok} or {cmd:WARN}).

{p 8 8 2}
{cmd:full}: a detailed block for each evaluation point, expanding all four
rank/condition/method combinations (p-order control, p-order treated, q-order
control, q-order treated) plus mass-point information.

{p 8 8 2}
{cmd:warnings}: shows only evaluation points with detected issues, each
accompanied by an explanation and remedial suggestion.

{p 4 4 2}
If the current {cmd:e(cmd)} is neither {cmd:rd2d} nor {cmd:rd2d_dist}, the
command exits with error 301.

{title:Options}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:output(}{it:summary}|{it:full}|{it:warnings}{cmd:)}}selects the
display mode; default {cmd:summary}{p_end}
{synopt:{cmd:subset(}{it:numlist}{cmd:)}}restricts display to the specified
evaluation-point indices (integers 1 to {cmd:e(neval)}){p_end}
{synoptline}

{title:Diagnostic indicators}

{p 4 4 2}
The diagnostics matrix {cmd:e(diagnostics)} is {it:neval} x 14 with columns:

{p 8 8 2}
{cmd:b1 b2}: boundary point coordinates

{p 8 8 2}
{cmd:rank_p0 cond_p0 fb_p0}: matrix rank, condition number, and fallback flag
for the control-side polynomial fit at order {it:p}

{p 8 8 2}
{cmd:rank_p1 cond_p1 fb_p1}: same for treated-side, order {it:p}

{p 8 8 2}
{cmd:rank_q0 cond_q0 fb_q0}: same for control-side, order {it:q}

{p 8 8 2}
{cmd:rank_q1 cond_q1 fb_q1}: same for treated-side, order {it:q}

{p 4 4 2}
The fallback flag is 0 for standard inversion ({cmd:invsym}) and 1 for
pseudo-inverse ({cmd:pinv}). The full-rank basis count for a 2D polynomial of
order {it:k} is ({it:k}+1)({it:k}+2)/2.

{title:Warning thresholds}

{p 4 4 2}
The command applies the following thresholds to flag problematic points:

{p 8 8 2}
{bf:Condition number > 1e6}: the design matrix is near-singular. Suggestion:
use the {cmd:stdvars} option to standardize running variables, or increase the
bandwidth.

{p 8 8 2}
{bf:Fallback = pinv}: the design matrix is rank-deficient and a pseudo-inverse
was used instead of the standard Cholesky inversion. Suggestion: increase the
bandwidth or reduce the polynomial order.

{p 8 8 2}
{bf:Effective N < 30}: very few observations fall within the bandwidth window.
Suggestion: increase the bandwidth.

{p 8 8 2}
{bf:Mass-point ratio > 0.2}: more than 20% of the effective sample consists of
mass points. Suggestion: use the {cmd:masspoints(adjust)} option.

{title:Stored results}

{p 4 4 2}
{cmd:rd2d_diagnostics} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(max_cond)}}maximum condition number across all displayed
evaluation points and all four fit sides{p_end}
{synopt:{cmd:r(n_fallback)}}number of evaluation points where at least one
side used pinv fallback{p_end}
{synopt:{cmd:r(min_effN)}}minimum effective sample size (Nh0 or Nh1) across
all displayed points{p_end}
{synopt:{cmd:r(has_warnings)}}1 if any threshold is exceeded, 0
otherwise{p_end}

{title:Examples}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.2) p(1) q(2)}{p_end}
{phang2}{cmd:. rd2d_diagnostics}{p_end}
{phang2}{cmd:. rd2d_diagnostics, output(full)}{p_end}
{phang2}{cmd:. rd2d_diagnostics, output(warnings)}{p_end}
{phang2}{cmd:. rd2d_diagnostics, subset(1 3)}{p_end}

{p 4 4 2}
Programmatic use of return values for automated quality gates:

{phang2}{cmd:. rd2d_diagnostics}{p_end}
{phang2}{cmd:. if r(has_warnings) == 1 {c -(}}{p_end}
{phang2}{cmd:.     di "Estimation has quality issues"}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{p 4 4 2}
After {cmd:rd2d_dist}:

{phang2}{cmd:. rd2d_dist y dist treat, at(0 0 0.3 0.7) h(0.8)}{p_end}
{phang2}{cmd:. rd2d_diagnostics}{p_end}
{phang2}{cmd:. rd2d_diagnostics, output(warnings)}{p_end}

{title:Interpretation guide}

{p 4 4 2}
{bf:Condition number}: measures how sensitive the solution is to small
perturbations in the data. Values below 1e4 are typically safe; values between
1e4 and 1e6 merit caution; values above 1e6 indicate near-singularity.

{p 4 4 2}
{bf:Rank}: the numerical rank of the design matrix. When rank is below the
full-rank basis count, the polynomial cannot be uniquely identified, and the
pseudo-inverse fallback is triggered.

{p 4 4 2}
{bf:Effective N}: the number of observations within the bandwidth window for
each side. Small effective samples inflate variance and may cause numerical
instability.

{p 4 4 2}
{bf:Mass-point ratio}: the proportion of the effective sample at identical
covariate locations. High ratios reduce effective degrees of freedom.

{p 4 4 2}
{bf:CB PSD adjustment}: when the confidence band correlation matrix has
negative eigenvalues (indicating numerical imprecision), eigenvalues are
clamped to zero. The {cmd:cb_psd_adjusted} flag records whether this adjustment
was applied.

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rd2d_summary}, {help rdbw2d}
