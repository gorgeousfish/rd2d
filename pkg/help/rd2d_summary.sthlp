{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d_summary}

{p 4 4 2}
{cmd:rd2d_summary} provides structured, publication-quality display of
estimation results stored by {cmd:rd2d} or {cmd:rd2d_dist}. It reads from
{cmd:e()} without modifying any stored object.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_summary}
[{cmd:,}
{cmd:output(}{it:main}|{it:bw}{cmd:)}
{cmd:cbuniform}
{cmd:subset(}{it:numlist}{cmd:)}
{cmd:aate(}{it:numlist}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_summary} formats and displays the matrices {cmd:e(results)} and
{cmd:e(bws)} left in memory by the most recent {cmd:rd2d} or {cmd:rd2d_dist}
estimation. The command is read-only: it does not alter scalars, macros, or
matrices in the {cmd:e()} object.

{p 4 4 2}
Two display modes are available. The default {cmd:output(main)} prints the
estimation results table with point estimates, test statistics, p-values, and
confidence intervals or uniform bands. The alternative {cmd:output(bw)} prints
bandwidth diagnostics with effective sample sizes for each evaluation point.

{p 4 4 2}
If the current {cmd:e(cmd)} is neither {cmd:rd2d} nor {cmd:rd2d_dist}, the
command exits with error 301.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:output(}{it:main}|{it:bw}{cmd:)}}selects the displayed table;
default {cmd:main}{p_end}
{synopt:{cmd:cbuniform}}replaces pointwise confidence intervals with uniform
confidence bands in the last column; requires prior estimation with
{cmd:cbands}{p_end}
{synopt:{cmd:subset(}{it:numlist}{cmd:)}}restricts display to the specified
evaluation-point indices (integers 1 to {cmd:e(neval)}); order and repetition
are preserved{p_end}
{synopt:{cmd:aate(}{it:numlist}{cmd:)}}computes and displays the Aggregated
Average Treatment Effect using the supplied weight vector; requires
{cmd:e(cov_q)} from a prior {cmd:cbands} estimation{p_end}
{synoptline}

{phang}
{opt output(main)} displays the estimation results table. Columns are described
in {it:Output format} below.
{p_end}

{phang}
{opt output(bw)} displays the bandwidth diagnostics table. For {cmd:rd2d},
columns are {cmd:b1 b2 h01 h02 h11 h12 Nh0 Nh1} (boundary point, control and
treated bandwidths by coordinate, effective sample sizes). For {cmd:rd2d_dist},
columns are {cmd:b1 b2 h0 h1 Nh0 Nh1} when boundary coordinates are present,
or {cmd:h0 h1 Nh0 Nh1} otherwise.
{p_end}

{phang}
{opt cbuniform} switches the final column from pointwise confidence intervals
to uniform confidence bands. The command verifies that {cmd:e(cbands)} equals
{cmd:on}; otherwise it exits with error 198.
{p_end}

{phang}
{opt subset(numlist)} filters the displayed rows to the listed indices. Indices
must be positive integers not exceeding {cmd:e(neval)}. The display preserves
the user-supplied order and allows repeated indices, so one may reorder or
duplicate rows for comparison.
{p_end}

{phang}
{opt aate(numlist)} supplies a weight vector whose length must equal the number
of displayed evaluation points (after any {cmd:subset()} restriction). Weights
are normalized internally so they sum to one. The AATE row appears below the
main table. If {cmd:e(cov_q)} is absent, the command exits with error 198.
{p_end}

{title:Output format}

{p 4 4 2}
{bf:Main table} ({cmd:output(main)}):

{p 8 8 2}
When boundary coordinates are available (typical for {cmd:rd2d} and for
{cmd:rd2d_dist} with explicit boundary points):

        {text:====================================================================}
        {text:  ID      b1      b2    Est.       z   P>|z|     95% CI}
        {text:====================================================================}
        {text:   1    0.000   1.000  0.4321   2.1045  0.0353   [ 0.0297,  0.8345]}
        {text:   2    0.500   0.500  0.5678   3.0120  0.0026   [ 0.1983,  0.9373]}
        {text:====================================================================}

{p 8 8 2}
When boundary coordinates are absent ({cmd:rd2d_dist} without stored boundary
points), the {cmd:b1} and {cmd:b2} columns are suppressed.

{p 4 4 2}
The final column shows the pointwise confidence interval at the stored
confidence level. With {cmd:cbuniform}, the label changes to
{cmd:95% Unif. CB} and endpoints come from {cmd:CB_lower}/{cmd:CB_upper} in
{cmd:e(results)}.

{p 4 4 2}
One-sided intervals: when the estimation used {cmd:side(left)}, the lower bound
is displayed as {cmd:-inf}. When {cmd:side(right)}, the upper bound is
displayed as {cmd:inf}. Stored matrix sentinels remain
{cmd:-c(maxdouble)}/{cmd:c(maxdouble)}.

{p 4 4 2}
{bf:Bandwidth table} ({cmd:output(bw)}):

{p 8 8 2}
For {cmd:rd2d} (location-based bandwidths are two-dimensional per side):

        {text:===================================================================================}
        {text:  ID      b1      b2     h01     h02     h11     h12     Nh0     Nh1}
        {text:===================================================================================}
        {text:   1    0.000   1.000   1.250   1.250   1.300   1.300     120     135}
        {text:===================================================================================}

{p 8 8 2}
For {cmd:rd2d_dist} (scalar bandwidths per side):

        {text:=================================================================}
        {text:  ID      b1      b2      h0      h1     Nh0     Nh1}
        {text:=================================================================}
        {text:   1    0.000   1.000   0.800   0.900      95     110}
        {text:=================================================================}

{title:AATE aggregation}

{p 4 4 2}
Given {it:J} displayed evaluation points with estimates
{it:tau_1}, ..., {it:tau_J} and a user-supplied weight vector
{it:w_1}, ..., {it:w_J}, the Aggregated Average Treatment Effect is

{p 8 8 2}
tau_AATE = sum_{j=1}^{J} w_j * tau_q(b_j)

{p 4 4 2}
where weights are normalized to sum to one. The variance is

{p 8 8 2}
Var(tau_AATE) = w' * C_q * w

{p 4 4 2}
where {it:C_q} is the submatrix of {cmd:e(cov_q)} corresponding to the
displayed evaluation points (respecting {cmd:subset()} selection). Inference
uses the Gaussian approximation: z = tau_AATE / se, with a two-sided p-value
and a symmetric confidence interval at the stored confidence level.

{p 4 4 2}
Requirements: (1) {cmd:e(cov_q)} must exist (computed by the prior {cmd:cbands}
estimation); (2) the number of weights must equal the number of displayed
points; (3) weights must sum to a positive value.

{p 4 4 2}
Equal weights produce the unweighted boundary average. Unequal weights allow
boundary-length or application-specific aggregation. The normalized weights
are printed above the AATE results for transparency.

{title:Examples}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.2) p(1) q(2) cbands}{p_end}
{phang2}{cmd:. rd2d_summary}{p_end}
{phang2}{cmd:. rd2d_summary, output(bw)}{p_end}
{phang2}{cmd:. rd2d_summary, cbuniform}{p_end}
{phang2}{cmd:. rd2d_summary, subset(2 1)}{p_end}
{phang2}{cmd:. rd2d_summary, aate(1 1)}{p_end}

{p 4 4 2}
Subset with repetition, useful for comparing a single point against the
aggregate:

{phang2}{cmd:. rd2d_summary, subset(1 1 2)}{p_end}
{phang2}{cmd:. rd2d_summary, subset(2) aate(1)}{p_end}

{p 4 4 2}
After {cmd:rd2d_dist}:

{phang2}{cmd:. rd2d_dist y dist treat, at(0 0 0.3 0.7) h(0.8)}{p_end}
{phang2}{cmd:. rd2d_summary}{p_end}
{phang2}{cmd:. rd2d_summary, output(bw)}{p_end}

{p 4 4 2}
The command does not modify stored results, so multiple calls with different
display options leave the estimation object unchanged for subsequent scripting.

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu's boundary RD method and software sources define the
estimation objects displayed by {cmd:rd2d_summary}. The covariance matrix
{cmd:e(cov_q)} used for AATE inference is the q-order joint covariance
documented in {help rd2d} and {help rd2d_dist}.

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rdbw2d}, {help rdbw2d_dist}
