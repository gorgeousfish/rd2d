{smcl}
{* *! version 1.1.0 30jun2026}{...}
{title:rd2d_bwsens}

{p 4 4 2}
{cmd:rd2d_bwsens} performs bandwidth sensitivity analysis after {cmd:rd2d} or
{cmd:rd2d_dist} estimation. It re-estimates the model over a grid of bandwidth
multipliers to assess robustness of treatment effect estimates, standard errors,
and confidence intervals to bandwidth choice.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_bwsens}
[{cmd:,}
{cmd:grid(}{it:numlist}{cmd:)}
{cmd:subset(}{it:numlist}{cmd:)}
{cmd:plot}
{cmd:saving(}{it:filename}{cmd:)}
{cmd:quiet}
{cmd:format(}{it:string}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_bwsens} is a post-estimation command that systematically varies the
bandwidth around its baseline value (multiplier = 1.0) and re-estimates the
treatment effect at each scaled bandwidth. The multiplier scales all bandwidth
dimensions uniformly: for {cmd:rd2d}, all four bandwidth components
(h01, h02, h11, h12) are scaled; for {cmd:rd2d_dist}, both h0 and h1 are scaled.

{p 4 4 2}
The output table shows, for each multiplier, the robust bias-corrected estimate
(Est.q), its standard error (Se.q), and the pointwise confidence interval. The
baseline result (multiplier = 1.0) is marked with an asterisk for reference.

{p 4 4 2}
This tool is intended for exploratory assessment of sensitivity. Published
results should report estimates at the MSE-optimal or user-specified bandwidth,
not the grid of sensitivity multipliers.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:grid(}{it:numlist}{cmd:)}}specifies the bandwidth multipliers;
default is {cmd:0.5 0.75 0.9 1.0 1.1 1.25 1.5 2.0}{p_end}
{synopt:{cmd:subset(}{it:numlist}{cmd:)}}restricts sensitivity analysis to the
specified evaluation-point indices (integers 1 to {cmd:e(neval)}){p_end}
{synopt:{cmd:plot}}generates a connected-line graph of the point estimates with
confidence intervals across multipliers{p_end}
{synopt:{cmd:saving(}{it:filename}{cmd:)}}exports the sensitivity plot to the
specified file; requires {cmd:plot}{p_end}
{synopt:{cmd:quiet}}suppresses the table display (results are still returned in
{cmd:r()}){p_end}
{synopt:{cmd:format(}{it:string}{cmd:)}}sets the numeric display format;
default is {cmd:%9.4f}{p_end}
{synoptline}

{phang}
{opt grid(numlist)} specifies a set of positive multipliers. Each multiplier {it:m}
produces scaled bandwidths h_scaled = m * h_baseline. A multiplier of 1.0 reproduces
the original estimation. Values below 1.0 test narrower bandwidths (more local,
higher variance); values above 1.0 test wider bandwidths (more bias, lower
variance). All elements must be strictly positive.
{p_end}

{phang}
{opt subset(numlist)} restricts the sensitivity display to specified evaluation
points. When multiple evaluation points exist, this allows focusing on one or a
few points of primary interest. Indices must be positive integers not exceeding
{cmd:e(neval)}.
{p_end}

{phang}
{opt plot} produces a graphical display. The x-axis shows the bandwidth multiplier
and the y-axis shows the treatment effect estimate. The baseline (m = 1) is marked
with a dashed vertical line. Confidence intervals are shown as range caps.
{p_end}

{phang}
{opt saving(filename)} saves the plot to disk. Common formats include {cmd:.png},
{cmd:.pdf}, and {cmd:.eps}. This option is silently ignored if {cmd:plot} is not
specified.
{p_end}

{phang}
{opt quiet} suppresses all table output. Use this when you only need the returned
matrix for further programmatic processing.
{p_end}

{phang}
{opt format(string)} controls the numeric format of displayed estimates (e.g.,
{cmd:%12.6f} for more decimal places). Default is {cmd:%9.4f}.
{p_end}

{title:Stored results}

{p 4 4 2}
{cmd:rd2d_bwsens} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(sens_results)}}matrix of sensitivity results; rows are grid
points, columns are multiplier followed by (Est, Se, CI_lo, CI_hi) for each
displayed evaluation point{p_end}
{synopt:{cmd:r(grid)}}string containing the grid of multipliers used{p_end}
{synopt:{cmd:r(ngrid)}}number of grid points{p_end}
{synopt:{cmd:r(ndisp)}}number of displayed evaluation points{p_end}
{synopt:{cmd:r(cmd)}}estimation command ({cmd:rd2d} or {cmd:rd2d_dist}){p_end}
{synoptline}

{title:Interpretation}

{p 4 4 2}
Stable results across multipliers indicate that the treatment effect estimate is
robust to bandwidth choice. Signs to look for:

{p 8 8 2}
{bf:Stability}: Point estimates remain similar across the grid. The confidence
intervals overlap substantially. This suggests the finding is not sensitive to
the particular bandwidth value.

{p 8 8 2}
{bf:Instability}: Estimates change substantially (e.g., sign reversals or large
magnitude changes) across multipliers. Wide confidence intervals at smaller
multipliers (fewer effective observations) are expected but sign changes warrant
caution.

{p 4 4 2}
Note: Sensitivity analysis does not replace formal bandwidth selection. The
MSE-optimal bandwidth from {cmd:rdbw2d} or {cmd:rdbw2d_dist} remains the
recommended choice for point estimation and inference. This tool is for
supplementary robustness assessment.

{title:Examples}

{p 4 4 2}
Basic usage after {cmd:rd2d}:

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) p(1) q(2) kernel(triangular)}{p_end}
{phang2}{cmd:. rd2d_bwsens}{p_end}

{p 4 4 2}
Custom grid with finer resolution around the baseline:

{phang2}{cmd:. rd2d_bwsens, grid(0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4)}{p_end}

{p 4 4 2}
Focus on a single evaluation point with plot:

{phang2}{cmd:. rd2d_bwsens, subset(1) plot saving(bw_sensitivity.png)}{p_end}

{p 4 4 2}
After {cmd:rd2d_dist}:

{phang2}{cmd:. rd2d_dist y dist1 dist2, p(1) q(2) kernel(epanechnikov)}{p_end}
{phang2}{cmd:. rd2d_bwsens}{p_end}

{p 4 4 2}
Suppress display and work with the returned matrix:

{phang2}{cmd:. rd2d_bwsens, quiet}{p_end}
{phang2}{cmd:. matrix list r(sens_results)}{p_end}

{title:Methodological notes}

{p 4 4 2}
The bandwidth sensitivity analysis follows the common practice in the regression
discontinuity literature of reporting estimates at multiples of the MSE-optimal
bandwidth (see Cattaneo, Idrobo, and Titiunik, 2020, {it:A Practical Introduction}
{it:to Regression Discontinuity Designs}). The default grid spans from half to
double the baseline bandwidth, which covers the range typically examined in
applied work.

{p 4 4 2}
All estimation settings except the bandwidth (kernel, polynomial order, VCE
method, cluster variable, etc.) are held fixed across the grid to isolate the
effect of bandwidth choice on the results.

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rdbw2d}, {help rdbw2d_dist},
{help rd2d_summary}, {help rd2d_plot}
