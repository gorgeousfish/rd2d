{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d_plot}

{p 4 4 2}
{cmd:rd2d_plot} produces post-estimation graphs for boundary regression
discontinuity treatment effect estimates stored by {cmd:rd2d} or {cmd:rd2d_dist}.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_plot}
[{cmd:,}
{cmd:type(}{it:effect|heat}{cmd:)}
{cmd:estimate(}{it:q|p}{cmd:)}
{cmd:interval(}{it:ci|cb|both|none}{cmd:)}
{cmd:subset(}{it:numlist}{cmd:)}
{cmd:label}
{cmd:saving(}{it:filename}{cmd:)}
{cmd:title(}{it:string}{cmd:)}
{cmd:subtitle(}{it:string}{cmd:)}
{cmd:scheme(}{it:schemename}{cmd:)}
{cmd:name(}{it:name}{cmd:)}
{cmd:replace}
{cmd:nodraw}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_plot} visualizes boundary treatment effect estimates from a preceding
{cmd:rd2d} or {cmd:rd2d_dist} estimation command. The command reads the stored
{cmd:e(results)} matrix and produces either an effect plot (point estimates with
uncertainty intervals indexed by evaluation point) or a heat map (bubble scatter
in the two-dimensional boundary coordinate space).

{p 4 4 2}
Two estimand paths are available for display: the bias-corrected q-order
estimate, which is the default inferential target in {cmd:rd2d}, and the
conventional p-order estimate. Interval overlays include pointwise confidence
intervals, uniform confidence bands, or both simultaneously.

{p 4 4 2}
The command operates entirely on stored {cmd:e()} results and does not
re-estimate. It preserves the current dataset via {cmd:preserve}/{cmd:restore}.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:type()}}plot type; default {cmd:effect}{p_end}
{synopt:{cmd:estimate()}}estimand path; default {cmd:q}{p_end}
{synopt:{cmd:interval()}}interval overlay; default {cmd:ci}{p_end}
{synopt:{cmd:subset()}}evaluation-point indices to plot{p_end}
{synopt:{cmd:label}}label x-axis with original point indices{p_end}
{synopt:{cmd:saving()}}save graph to file{p_end}
{synopt:{cmd:title()}}override default title{p_end}
{synopt:{cmd:subtitle()}}override default subtitle{p_end}
{synopt:{cmd:scheme()}}Stata graph scheme{p_end}
{synopt:{cmd:name()}}graph window name{p_end}
{synopt:{cmd:replace}}allow overwriting of saved file or named graph{p_end}
{synopt:{cmd:nodraw}}suppress graph display{p_end}
{synoptline}

{phang}
{opt type(effect|heat)} selects the graph geometry.

{p 8 8 2}
{cmd:effect} produces a coefficient plot: evaluation points on the horizontal
axis, point estimates as markers, and uncertainty intervals as vertical bars or
shaded regions. The zero reference line is drawn as a dashed rule. This is the
natural display for inspecting boundary heterogeneity across evaluation points.

{p 8 8 2}
{cmd:heat} produces a bubble scatter in the ({it:b1}, {it:b2}) coordinate
plane. Marker area is proportional to the absolute value of the treatment
effect; positive effects are rendered in blue, negative effects in orange.
This display requires that the stored results contain nondegenerate boundary
coordinates. For {cmd:rd2d_dist} results that carry only distance information,
{cmd:type(heat)} exits with an error.
{p_end}

{phang}
{opt estimate(q|p)} selects the estimand path.

{p 8 8 2}
{cmd:q} (default) displays the robust bias-corrected estimate and its
associated inferential quantities. This is the recommended reporting path for
boundary RD inference under the Cattaneo, Titiunik, and Yu (2025) framework,
where the q-order polynomial removes first-order smoothing bias.

{p 8 8 2}
{cmd:p} displays the conventional p-order polynomial estimate with pointwise
confidence intervals computed from {cmd:Se_p}. Confidence bands are not defined
for the conventional path; if {cmd:interval(cb)} or {cmd:interval(both)} is
combined with {cmd:estimate(p)}, the interval is silently changed to {cmd:ci}.
{p_end}

{phang}
{opt interval(ci|cb|both|none)} controls uncertainty overlays.

{p 8 8 2}
{cmd:ci} (default) draws pointwise confidence intervals at the estimation
level stored in {cmd:e(level)}. For two-sided intervals, capped range bars
span the lower and upper endpoints. For one-sided intervals, a spike extends
from the estimate toward the finite endpoint.

{p 8 8 2}
{cmd:cb} draws uniform confidence bands as a shaded region. This option
requires that the preceding estimation command was run with {cmd:cbands};
otherwise {cmd:rd2d_plot} exits with an error directing the user to re-estimate
or select {cmd:interval(ci)}.

{p 8 8 2}
{cmd:both} overlays pointwise intervals and uniform bands simultaneously. The
band region appears behind the interval bars, allowing visual comparison of
pointwise and uniform coverage.

{p 8 8 2}
{cmd:none} suppresses all interval overlays and shows only point estimates.
{p_end}

{phang}
{opt subset(numlist)} restricts the plot to a subset of evaluation points. Each
integer in the numlist indexes a row in {cmd:e(results)}. Values exceeding the
number of evaluation points stored in {cmd:e(neval)} produce an error. The
original row indices are preserved as axis labels, maintaining correspondence
with the printed estimation table.
{p_end}

{phang}
{opt label} requests that the x-axis tick marks display the original
evaluation-point indices. For {cmd:type(effect)} this is the default behavior;
the option is relevant only when combined with other customizations.
{p_end}

{phang}
{opt saving(filename)} saves the graph to {it:filename}. When combined with
{cmd:replace}, an existing file is overwritten.
{p_end}

{phang}
{opt title(string)} overrides the automatic title. Without this option, the
default title is "Bias-corrected treatment effect estimates" for
{cmd:estimate(q)} or "Conventional treatment effect estimates" for
{cmd:estimate(p)} in effect plots, and "Treatment effects at boundary" for
heat maps.
{p_end}

{phang}
{opt subtitle(string)} overrides the automatic subtitle. The default subtitle
reports the kernel, VCE type, and confidence level from stored results.
{p_end}

{phang}
{opt scheme(schemename)} applies a Stata graph scheme.
{p_end}

{phang}
{opt name(name)} assigns a window name to the graph. Combined with
{cmd:replace}, an existing graph window of the same name is replaced.
{p_end}

{phang}
{opt replace} permits overwriting of an existing saved file or named graph
window.
{p_end}

{phang}
{opt nodraw} creates the graph object without rendering it on screen.
{p_end}

{title:Remarks}

{p 4 4 2}
{cmd:rd2d_plot} requires active {cmd:e()} results from {cmd:rd2d} or
{cmd:rd2d_dist}. If neither command has been run in the current estimation
frame, the command exits with error 198.

{p 4 4 2}
The distinction between {cmd:interval(ci)} and {cmd:interval(cb)} reflects a
fundamental difference in coverage guarantee. Pointwise intervals cover each
individual evaluation point at the stated level. Uniform bands cover all
evaluation points simultaneously at the stated level, providing a stronger
simultaneous-inference statement appropriate for assessing boundary
heterogeneity across multiple locations.

{p 4 4 2}
When the preceding estimation used {cmd:side(left)} or {cmd:side(right)}, the
stored intervals are one-sided. The plot adapts: for a left-sided interval the
lower bound is open (negative infinity), and the displayed spike runs from the
estimate to the finite upper endpoint. Analogous logic applies to right-sided
intervals.

{p 4 4 2}
For {cmd:type(heat)}, the command separates evaluation points by sign and
renders them in distinct colors. Marker area encodes the magnitude of the
treatment effect. This representation is informative when multiple boundary
points span a spatially extended boundary segment and the analyst seeks to
identify regions of large or small effects.

{title:Examples}

{phang}Basic effect plot after rd2d:{p_end}
{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0  0 1  0 -1) h(1.5) cbands}{p_end}
{phang2}{cmd:. rd2d_plot}{p_end}

{phang}Display conventional estimates without intervals:{p_end}
{phang2}{cmd:. rd2d_plot, estimate(p) interval(none)}{p_end}

{phang}Uniform confidence bands only:{p_end}
{phang2}{cmd:. rd2d_plot, interval(cb)}{p_end}

{phang}Both pointwise CI and uniform CB overlaid:{p_end}
{phang2}{cmd:. rd2d_plot, interval(both)}{p_end}

{phang}Subset of evaluation points:{p_end}
{phang2}{cmd:. rd2d_plot, subset(1 3) type(effect)}{p_end}

{phang}Heat map of treatment effects at boundary coordinates:{p_end}
{phang2}{cmd:. rd2d_plot, type(heat)}{p_end}

{phang}Save graph to file:{p_end}
{phang2}{cmd:. rd2d_plot, saving(rd2d_effects.gph) replace}{p_end}

{title:References}

{p 4 8 2}
Cattaneo, M. D., R. Titiunik, and R. Yu. 2025. Boundary regression
discontinuity design: location and distance approaches with bivariate score.
The boundary-coordinate target and uniform inference framework implemented by
{cmd:rd2d} and visualized by {cmd:rd2d_plot} are defined in that source.

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rdbw2d}, {help rdbw2d_dist}
