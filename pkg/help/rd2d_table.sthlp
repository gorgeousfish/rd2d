{smcl}
{* *! version 1.1.0 24jun2026}{...}
{title:rd2d_table}

{p 4 4 2}
{cmd:rd2d_table} repacks estimation results from {cmd:rd2d} or {cmd:rd2d_dist}
into a standard Stata eclass object ({cmd:e(b)}, {cmd:e(V)}), enabling direct
use with {cmd:esttab}, {cmd:estout}, {cmd:estimates table}, and other
post-estimation tabulation tools.

{title:Syntax}

{p 8 8 2}
{cmd:rd2d_table}
[{cmd:using} {it:filename}]
[{cmd:,}
{cmd:replace}
{cmd:format(}{it:string}{cmd:)}
{cmd:estimate(}{it:q}|{it:p}{cmd:)}
{cmd:subset(}{it:numlist}{cmd:)}
{cmd:tex}
{cmd:csv}
{cmd:nostar}
{cmd:level(}{it:#}{cmd:)}
{cmd:title(}{it:string}{cmd:)}]

{title:Description}

{p 4 4 2}
{cmd:rd2d_table} transforms the multi-point estimation results stored by
{cmd:rd2d} or {cmd:rd2d_dist} into a single-equation eclass representation.
Each evaluation point becomes one coefficient in {cmd:e(b)}, with its
corresponding variance (or full covariance, when {cmd:e(cov_q)} or
{cmd:e(cov_p)} is available) stored in {cmd:e(V)}.

{p 4 4 2}
After running {cmd:rd2d_table}, the posted estimation is compatible with all
standard Stata post-estimation commands that read {cmd:e(b)} and {cmd:e(V)},
including {cmd:estimates store}, {cmd:estimates table}, {cmd:esttab},
{cmd:estout}, and {cmd:lincom}.

{p 4 4 2}
If a {cmd:using} filename is specified and {cmd:esttab} (from the {cmd:estout}
package) is installed, the command directly exports the table. Without
{cmd:using}, the command displays results via {cmd:ereturn display}.

{p 4 4 2}
If the current {cmd:e(cmd)} is neither {cmd:rd2d} nor {cmd:rd2d_dist}, the
command exits with error 301.

{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:estimate(}{it:q}|{it:p}{cmd:)}}selects bias-corrected ({cmd:q},
default) or conventional ({cmd:p}) estimates from {cmd:e(results)}{p_end}
{synopt:{cmd:subset(}{it:numlist}{cmd:)}}restricts posting to the specified
evaluation-point indices (integers 1 to {cmd:e(neval)}){p_end}
{synopt:{cmd:format(}{it:string}{cmd:)}}numeric display format for esttab
export; default {cmd:%9.4f}{p_end}
{synopt:{cmd:replace}}allows overwriting an existing output file{p_end}
{synopt:{cmd:tex}}forces LaTeX output format when using a filename{p_end}
{synopt:{cmd:csv}}forces CSV output format when using a filename{p_end}
{synopt:{cmd:nostar}}suppresses significance stars in esttab output{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level for display; default inherits
from the prior estimation{p_end}
{synopt:{cmd:title(}{it:string}{cmd:)}}table title passed to esttab{p_end}
{synoptline}

{phang}
{opt estimate(q)} (the default) posts the bias-corrected robust estimates and
standard errors from columns 5-6 of {cmd:e(results)}. When {cmd:e(cov_q)} is
available (from a prior {cmd:cbands} estimation), the full covariance matrix is
used in {cmd:e(V)}. Otherwise, {cmd:e(V)} is diagonal with squared standard
errors.
{p_end}

{phang}
{opt estimate(p)} posts the conventional polynomial estimates from columns 3-4
of {cmd:e(results)}. Uses {cmd:e(cov_p)} for the full covariance when available.
{p_end}

{phang}
{opt subset(numlist)} filters the posted coefficients to the listed indices.
Indices must be positive integers not exceeding {cmd:e(neval)}.
{p_end}

{phang}
{opt using filename} triggers export via {cmd:esttab}. If {cmd:esttab} is not
installed, the command exits with error 199 and suggests installation via
{cmd:ssc install estout}. The file format is inferred from the extension
({cmd:.tex} or {cmd:.csv}) unless overridden by {cmd:tex} or {cmd:csv}.
{p_end}

{title:Coefficient naming}

{p 4 4 2}
Each posted coefficient is named {cmd:tau_}{it:b1}{cmd:_}{it:b2} where
{it:b1} and {it:b2} are the boundary coordinates with dots replaced by
{cmd:d} and minus signs by {cmd:n}. For example, boundary point (0, 27.5)
becomes {cmd:tau_0_27d5}. When boundary coordinates are absent
({cmd:rd2d_dist} without explicit boundary points), names are
{cmd:tau_1}, {cmd:tau_2}, etc.

{title:Examples}

{pstd}Basic usage: post and display{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.2) p(1) q(2)}{p_end}
{phang2}{cmd:. rd2d_table}{p_end}

{pstd}Export to LaTeX:{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.2) p(1) q(2)}{p_end}
{phang2}{cmd:. rd2d_table using results.tex, replace tex}{p_end}

{pstd}Export a subset to CSV:{p_end}

{phang2}{cmd:. rd2d_table using results.csv, subset(1 3 5) csv replace}{p_end}

{pstd}Use with estimates store for multi-model tables:{p_end}

{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.2) p(1) q(2)}{p_end}
{phang2}{cmd:. rd2d_table}{p_end}
{phang2}{cmd:. estimates store model1}{p_end}
{phang2}{cmd:. rd2d y x1 x2 treat, at(0 0 0.5 0.5) h(1.5) p(1) q(2)}{p_end}
{phang2}{cmd:. rd2d_table}{p_end}
{phang2}{cmd:. estimates store model2}{p_end}
{phang2}{cmd:. esttab model1 model2, se}{p_end}

{pstd}Conventional estimates:{p_end}

{phang2}{cmd:. rd2d_table, estimate(p)}{p_end}

{pstd}After rd2d_dist:{p_end}

{phang2}{cmd:. rd2d_dist y dist treat, at(0 0 0.3 0.7) h(0.8)}{p_end}
{phang2}{cmd:. rd2d_table}{p_end}

{title:Stored results}

{p 4 4 2}
{cmd:rd2d_table} posts as an {cmd:eclass} command and stores the following:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}total number of observations{p_end}
{synopt:{cmd:e(N0)}}control observations (if available){p_end}
{synopt:{cmd:e(N1)}}treated observations (if available){p_end}
{synopt:{cmd:e(neval)}}number of evaluation points posted{p_end}
{synopt:{cmd:e(p)}}polynomial order{p_end}
{synopt:{cmd:e(q)}}bias-correction order{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rd2d_table}{p_end}
{synopt:{cmd:e(cmd_source)}}original estimation command ({cmd:rd2d} or
{cmd:rd2d_dist}){p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(vce)}}variance-covariance method{p_end}
{synopt:{cmd:e(estimate_order)}}{cmd:q} or {cmd:p}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}1 x J coefficient vector (point estimates){p_end}
{synopt:{cmd:e(V)}}J x J variance-covariance matrix{p_end}
{synopt:{cmd:e(results)}}J x K submatrix of original {cmd:e(results)}{p_end}

{title:Remarks}

{p 4 4 2}
{cmd:rd2d_table} overwrites the current {cmd:e()} object. To preserve the
original estimation, use {cmd:estimates store} before calling {cmd:rd2d_table},
or re-run the estimation afterward.

{p 4 4 2}
When {cmd:e(cov_q)} is available (from estimation with {cmd:cbands}), the
posted {cmd:e(V)} contains the full joint covariance. This enables correct
inference for linear combinations across evaluation points via {cmd:lincom}.

{title:References}

{p 4 8 2}
Cattaneo, Titiunik, and Yu (2025). The estimation objects reformatted by
{cmd:rd2d_table} are defined by the boundary RD method and software sources.
See {help rd2d} and {help rd2d_dist} for the original estimation commands.

{title:Also see}

{p 4 4 2}
{help rd2d}, {help rd2d_dist}, {help rd2d_summary}, {help rd2d_aggregate},
{help esttab} (if installed)
