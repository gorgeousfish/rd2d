# Changelog

## [1.2.0] - 2026-07-01

### Changed (BREAKING)
- `cbands` is now enabled by default for both `rd2d` and `rd2d_dist`, matching the R package default
- Users who do not want uniform confidence bands must now specify `nocbands`
- The old `cbands` option is retained for backward compatibility (specifying it is now a no-op)

### Added
- `nocbands` option to suppress uniform confidence bands computation

## [1.1.0] - 2026-06-24

### Added
- Shared masspoints warning utility (`_rd2d_masspoints_warn.ado`) for consistent diagnostics across all commands
- Example do-files: `rd2d_cluster_example.do` (cluster VCE), `rd2d_cbands_example.do` (confidence bands), `rd2d_imse_example.do` (IMSE bandwidth selection)
- Extended Monte Carlo validation test (`verify_phase5_monte_carlo_extended.do`) with 3 DGPs x 2 sample sizes
- Cross-command consistency test (`verify_phase6_cross_command_consistency.do`)
- Mata function unit tests (`verify_phase6_mata_unit_tests.do`)
- `stata.toc` file for `net install` support
- Mata function interface documentation for all internal functions
- Data file codebook (`data/README.md`)

### Changed
- Version bump to 1.1.0 across all package files
- Error messages now display actual invalid values and kernel aliases across all 4 commands
- Table display uses centralized layout parameters (line_width, hline_rule, layout mode)
- Parameter validation messages standardized with valid options and received values
- `Nh0`/`Nh1` stored as integers in results matrices
- Help files expanded with `deriv()`, `tangvec()`, `cluster()` documentation
- Stored results sections completed with all return values documented

### Fixed
- `masspoints(off)` control flow error in `rd2d.ado`
- Missing `bwsource` local when user provides `h()` directly
- `tangvec()` warning format inconsistency between `rd2d.ado` and `rdbw2d.ado`

### Added (post-estimation commands)

- `rd2d_plot`: post-estimation plotting with effect plots and heat maps.
- `rd2d_summary`: structured display of estimation results with bandwidth diagnostics.
- `rd2d_aggregate`: aggregation of boundary treatment effects (WBATE, AATE, LBATE per Theorems 4-6).
- `rd2d_table_example.do`: example script for LaTeX/CSV table construction from e(results).
- `verify_phase6_boundary_conditions.do`: regression tests for one-sided CI signs, PSD repair, cluster edge cases, extreme bandwidths, kink, and masspoints diagnostics.
- Help files for all three post-estimation commands (rd2d_plot.sthlp, rd2d_summary.sthlp, rd2d_aggregate.sthlp).

## [1.0.0] - 2026-05-10

Initial Stata package v1.0.0 release surface for boundary regression
discontinuity designs.

### Added

- `rd2d` for location-based bivariate local polynomial estimation and
  inference.
- `rdbw2d` for location-based bandwidth selection.
- `rd2d_dist` for signed-distance local polynomial estimation and inference.
- `rdbw2d_dist` for signed-distance bandwidth selection.
- `rdbw2d` now accepts `cluster()` as an alias for `c()` when supplying the
  selector cluster identifier.
- `rd2d_dist` now accepts side-specific manual `h(h0 h1)` bandwidths and
  per-distance-column manual bandwidth pairs, matching the R public surface.
- Robust bias-corrected pointwise intervals for smooth-boundary inference.
- Optional Gaussian-simulation uniform confidence bands through `cbands`.
- HC0-HC3 variance estimators, HC0/HC1 cluster paths, mass-point diagnostics,
  generalized-inverse fallback diagnostics, and confidence-band PSD repair
  diagnostics.
- English README, Stata help files, reproducible location and distance
  examples, and net-install package metadata.

### Verification

- Formula-level WLS fixtures check the Stata local polynomial helpers without
  using R as truth.
- R parity and generated-example E2E tests cover bandwidths, estimates,
  standard errors, confidence intervals, effective sample sizes, and
  manual-bandwidth branches.
- Monte Carlo fast tests report bias, SD, RMSE, pointwise coverage, uniform-band
  coverage, interval length, band length, MCSE, and effective sample sizes.
- The optional paper-style Monte Carlo audit is now included in the tolerance
  registry with an explicit 5-MCSE coverage window for its default `reps=12`
  configuration.
- Package install smoke stages `rd2d.pkg` and `stata.toc`, installs into an
  isolated PLUS directory twice, checks all public commands and help files, and
  reruns both packaged examples after `net get`.
- Final release check passed on 2026-05-07 with
  `make check-release-final STATA=/usr/local/bin/stata-mp` in 252.16 seconds,
  covering 17 source anchors, 4 formula fixture cases, 4336 R-baseline numeric
  comparisons, documented R-reference regressions, 4 public commands, 22
  requirements, 12 package files, generated-example E2E and Monte Carlo
  validation, full Stata smoke, install smoke, and failure count 0.
- Post-release reference checks cover local R documentation repairs for
  `rd2d.dist()`/`rdbw2d.dist()` `kink` defaults, the location `rd2d()`
  `Se.p`/`Se.q` result-label documentation, zero-distance mass-point
  side-count behavior, the `rd2d_dgp.R` linear-DGP coefficient used by the
  paper table, and the `rdbw2d()` derivative `msetwo`/`imsetwo` denominator
  convention.
- A post-release source/live audit now verifies all 7 documented kernel aliases
  (`uni`, `unif`, `tri`, `triag`, `epa`, `epan`, and `gau`) across all four
  public commands.
- A post-release source/live audit aligned the location `scaleregul()` default
  with the paper default `3` for `rd2d` and `rdbw2d`, while preserving explicit
  `scaleregul(1)` for cross-language reference comparisons.
- A post-release release-surface audit documented the existing `rdbw2d`
  `deriv()` and `tangvec()` selector targets in help, README, package mirrors,
  and release checks.
- A post-release release-surface audit documented and now verifies the live
  automatic selector defaults: `bwselect(mserd)` for all public commands and
  `method(dpi)` for location commands.
- A post-release release-surface audit documented and now verifies the live
  kernel-shape defaults: `kernel(triangular)` for all public commands,
  `ktype(prod)` for location commands, and `kink(off)` for distance commands.
- The same default-surface audit documented adjacent omitted-option defaults:
  `p(1)`, `masspoints(check)`, location `scalebiascrct(1)`, estimation
  `rbc(on)`, and estimation `side(two)`.
- A post-release source/live audit repaired malformed manual `h()` validation
  for `rd2d` and `rd2d_dist`: nonnumeric, nonpositive, and missing bandwidth
  tokens now return controlled `r(198)` errors instead of falling through to
  Stata expression parsing; release checks keep the missing `h(.)` probes.
- A post-release control-surface audit added focused numeric-option checks
  for location commands: malformed `at()`, `deriv()`, and `tangvec()` tokens
  plus zero or missing `tangvec()` direction vectors now have permanent
  `r(198)` regression coverage for `rdbw2d` and `rd2d`.
- A post-release release-surface audit closed missing scale-constant
  validation: all public commands now reject missing `scaleregul(.)`, and
  location commands reject missing `scalebiascrct(.)`, with controlled
  `r(198)` failures; release checks keep those numeric-option probes in place.
- A distance control-surface audit expanded the numeric-option checks to
  keep missing `cqt(.)` validation locked for `rdbw2d_dist` and `rd2d_dist`.
- A location control-surface audit expanded the numeric-option checks from
  32 to 34 probes so missing `at(. 0)` coordinates are locked for `rdbw2d`
  and `rd2d`.
- A bandwidth-control audit expanded the numeric-option checks from 34 to 38
  probes so missing `bwcheck(.)` inputs are locked for all four public
  commands.
- An estimation-control boundary audit expanded the numeric-option checks from
  38 to 46 probes so out-of-range `level(0)` / `level(100)` and nonpositive
  `repp(0)` / `repp(-1)` inputs are locked for both estimation commands.
- A q-order validation audit expanded the numeric-option checks from 46 to
  50 probes so negative and missing estimation `q()` inputs are rejected with controlled
  `r(198)` for both `rd2d` and `rd2d_dist` instead of being treated as omitted
  defaults.
- A follow-on q-order validation audit expanded the numeric-option checks
  from 50 to 52 probes so explicit `q()` values below `p()` remain locked as
  controlled `r(198)` failures for both estimation commands.
- A signed-distance input audit expanded the numeric-option checks from 52
  to 54 probes so `rdbw2d_dist` and `rd2d_dist` reject inconsistent
  multi-column distance side signs with controlled `r(198)` before bandwidth,
  estimation, or covariance calculations begin.
- A polynomial-order validation audit expanded the numeric-option checks from
  54 to 66 probes so missing `p(.)`, fractional `p(1.5)`, and negative
  `p(-1)` inputs remain locked as controlled `r(198)` failures across
  `rd2d`, `rdbw2d`, `rd2d_dist`, and `rdbw2d_dist`.
- A location treatment-indicator audit expanded the numeric-option checks
  from 66 to 68 probes so `rdbw2d` and `rd2d` reject non-binary treatment
  indicators with controlled `r(198)` instead of silently treating all
  nonzero values as treated.
- A location target-shape audit expanded the numeric-option checks from 68 to
  69 probes so `rd2d` rejects `tangvec()` with multiple `at()` points using a
  controlled `r(198)` failure while preserving the single-point `tangvec()`
  path.
- A cluster-input audit repaired the Stata public `cluster()` / `c()` contract:
  all four public commands now accept string cluster identifiers by mapping
  them to temporary numeric groups inside the marked estimation sample, while
  preserving the user-supplied cluster variable name in estimation returns.
  The Stata complete-case cluster check now compares numeric-versus-string
  equality on generated package-reference location and signed-distance paths at
  tolerance `1e-12`.
- A post-release test-index audit aligned `test/README.md` with the expanded
  validation checks and keeps the README aligned with missing `h(.)`, missing
  `scaleregul(.)`, missing `scalebiascrct(.)`, and missing distance `cqt(.)`
  coverage.
- A post-release cluster-HC1 audit aligned inference and bandwidth-selector
  constants with the single cluster finite-sample multiplier, without also
  applying the non-cluster HC1 residual multiplier. Focused verification now
  reports `rd2d` `N=16`, clusters=8, `Se.q=0.020127`; `rd2d_dist` `N=12`,
  clusters=6, `Se.q=0.189710`; tolerance `1e-11`.
- Generated-example and packaged-example graphs now have a release check that
  confirms the exported PNG files are nonblank `1600x960` displays with visible
  data regions, readable guide text, boundary or zero-distance reference lines,
  fitted curves, and side-specific observations in the plot region.  The check
  keeps title and legend marks from being counted as data evidence and confirms
  the generated-example location display includes the selected boundary
  evaluation markers.
- A post-release output-surface audit now locks large public display contracts:
  100-point location `cbands` output with `at99`/`at100` rows and `100 x 100`
  `e(V)`, 100-column distance `cbands` output with full stored row/column
  names, 100-column `rdbw2d_dist` long-label suffixes, signed-distance
  bandwidth values around `312.6`, large location coordinates printed on one
  row, narrow `linesize(63)` selector/estimation tables, and compact
  `linesize(50)` tables with shortened `CI.lo`/`CI.hi` and `CB.lo`/`CB.hi`
  endpoint labels that keep bandwidth columns and confidence-band endpoints
  readable while preserving stored matrix names and columns. One-sided open
  endpoints now print as `-inf`/`inf` while stored matrices keep the
  `-c(maxdouble)`/`c(maxdouble)` sentinels.
- A post-release stored-result semantics audit now documents and verifies
  selector `r(mseconsts)` `Nh` columns as bandwidth-constant in-band counts,
  automatic estimation `e(bws)` `Nh0`/`Nh1` columns as final p-fit in-band
  counts, location selector absence of stored `bwmin`/`bwmax` columns, and
  distance selector `bwmin0`/`bwmin1`/`bwmax0`/`bwmax1` raw-support versus
  unique-support clamp semantics under `bwcheck()`.

### Known R Reference Issues And Repairs

- The R `rd2d_cb()` right-sided one-sided band path can report a negative
  simulation critical value; the Stata package uses positive one-sided
  critical values.
- The R `rdbw2d.dist()` masspoint diagnostics classify zero-distance unique
  support on the control side, even though the same function treats
  nonnegative distances as treated observations.
- The local R `rd2d()` roxygen/Rd result documentation previously described
  live `Se.p`/`Se.q` standard-error columns as `Var.p`/`Var.q`; the local
  documentation has been repaired.
- The local R `rd2d.dist()` result table and roxygen/Rd documentation
  previously labeled live standard-error columns as `Var.p`/`Var.q`; the
  local R package now exposes them as `Se.p`/`Se.q`.
- The local R source/data bundle does not exactly reproduce the paper's
  displayed `rd2d.dist()` example table.
- The R `rdbw2d.dist()` implementation default for `bwcheck` differs from the
  paper and R documentation.
- The R `rd2d.dist()` q-order distance standard error differs from the
  paper-consistent raw-basis WLS variance oracle used by the Stata package.
- The R `rdbw2d()` `stdvars=TRUE` path returns `b1`/`b2` in standardized
  coordinates while bandwidths are returned on the original scale.
- The R `rdbw2d()` runtime default for location `scaleregul` is `1`, while the
  paper, R documentation, and R `rd2d()` wrapper default to `3`; Stata uses the
  paper-aligned location default `3`.
- The bundled R `rd2d_dgp.R` linear-DGP control-side `x2` coefficient differs
  from the paper table; Stata Monte Carlo validation uses the paper coefficient.
- The local R roxygen and Rd documentation for `cbands` and
  `rd2d.dist()`/`rdbw2d.dist()` `kink` defaults were repaired to match the
  current function signatures and paper defaults; the location `rd2d()` and
  distance `rd2d.dist()` `Se.p`/`Se.q` result-label documentation were also
  repaired. The regression hooks keep those repairs aligned.

### Limitations

- The Stata package does not call R at runtime.
- Distance commands do not accept separate boundary-coordinate metadata, so
  distance result matrices include `b1` and `b2` placeholder columns whose
  entries are missing by design.
- With `cluster()`, HC2/HC3 requests are reset to HC1 because the implemented
  cluster path is HC0/HC1.
- Distance `kink(on)` uses the undersmoothing convention and should not be
  described as robust bias correction.
