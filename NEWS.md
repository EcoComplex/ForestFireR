# NEWS

## ForestFireR 0.6.0

### New features

* Added `eta_nat`, the native counterpart to `eta_inv`: a dimensionless
  post-fire regrowth probability for the native species, converted to the
  `Lrg_01` rate constant the same way `eta_inv` converts to `Lrg_02` (see
  `?xi2lambda`). The underlying C++ engine already supported `Lrg_01`
  symmetrically across all three sub-engines (spatial, mean-field ODE,
  well-mixed stochastic) -- only the R-level parameter was missing.
  Default `eta_nat = 0` matches the paper's own parametrization exactly,
  so this is backward compatible: no existing call changes behavior.
* Added to all four simulate functions (`simulate_spatial()`,
  `simulate_mean_field()`, `simulate_mean_field_stochastic()`,
  `simulate_spatial_from_grid()`) and to `simulate_spatial_movie()`;
  added `eta_nat = 0` to `default_rates()`'s returned list.
* Updated roxygen docs for `xi2lambda()`/`eta2lambda()`, `default_rates()`,
  `simulate_spatial()` (and the functions that `@inheritParams` from it),
  and `FF_STATE` to describe both species' regrowth channels
  symmetrically, and regenerated `man/*.Rd`.
* Added `test-native-regrowth.R`: isolates each species' post-fire
  regrowth channel in turn (fire, replacement and the other species'
  `eta` all off) starting from an all-empty-postfire grid, and checks the
  a priori known outcome -- only the species whose `eta` is nonzero can
  ever appear. Also confirmed `eta_nat = eta_inv = 0` leaves the grid
  untouched. Extended `test-xi2lambda.R`'s `default_rates()` check to
  cover `xi_nat`, `xi_inv`, `eta_nat`, `eta_inv`.
* Note for anyone isolating regrowth dynamics like the new test does:
  `eta2lambda()`/`xi2lambda()` convert both `eta_nat` and `eta_inv` using
  `L_01` (not `L_02`) as the reference rate (see `?xi2lambda`'s Details),
  so `L_01 = 0` zeroes `Lrg_01` too, not just ordinary native
  colonization -- keep `L_01` nonzero when you want either regrowth
  channel active.

### Documentation site

* Added a `pkgdown` site (`_pkgdown.yml`) instead of a vignette, per
  <https://r-pkgs.org/vignettes.html>'s recommendation. Verified it builds
  cleanly (`pkgdown::build_site()`, all 12 reference pages, no errors)
  before committing. Deployment is automated via
  `.github/workflows/pkgdown.yaml` (the standard `r-lib/actions` template):
  every push to `main` rebuilds the site and publishes it to the
  `gh-pages` branch -- the site itself is not committed to `main` (see
  `.gitignore`). One manual step remains: in the GitHub repo's Settings ->
  Pages, set the source to the `gh-pages` branch (it's created
  automatically the first time the workflow runs after this is pushed).
  Once enabled, the site is served at
  <https://EcoComplex.github.io/ForestFireR/>.

## ForestFireR 0.5.0

### Testing

* Added a `testthat` (edition 3) unit test suite (`tests/testthat/`, the
  package's first), 8 files and 130+ expectations. All check exact, a
  priori known outcomes rather than statistical properties:
  - `test-null-dynamics.R`: with every rate at zero, the grid comes out
    byte-for-byte identical to how it went in and `time_sim = 0`.
  - `test-single-species-colonization.R`: with only one species present
    and no fire, colonization deterministically fills every empty site
    with that species (`n1 = 1` or `n2 = 1`, exactly).
  - `test-ignition-disabled.R`: with `Lig_13 = Lig_23 = 0`, fire and
    post-fire empty ground never appear, even with fire-spread rates set
    high.
  - `test-isolated-fire-burn.R`: with fire spread and new ignition both
    off, an isolated burning cell burns out to exactly one
    empty-post-fire site and nothing else in the grid changes.
  - `test-reproducibility.R`: a fixed seed reproduces the exact same run
    (spatial and well-mixed stochastic engines), and different seeds do
    not coincide.
  - `test-mass-conservation.R`: `n1+n2+n3+n4+n5 = 1` at every recorded
    point of the trajectory, for all three engines (spatial, mean-field
    ODE, well-mixed stochastic).
  - `test-no-na-regression.R`: a regression guard against the NA-output
    class of bug previously found in Fig. 6 (see the repo's
    `fig6_lambda_na_bug` notes) -- asserts finite, non-negative,
    conserved output under Table 1's own default rates.
  - `test-xi2lambda.R`: exact closed-form checks of the `xi`/`eta` ->
    rate conversion and of `default_rates()`'s Table 1 values.
* Added `testthat (>= 3.0.0)` to `Suggests` and
  `Config/testthat/edition: 3` to `DESCRIPTION`.

### Investigated, not a bug

* Confirmed that `simulate_*()`'s returned `time_sim` can exceed the
  requested `T`, sometimes by a non-trivial margin. This is expected
  Gillespie behavior, not a bug: the loop must apply the event that
  crosses the requested horizon in order to report a valid state, and
  that event's own waiting time can be large whenever the total
  propensity is low at that point in the run. Tests that check timing
  assert `time_sim >= T` (barring extinction/absorption), never
  closeness to `T`.

## ForestFireR 0.4.01

### Bug fixes

* Fixed `InitialConditionNonHomogeneous()` (used by `generate_landscape_layers()`
  for multi-species initial landscapes, e.g. Fig. 8 and Fig. 9 of the
  supplementary material). Pair classification during the clustering
  growth process relied on numeric-mean thresholds
  (`threshold1 = (state1+state1)/2`, etc.) that only classify pairs
  correctly when `state1 < state2`. For native/invader layers laid over an
  `empty_postfire` background (`background = 5 > pattern states 1, 2`),
  this assumption did not hold, corrupting cluster counts for every
  clustering parameter `p` and silently producing landscapes with the
  wrong spatial structure (e.g. very different `p` values could produce
  nearly identical patterns at `t = 0`).
* Replaced the threshold-based classification with direct state-equality
  comparisons (`S[i][j] == state1`, etc.) in both the pair-counting loop
  and the two move-evaluation blocks (`p1`, `p2`) of
  `InitialConditionNonHomogeneous()`. Verified by cluster-count regression
  tests against the single-layer reference implementation, and by a full
  paper-scale re-run of Fig. 9 showing the expected monotonic ordering by
  `p` and the expected cluster shapes.
* Only initial-landscape generation for **two or more layers** was
  affected (Fig. 8 and Fig. 9's reduced-dynamics scenarios); single-layer
  landscapes (Fig. 10, Fig. 11) were not affected.

### Investigated, not a bug

* At `eta_inv = 0` (fully symmetric parameters, no fire, no replacement),
  `rho_nat` for `p = 1` sits well below the ~0.5 expected by symmetry at
  `T = 300` (Fig. 8/9's simulation horizon), while lower-`p` curves are
  already close to 0.5 -- visible in the original published Fig. 8 too.
  Investigated and confirmed this is NOT an asymmetry between native and
  invader: at `T = 300`, `p = 1` gives `n1 ~= n2 ~= 0.30` with ~40% of the
  grid still unconverted empty-post-fire background (`n5 ~= 0.40`);
  letting the same runs continue to `T ~= 1000` brings both `n1` and `n2`
  to ~0.50 with `n5 = 0`. Root cause: colonization of empty-post-fire
  cells is driven by local perimeter (neighbour count), and a fully
  clustered (`p = 1`) landscape has far less perimeter per occupied cell
  than a fragmented one, so it consumes the remaining background much
  more slowly. `T = 300` is simply not long enough for the `p = 1` curve
  to reach the steady state that lower-`p` curves already reached --
  native and invader are affected equally (`n1` and `n2` track each other
  throughout). Not a code defect; a modeling/figure-design consideration
  (whether Fig. 8/9's fixed `T` should be treated as a snapshot of
  clustering-dependent transient dynamics, or whether the `p = 1` curve
  needs a longer horizon to compare like-for-like against the others).

## ForestFireR 0.4.0

* Initial versioned release used for the manuscript figures (mean-field
  and spatial Gillespie simulation engine, R wrapper functions using the
  paper's own notation, figure-reproduction scripts).
