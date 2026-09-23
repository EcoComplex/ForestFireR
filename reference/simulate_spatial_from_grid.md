# Simulate the spatially-explicit model from a supplied landscape

Like
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
but starts from an initial grid you supply (typically from
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md))
instead of building one internally from `density2`/`p`. Use this for any
initial condition
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
can't express on its own, typically alongside some of the rate constants
zeroed out to isolate a subset of the dynamics. Two useful recipes: to
study post-fire regrowth/colonization in isolation, disable fire
entirely (`xi_nat = 0, xi_inv = 0, Lig_23 = 0`) and vegetation
replacement (`L_12 = 0, L_21 = 0`); to study fire spread in isolation
(e.g. starting from a small ignition patch over an all-invader domain),
disable colonization/regrowth (`L_01 = 0, L_02 = 0, eta_inv = 0`) and
native-side reactions (irrelevant when there's no native vegetation in
the initial condition), keeping only fire spread/burnout.

## Usage

``` r
simulate_spatial_from_grid(
  T,
  initial_grid,
  xi_nat = 0.5,
  xi_inv = 0.6,
  eta_nat = 0,
  eta_inv = 0.6,
  L_01 = 0.03,
  L_02 = 0.03,
  L_12 = 0.005,
  L_21 = 0.01,
  L_30 = 1e+06,
  Lig_13 = 0,
  Lig_23 = 1e-04,
  periodic = TRUE,
  seed = NULL,
  record_dt = -1,
  record_grid = FALSE,
  check_extinction = TRUE,
  capture_fire_snapshots = FALSE,
  fire_snapshot_min_gap = 0.001
)
```

## Arguments

- T:

  Simulation time horizon.

- initial_grid:

  An `L`x`L` integer matrix using
  [FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)'s
  codes (e.g. from
  [`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)).

- xi_nat, xi_inv:

  Fire-spread probabilities for native/invader, in `[0, 1)` (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).
  `xi = 0.5` is the spread/no-spread threshold (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)'s
  Details); the defaults here (`xi_nat = 0.5`, `xi_inv = 0.6`) put the
  native species right at that threshold and the invader moderately
  above it.

- eta_nat, eta_inv:

  Native/invader post-fire regrowth probabilities (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).
  Default `eta_nat = 0` matches the paper's own parametrization (no
  native post-fire regrowth advantage); default `eta_inv = 0.6`
  represents a moderate post-fire regrowth advantage for the invader.

- L_01, L_02, L_12, L_21, L_30, Lig_13, Lig_23:

  Base rate constants; see
  [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
  for Table 1 defaults. Override individual entries to explore parameter
  space, e.g. `L_21 = 0.02`.

- periodic:

  Use periodic boundary conditions (paper default `TRUE`).

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

- record_dt:

  Trajectory sampling interval, in the same time units as `T`. `-1`
  (default) records only the final state (fastest, use this for
  parameter sweeps). `0` records every single event (dense, memory
  heavy, only for small/short runs). A positive value samples the
  trajectory roughly every `record_dt` time units (use this whenever you
  need the full time series rather than just the endpoint, e.g.
  `record_dt = T/200` for a smooth-looking curve).

- record_grid:

  If `TRUE`, also return the initial and final grid snapshots (`L`x`L`
  integer matrices) – use this whenever you need to inspect or plot the
  spatial pattern, not just the aggregate densities.

- check_extinction:

  The engine's usual early stop – treat the native species as extinct
  and end the run as soon as its density drops below `0.0001` – assumes
  `initial_grid` actually contains native vegetation to begin with. Set
  `FALSE` for an `initial_grid` that has none by design (e.g. a
  fire-over-invader domain with no native vegetation), or that early
  stop triggers on the very first step. Leave at the default `TRUE`
  whenever native vegetation is present (matches
  [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
  which always uses the check).

- capture_fire_snapshots:

  If `TRUE`, also return every grid snapshot captured while the fire
  compartment was actively non-empty (`n3 > 0`), throttled to at most
  one capture per `fire_snapshot_min_gap` time units. Use this to find
  an illustrative "fire in progress" grid near a target illustration
  time: `record_grid`'s start/end snapshots essentially never show fire,
  because the default fire-spread/fire-extinction rates (~1e6/year,
  hours-scale) resolve any given outbreak in a sliver of simulated time
  invisible to a snapshot at an arbitrary fixed time – an outbreak has
  to be searched for explicitly, not hoped for. After the run, pick the
  returned `fire_snapshot_times` entry closest to your target time and
  use the matching `fire_snapshot_grids` entry; there is no guarantee
  one exists near any particular time (ignition is itself a rare Poisson
  process – `lambda^ig = 1e-4`/year per invader-occupied site by
  default), so check `length(fire_snapshot_times)` and how close the
  nearest one actually is before trusting it as illustrative of that
  moment.

- fire_snapshot_min_gap:

  Minimum simulated-time gap (years) between two fire snapshots; keeps a
  single outbreak from filling the returned list with near-duplicate
  frames. Default `0.001` (~9 hours) is well under a typical outbreak's
  total duration but well above the time between individual Gillespie
  events during one.
