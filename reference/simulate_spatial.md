# Simulate the spatially-explicit stochastic model

Runs the model's Gillespie stochastic simulation algorithm on an `L`x`L`
periodic (or bounded) grid, starting from an invader landscape generated
the same way as
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md).
This is the direct R-friendly interface to the model's continuous-time
Markov chain (master equation) spatial dynamics, as described in the
paper.

## Usage

``` r
simulate_spatial(
  T,
  L = 100,
  density2 = 0.1,
  p = 1,
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
  capture_fire_snapshots = FALSE,
  fire_snapshot_min_gap = 0.001
)
```

## Arguments

- T:

  Simulation time horizon (years).

- L:

  Grid side length (grid has `L*L` sites). Runtime scales roughly with
  `L^2 * T`; for interactive exploration keep `L` around 50-100, reserve
  `L = 200` (paper resolution) for longer background runs (a single
  `L = 200`, `T = 100` run with default rates takes on the order of a
  couple of minutes in a typical environment).

- density2:

  Initial invader density.

- p:

  Spatial heterogeneity of the initial landscape (see
  [`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)).

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

## Value

A list with the final `time_sim` and densities `n0`..`n5` (`n0` = total
empty, `n1` = native, `n2` = invader, `n3` = fire, `n4` =
empty-after-vegetation-death, `n5` = empty-after-fire), plus
`trajectory` (a data frame, if `record_dt >= 0`), `initial_grid` /
`final_grid` (if `record_grid = TRUE`), and `fire_snapshot_times` (a
numeric vector) / `fire_snapshot_grids` (a list of `L`x`L` integer
matrices, same length and order) if `capture_fire_snapshots = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
r <- simulate_spatial(T = 20, L = 60, density2 = 0.05, p = 1,
                       xi_inv = 0.6, eta_inv = 0.6, record_dt = 0.5)
plot(r$trajectory$time, r$trajectory$n1, type = "l")
} # }
```
