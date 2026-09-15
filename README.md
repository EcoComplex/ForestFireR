# ForestFireR

R package implementing the forest-fire invasion model from *"Dynamics of
invasion in forest-fire models"*: a spatial stochastic (Gillespie) and
deterministic mean-field model of fire-driven biological invasion in a
two-species, fire-vegetation system. The simulation engine is written in
C++ and exposed to R through Rcpp, with an R-friendly interface built
around a compact dimensionless notation (`xi`, `eta`, `p`) for fire
spread, post-fire regrowth and spatial heterogeneity.

## What this is

The package gives you a documented, R-friendly interface to the model:
`generate_landscape()`/`generate_landscape_layers()` to build initial
spatial configurations, `simulate_spatial()` for the full
spatially-explicit Gillespie simulation, `simulate_mean_field()` for the
deterministic well-mixed limit, and `simulate_mean_field_stochastic()`
for its stochastic well-mixed counterpart, plus recording options
(trajectories, grid snapshots, fire-in-progress snapshots) and
reproducible RNG seeding.


## Features

- Optional trajectory recording (`record_dt` argument) for both the
  spatial and mean-field models, instead of only final densities.
- Optional initial/final grid snapshots (`record_grid = TRUE`) for the
  spatial model.
- A well-mixed *stochastic* counterpart (`simulate_mean_field_stochastic()`)
  to the deterministic mean-field model, for separating the effect of
  spatial structure from demographic stochasticity.
- `generate_landscape()` to build/inspect an initial spatial configuration
  on its own.
- R-facing functions use a compact dimensionless notation (`xi_nat`,
  `xi_inv`, `eta_inv`, `p`) with `xi2lambda()`/`eta2lambda()` handling the
  conversion to the underlying rate constants, and `default_rates()`
  gives Table 1's baseline values.
- `generate_landscape_layers()` + `simulate_spatial_from_grid()`:
  `generate_landscape()` can only place one pattern species over one
  background (its `density2`/`p` arguments, e.g. invader over native).
  Some analyses need initial conditions that don't fit that shape -- two
  species placed *independently* over a shared empty background (e.g. 5%
  native + 5% invader over an all-`empty_postfire` domain), or a small
  pattern placed over an already-uniform domain (e.g. 1% fire over an
  all-invader domain). `generate_landscape_layers()` builds these by
  filling the whole grid with one state and then growing a *sequence* of
  layers onto it, each layer's candidates restricted to whatever
  background the layers before it left behind; `simulate_spatial_from_grid()`
  then runs the usual spatial dynamics starting from that grid instead of
  building one internally. See `FF_STATE` for the grid's state codes
  (`native`/`invader`/`fire`/`empty`/`empty_postfire`) and each function's
  documentation for worked examples. Two things to know before using this:
  - **Native-extinction early stop.** The engine's Gillespie loop treats
    native density dropping below `0.0001` as extinction and stops the run
    there -- correct whenever the initial landscape has native vegetation
    to begin with (as `generate_landscape()`'s `density2`/`p` initial
    conditions always do), but wrong for an `initial_grid` built with no
    native vegetation at all (e.g. a fire-over-invader domain): `n1` is `0`
    from the first step, so the run would otherwise stop after a single
    event. Pass `check_extinction = FALSE` to `simulate_spatial_from_grid()`
    whenever `initial_grid` has no native vegetation by design.
  - **Layer density budget.** A layer's target (`density * L*L` cells) can
    exceed what's actually left in its `background` state (e.g. two
    layers' densities summing to more than what `fill_state` provided).
    When that happens, the layer stops early rather than searching
    indefinitely for an unreachable target (bounded retries, then a
    relaxed fallback selection), and prints a one-time warning naming how
    many cells it actually managed to place.
- `capture_fire_snapshots` (on `simulate_spatial()` and
  `simulate_spatial_from_grid()`): the default fire-spread and
  fire-extinction rates are ~1e6/year (hours-scale) against years-scale
  vegetation dynamics, so any single active-fire episode resolves in a
  sliver of simulated time that a fixed-time grid capture (`record_grid`'s
  start/end snapshots) essentially never lands on. With
  `capture_fire_snapshots = TRUE`, the engine scans for `n3 > 0` after
  every accepted Gillespie event and records that grid (throttled to at
  most one capture per `fire_snapshot_min_gap` simulated-time units, so a
  single outbreak doesn't fill the return value with near-duplicate
  frames). The result gains `fire_snapshot_times` (numeric vector) and
  `fire_snapshot_grids` (list of `L`x`L` integer matrices, same length and
  order) -- pick whichever entry is closest to a target illustration time
  after the run completes. Overhead is small (~1% wall-clock in our
  testing) since `n3` is already recomputed after every event regardless.

## Validation

With the paper's stated baseline parameters (`xi_inv = eta_inv = 0.6`,
`xi_nat = 0.5`, 10% initial invader density, heterogeneous distribution, T
up to 200 years), both `simulate_mean_field()` and `simulate_spatial()`
reproduce the qualitative pattern the paper describes: native density declines
while invader density grows over the 200-year window (mean-field: native
0.90 -> 0.08, invader 0.10 -> 0.61). Sweeping `xi_inv`/`eta_inv` (see
`examples/03_phase_diagram_xi_eta_sweep.R`) reproduces a sharp
continuous-looking transition to native extinction, consistent with the
paper's description of the model's main qualitative result. The
heterogeneity example (`examples/02_landscape_heterogeneity.R`) reproduces
the paper's finding that a segregated/clustered invader (`p = 1`) drives
faster native decline than a randomly dispersed one (`p = 0`).

Mass conservation (`n1+n2+n3+n4+n5 = 1`) was checked to hold to floating
precision for both the spatial and mean-field engines.

## Performance

The spatial model's runtime is dominated by the number of Gillespie events,
which scales with grid size and, once a fire is actively burning, with
`lambda_F0` itself (very short waiting times between fire-related events).
Rough measurements in a modest sandboxed environment:

| L   | T   | wall time |
|-----|-----|-----------|
| 50  | 20  | ~0.1 s    |
| 100 | 10  | ~0.75 s   |
| 200 | 20  | ~19 s     |
| 200 | 100 | ~2+ min (longer if/when fire is actively spreading) |

For interactive exploration, keep `L` around 50-100. Reserve `L = 200`
(paper resolution) for longer, non-interactive/background runs, and use
`record_dt = -1` (final state only) for parameter sweeps where you don't
need the full trajectory.

## Citation

This package implements the model from:

> de la Fuente, R. and Saravia, L.A. "Dynamics of
> invasion in forest-fire models" (manuscript; see the
> https://github.com/EcoComplex/InvasionForestFires repository for the current submitted version).

If you use this package in published work, please cite that paper rather
than (or in addition to) this repository.

## Layout

- `src/forest_fire_model.cpp` -- the C++ simulation engine (Gillespie
  spatial dynamics and mean-field ODE integration) plus the Rcpp export
  layer.
- `R/rates.R` -- `xi2lambda()`, `eta2lambda()`, `default_rates()`.
- `R/simulate.R` -- `generate_landscape()`, `simulate_spatial()`,
  `simulate_mean_field()`, `simulate_mean_field_stochastic()`, `set_seed()`.
- `examples/` -- runnable scripts illustrating trajectories, landscape
  heterogeneity, and a parameter-sweep phase diagram.
