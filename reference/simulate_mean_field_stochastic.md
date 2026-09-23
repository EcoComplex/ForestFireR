# Simulate the well-mixed stochastic (non-spatial) model

Runs the well-mixed Gillespie simulation (mass-action, no spatial
structure) at a given population size `N`, the stochastic counterpart to
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md).
Useful for isolating the effect of spatial structure from the effect of
demographic stochasticity.

## Usage

``` r
simulate_mean_field_stochastic(
  T,
  N = 10000,
  n1_0 = 0.9,
  n2_0 = 0.05,
  n3_0 = 0,
  n4_0 = 0.05,
  n5_0 = 0,
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
  seed = NULL,
  record_dt = NULL
)
```

## Arguments

- T:

  Simulation time horizon.

- N:

  Well-mixed population size (analogous to `L*L` in the spatial model,
  but with no spatial structure).

- n1_0, n2_0, n3_0, n4_0, n5_0:

  Initial densities (should sum to 1).

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

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

- record_dt:

  Trajectory sampling interval. `-1` = final state only, `0` = every
  event, `>0` = sample interval. Default `T/200`.

## Value

A data frame with columns `time`, `n0`..`n5`. Note: the CTMC can reach a
true absorbing state (e.g. all-native with no empty sites and no
residual ignition source) before `T`, in which case the trajectory
simply stops there.
