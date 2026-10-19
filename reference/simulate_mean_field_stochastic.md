# Simulate the well-mixed stochastic (non-spatial) model

Simulates the continuous-time Markov chain of the well-mixed
(mass-action, no spatial structure) system with `N` individuals, the
stochastic counterpart of
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md).
Pairwise reactions have propensity \\\lambda N_i N_j / N\\ and
first-order ones \\\lambda N_i\\, so the process has the mean-field
equations as its large-`N` limit. Run against
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
at \\N = L^2\\ it separates the effect of spatial structure from that of
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
  Default `eta_nat = 0` leaves the native species without a post-fire
  regrowth advantage; default `eta_inv = 0.6` gives the invader a
  moderate one.

- L_01, L_02, L_12, L_21, L_30, Lig_13, Lig_23:

  Base rate constants; see
  [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
  for the baseline values. Override individual entries to explore
  parameter space, e.g. `L_21 = 0.02`.

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

- record_dt:

  Trajectory sampling interval. `-1` = final state only, `0` = every
  event, `>0` = sample interval. Default `T/200`.

## Value

A data frame with columns `time`, `n0`..`n5`.

## Details

Unlike the deterministic limit, where a residual burning density keeps
the fire feedback running at all times, here fire is a discrete event:
ignition happens at rate \\\lambda^{ig} N\_{inv}\\ and an outbreak
usually burns out at once. The chain can also reach a true absorbing
state before `T` – all-native with no unoccupied sites and no ignition
source, say – in which case the trajectory stops there. Both are reasons
a finite-`N` run can end somewhere the ODE does not.

## See also

[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md),
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)

## Examples

``` r
ws <- simulate_mean_field_stochastic(T = 200, N = 10000, n1_0 = 0.9,
                                     n2_0 = 0.05, n4_0 = 0.05, seed = 1)
plot(ws$time, ws$n1, type = "l", col = "forestgreen", ylim = c(0, 1),
     xlab = "time (years)", ylab = "density")
lines(ws$time, ws$n2, col = "firebrick")

# Demographic noise alone, against the deterministic limit
mf <- simulate_mean_field(T = 200, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05)
lines(mf$time, mf$n1, col = "forestgreen", lty = 2)
lines(mf$time, mf$n2, col = "firebrick", lty = 2)
```
