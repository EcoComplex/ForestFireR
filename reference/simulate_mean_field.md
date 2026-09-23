# Simulate the deterministic mean-field model

Integrates the mean-field ODE system describing the well-mixed limit of
the model (see the paper's mean-field derivation for the full
equations), and returns the full trajectory as a data frame. Internally
uses a semi-implicit (IMEX) step for the fire compartment so that it
stays numerically stable even under the paper's default `L_30 = 1e6` (an
explicit fixed-step RK4 on this term alone diverges for any practical
step size) – see the package's implementation notes
(`vignette("forest-fire-model")` or `NEWS.md`) for details. The
underlying vegetation/fire equations are unchanged from the original
model.

## Usage

``` r
simulate_mean_field(
  T,
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
  record_dt = NULL
)
```

## Arguments

- T:

  Simulation time horizon.

- n1_0, n2_0, n3_0, n4_0, n5_0:

  Initial densities for native, invader, fire,
  empty-after-vegetation-death and empty-after-fire respectively (should
  sum to 1).

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

- record_dt:

  Trajectory sampling interval (must be `> 0`; default `T/200`).

## Value

A data frame with columns `time`, `n0`..`n5`.

## Examples

``` r
if (FALSE) { # \dontrun{
mf <- simulate_mean_field(T = 200, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05)
plot(mf$time, mf$n1, type = "l", ylim = c(0,1))
lines(mf$time, mf$n2, col = "red")
} # }
```
