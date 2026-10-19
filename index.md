# ForestFireR

R package simulating a two-species vegetation system in which fire
mediates the competition between a native species and an invader. The
same continuous-time transition rules are simulated three ways: on an
explicit lattice with the Gillespie algorithm, as the deterministic
well-mixed (mean-field) limit, and as the well-mixed stochastic process
at a finite population size. The engine is written in C++ and exposed to
R through Rcpp, with an interface built around dimensionless parameters
(`xi`, `eta`, `p`) for fire spread, post-fire regrowth and the
clustering of the initial landscape.

The model is documented inside the package:
[`help("ForestFireR-model")`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)
defines the five states, every transition and its rate, the
dimensionless parametrization, the baseline parameter set and the
fire-mediated invasion mechanism it encodes, the mean-field equations,
and what differs between the three engines.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("EcoComplex/ForestFireR")
```

## Quick start

``` r

library(ForestFireR)

# Deterministic well-mixed limit
mf <- simulate_mean_field(T = 200, n1_0 = 0.9, n2_0 = 0.1)
plot(mf$time, mf$n1, type = "l", col = "forestgreen", ylim = c(0, 1),
     xlab = "time (years)", ylab = "density")
lines(mf$time, mf$n2, col = "firebrick")

# Lattice simulation, 10 per cent invader clustered in one patch
sp <- simulate_spatial(T = 200, L = 100, density2 = 0.1, p = 1,
                       seed = 1, record_dt = 1)
c(native = sp$n1, invader = sp$n2)

# Baseline parameters, and the probability-to-rate conversion
default_rates()
xi2lambda(0.6, 1e6)
```

## What the package provides

| function | purpose |
|----|----|
| [`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md) | initial lattice: one species placed over another with density and clustering `p` |
| [`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md) | initial lattice built from a sequence of placements over an arbitrary fill state |
| [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md) | lattice Gillespie simulation |
| [`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md) | same, starting from a supplied grid |
| [`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md) | deterministic well-mixed limit |
| [`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md) | well-mixed stochastic process at population size `N` |
| [`simulate_spatial_movie()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_movie.md) | render a run as a video or animated GIF |
| [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md), [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md), [`eta2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md) | baseline parameters and the dimensionless conversion |
| `FF_STATE` | the grid’s state codes |

Recording options: `record_dt` for trajectories, `record_grid` for
initial and final grids, `capture_fire_snapshots` for grids captured
while a fire is actually burning (a single outbreak resolves in about an
hour of simulated time against decades of vegetation dynamics, so a
fixed-time snapshot essentially never lands on one). Every simulate
function takes its own `seed`; the engine’s generator is
`std::mt19937_64`, independent of R’s.

Two things to know when building initial conditions with
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md):

- **Native-extinction early stop.** The Gillespie loop treats native
  density below `0.0001` as extinction and ends the run there. That is
  right whenever the landscape starts with native vegetation, and wrong
  for a grid built without any (a fire-over-invader domain, say), where
  `n1` is `0` from the first step. Pass `check_extinction = FALSE` to
  [`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)
  for those.
- **Layer density budget.** A layer’s target (`density * L*L` cells) can
  exceed what is left in its background state. The layer then stops
  early rather than searching for an unreachable target, and prints how
  many cells it actually placed.

## Worked scripts

Installed with the package:

``` r

list.files(system.file("examples", package = "ForestFireR"))
#> "01_trajectories_meanfield_vs_spatial.R"  # the three engines side by side
#> "02_landscape_heterogeneity.R"            # clustering of the initial landscape
#> "03_phase_diagram_xi_eta_sweep.R"         # (xi_inv, eta_inv) phase diagram
```

Every function’s help page also carries runnable examples, including the
reduced-dynamics recipes (post-fire regrowth alone, fire spread alone).

## Behaviour of the implementation

Starting from 10 per cent invader and the baseline parameters
(`xi_nat = 0.5`, `xi_inv = eta_inv = 0.6`, `eta_nat = 0`) over 200
years:

- **Mean field** (`n1_0 = 0.9`, `n2_0 = 0.1`): native declines from 0.90
  to 0.079 while the invader grows from 0.10 to 0.613.
- **Lattice** (`L = 100`, `density2 = 0.1`, `p = 1`, `seed = 1`): native
  ends at 0.050 and the invader at 0.675; the same run with a scattered
  invader (`p = 0`) ends at 0.057 native and 0.336 invader, the rest of
  the grid still unoccupied post-fire ground.
- **Well mixed, stochastic** (`N = 10000`): the invasion does not take
  off (native 0.971, invader 0.029 at `T = 200`). Fire here is a
  discrete, rare event, whereas the deterministic limit carries a
  residual burning density of order `1e-11` at all times, which keeps
  the feedback running; without fire the native species out-replaces the
  invader by construction.

Sweeping `xi_inv` and `eta_inv`
(`inst/examples/03_phase_diagram_xi_eta_sweep.R`) maps the transition to
native extinction over the plane of the two fire-mediated advantages.

Mass conservation (`n1+n2+n3+n4+n5 = 1`) holds to floating-point
precision in all three engines: the largest deviation over a recorded
`L = 100`, `T = 200` trajectory is `1.1e-16`. The `testthat` suite
checks exact, a priori known outcomes (null dynamics, single-species
colonization, isolated fire burnout, ignition switched off,
reproducibility, mass conservation, the `xi`/`eta` conversion, and the
random stream itself).

## Performance

Measured on an 8-core laptop with the default sum-tree engine
(`engine = "tree"`), median of three runs:

| run | wall time |
|----|----|
| `simulate_spatial(T = 200, L = 100)` | 0.05 s |
| `simulate_spatial(T = 200, L = 200)` | 0.24 s |
| `simulate_spatial(T = 200, L = 200, xi_inv = 0.9, eta_inv = 0.9)` | 0.37 s |
| `simulate_spatial(T = 1000, L = 200)` | 0.55 s |
| `simulate_mean_field(T = 300)` | 0.004 s |
| `simulate_mean_field_stochastic(T = 200, N = 10000)` | 0.001 s |

Runtime grows with the number of Gillespie events, so with grid size,
time horizon and with parameters that keep fires burning.
`record_dt = -1` (final state only) is the cheapest option for parameter
sweeps. `engine = "linear"` reproduces runs made with ForestFireR \<
0.7.0 and is orders of magnitude slower on large grids.

## Citation

``` r

citation("ForestFireR")
```

The package implements the model analysed in Saravia, L.A. and de la
Fuente, R., *Dynamics of invasion in forest-fire models* (manuscript;
analysis code at <https://github.com/EcoComplex/InvasionForestFires>).
Cite the package for the implementation and the article for the model.

## Layout

- `R/model.R` – the model reference page
  ([`help("ForestFireR-model")`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)).
- `R/rates.R` –
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md),
  [`eta2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md),
  [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md).
- `R/simulate.R` – landscape generation, the three simulation engines,
  `FF_STATE`.
- `R/movie.R` –
  [`simulate_spatial_movie()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_movie.md).
- `src/forest_fire_model.cpp` – the C++ engine (lattice Gillespie
  dynamics, well-mixed Gillespie dynamics, mean-field integration) and
  the Rcpp export layer.
- `inst/examples/` – worked scripts.
- `tests/testthat/` – unit tests.
