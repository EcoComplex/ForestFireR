# The fire-driven invasion model simulated by ForestFireR

Reference definition of the model this package simulates: two vegetation
types, a native species and an invader, competing for space in a system
where fire is the process that lets the invader win. The same transition
rules are simulated in three ways – on an explicit lattice
([`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)),
as the deterministic well-mixed limit
([`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)),
and as the well-mixed stochastic process at a finite population size
([`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)).
This page defines the states, the transitions and every parameter the
package's functions take.

## States

Each site of an `L` x `L` lattice is in one of five states, coded as in
[FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md):

|      |                  |                 |                                   |
|------|------------------|-----------------|-----------------------------------|
| code | name             | symbol          | meaning                           |
| 1    | `native`         | \\V^{nat}\\     | native vegetation                 |
| 2    | `invader`        | \\V^{inv}\\     | invading vegetation               |
| 3    | `fire`           | \\F\\           | burning                           |
| 4    | `empty`          | \\\emptyset\\   | unoccupied, never burned          |
| 5    | `empty_postfire` | \\\emptyset^F\\ | unoccupied, left behind by a fire |

Both unoccupied states can be colonised by a vegetated neighbour, but
only \\\emptyset^F\\ additionally has the long-range post-fire regrowth
channel (see `eta_nat`/`eta_inv` below). State \\\emptyset\\ is produced
only by vegetation death, whose rates the R interface fixes at zero, so
it appears only if the initial landscape contains it.

The simulate functions report the densities `n1`..`n5` (the fraction of
sites in each state, in the order of the table) plus `n0 = n4 + n5`, the
total unoccupied fraction. The five state densities sum to 1 at every
recorded point.

## Transitions

Sites change state in continuous time. On the lattice, \\v_X\\ is the
number of the site's four nearest neighbours (von Neumann neighbourhood,
periodic boundaries by default) that are in state \\X\\:

|  |  |  |
|----|----|----|
| transition | rate | set through |
| \\\emptyset \to V^{nat}\\ | \\\lambda\_{\emptyset nat}\\ v\_{nat}\\ | `L_01` |
| \\\emptyset \to V^{inv}\\ | \\\lambda\_{\emptyset inv}\\ v\_{inv}\\ | `L_02` |
| \\\emptyset^F \to V^{nat}\\ | \\\lambda\_{\emptyset nat}\\ v\_{nat} + \lambda^{rg}\_{nat}\\ | `L_01`, `eta_nat` |
| \\\emptyset^F \to V^{inv}\\ | \\\lambda\_{\emptyset inv}\\ v\_{inv} + \lambda^{rg}\_{inv}\\ | `L_02`, `eta_inv` |
| \\V^{nat} \to V^{inv}\\ | \\\lambda\_{nat\\inv}\\ v\_{inv}\\ | `L_12` |
| \\V^{inv} \to V^{nat}\\ | \\\lambda\_{inv\\nat}\\ v\_{nat}\\ | `L_21` |
| \\V^{nat} \to F\\ | \\\lambda^{ig}\_{nat} + \lambda^{sp}\_{nat}\\ v_F\\ | `Lig_13`, `xi_nat` |
| \\V^{inv} \to F\\ | \\\lambda^{ig}\_{inv} + \lambda^{sp}\_{inv}\\ v_F\\ | `Lig_23`, `xi_inv` |
| \\F \to \emptyset^F\\ | \\\lambda^{F\emptyset}\\ | `L_30` |
| \\V \to \emptyset\\ | \\\lambda\_{nat\\\emptyset}\\, \\\lambda\_{inv\\\emptyset}\\ | fixed at 0 |

Colonisation, competitive replacement and fire spread are local: their
rates are proportional to the matching neighbour count, so an isolated
site is inert to them. Spontaneous ignition, fire burnout and long-range
post-fire regrowth do not depend on the neighbourhood.

## Dimensionless parametrisation

Fire spread and post-fire regrowth are supplied as probabilities rather
than as rate constants, which keeps them comparable across parameter
sets:

- \\\xi_i = \lambda^{sp}\_i / (\lambda^{sp}\_i + \lambda^{F\emptyset})\\
  (`xi_nat`, `xi_inv`) is the probability that a burning neighbour sets
  a site of type \\i\\ alight before that neighbour burns out. \\\xi =
  0.5\\ is the point where spread and burnout are equally fast.

- \\\eta_i = \lambda^{rg}\_i / (\lambda^{rg}\_i + \lambda\_{\emptyset
  nat})\\ (`eta_nat`, `eta_inv`) is the probability that a post-fire
  site is taken by species \\i\\'s long-range regrowth rather than by
  ordinary colonisation. Both species use \\\lambda\_{\emptyset nat}\\
  (`L_01`) as the reference rate.

[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
and
[`eta2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
invert both relations, \\\lambda = \xi/(1-\xi)\\\lambda\_{ref}\\; every
simulate function takes `xi_*`/`eta_*` and converts internally.

The initial landscape adds one more dimensionless parameter, the
clustering `p` of the species placed on it (see
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)).

## Baseline parameters and the invasion mechanism

[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
returns the baseline set, with time in years: \\\lambda\_{\emptyset nat}
= \lambda\_{\emptyset inv} = 0.03\\, \\\lambda\_{nat\\inv} = 0.005\\,
\\\lambda\_{inv\\nat} = 0.01\\, \\\lambda^{F\emptyset} = 10^6\\,
\\\lambda^{ig}\_{nat} = 0\\, \\\lambda^{ig}\_{inv} = 10^{-4}\\,
\\\xi\_{nat} = 0.5\\, \\\xi\_{inv} = 0.6\\, \\\eta\_{nat} = 0\\,
\\\eta\_{inv} = 0.6\\.

These values make the competition asymmetric in the native species'
favour: the native replaces the invader twice as fast as the reverse
(\\\lambda\_{inv\\nat} = 2 \lambda\_{nat\\inv}\\), and the two colonise
unoccupied space equally well, so without fire the invader cannot
spread. All of its advantages act through fire: it is the only species
that ignites spontaneously, it carries fire better (\\\xi\_{inv} \>
\xi\_{nat}\\), and it is the only one that regrows into post-fire sites
through the long-range channel (\\\eta\_{inv} \> \eta\_{nat} = 0\\).
Fire converts native vegetation into \\\emptyset^F\\, which the invader
claims preferentially, which raises the amount of fire-prone vegetation:
invasion and fire reinforce each other. \\\xi\_{inv}\\ and
\\\eta\_{inv}\\ are the strength of that feedback, and are the
parameters most analyses vary.

Fire is also much faster than vegetation: with \\\lambda^{F\emptyset} =
10^6\\/year a site burns for about half an hour, against the decades
over which vegetation turns over. A single outbreak therefore resolves
in a sliver of simulated time (see `capture_fire_snapshots` in
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
for how to observe one).

## Mean-field limit

[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
integrates the well-mixed limit of the same transitions, in which each
neighbour count \\v_X\\ is replaced by the global density \\n_X\\:

\$\$\dot n_1 = \lambda\_{\emptyset nat} n_1 (n_4 + n_5) +
\lambda^{rg}\_{nat} n_5 + \lambda\_{inv\\nat} n_1 n_2 -
\lambda\_{nat\\inv} n_1 n_2 - (\lambda^{ig}\_{nat} + \lambda^{sp}\_{nat}
n_3) n_1\$\$ \$\$\dot n_2 = \lambda\_{\emptyset inv} n_2 (n_4 + n_5) +
\lambda^{rg}\_{inv} n_5 + \lambda\_{nat\\inv} n_1 n_2 -
\lambda\_{inv\\nat} n_1 n_2 - (\lambda^{ig}\_{inv} + \lambda^{sp}\_{inv}
n_3) n_2\$\$ \$\$\dot n_3 = (\lambda^{ig}\_{nat} + \lambda^{sp}\_{nat}
n_3) n_1 + (\lambda^{ig}\_{inv} + \lambda^{sp}\_{inv} n_3) n_2 -
\lambda^{F\emptyset} n_3\$\$ \$\$\dot n_4 = -(\lambda\_{\emptyset nat}
n_1 + \lambda\_{\emptyset inv} n_2) n_4\$\$ \$\$\dot n_5 =
\lambda^{F\emptyset} n_3 - (\lambda\_{\emptyset nat} n_1 +
\lambda\_{\emptyset inv} n_2) n_5 - (\lambda^{rg}\_{nat} +
\lambda^{rg}\_{inv}) n_5\$\$

(the vegetation death terms, zero in the R interface, are omitted). Fire
enters as a continuous density: \\n_3\\ settles at
\\(\lambda^{ig}\_{nat} n_1 + \lambda^{ig}\_{inv}
n_2)/(\lambda^{F\emptyset} - \lambda^{sp}\_{nat} n_1 -
\lambda^{sp}\_{inv} n_2)\\, of order \\10^{-11}\\ under the baseline
parameters – tiny, but always burning, so the feedback runs
continuously. In the lattice and finite-`N` versions fire is instead a
discrete, rare event. That difference, rather than the presence or
absence of space alone, is why the three engines can disagree.

## The three engines

- [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
  and
  [`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md):
  the lattice continuous-time Markov chain, simulated exactly with the
  Gillespie algorithm. Boundaries are periodic unless
  `periodic = FALSE`.

- [`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md):
  the deterministic well-mixed limit above.

- [`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md):
  the continuous-time Markov chain of the same well-mixed system with
  `N` individuals. Pairwise reactions have propensity \\\lambda N_i N_j
  / N\\ and first-order reactions \\\lambda N_i\\, so the process has
  the mean-field equations as its `N` -\> infinity limit, with
  demographic stochasticity and discrete fire at finite `N`. Comparing
  it against
  [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
  at \\N = L^2\\ separates the effect of spatial structure from that of
  demographic noise.

## Initial landscapes

[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
places one species over a background of another with a prescribed
density and clustering `p`;
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)
builds landscapes from a sequence of such placements over an arbitrary
fill state. Both grow the pattern cell by cell while steering the
neighbour pair statistics towards the targets set by `p`, in the spirit
of the pair-correlation landscape generators of Hiebeler (2000).

## Implementation notes

- Gillespie event selection keeps the per-site total rates in a binary
  sum tree (`engine = "tree"`, O(log N) per event). `engine = "linear"`
  is the cumulative scan over all sites, kept for reproducing runs made
  with ForestFireR \< 0.7.0.

- Random numbers come from the engine's own `std::mt19937_64`,
  independent of R's generator, and seeded per call from each function's
  `seed` argument (or from the clock when that is `NULL`). Set
  `options(ForestFireR.rng = "legacy")` to select the `rand()` generator
  of versions up to 0.7.x.

- The mean-field system is stiff: the fire compartment is updated with
  its exact linear solution over each step while the vegetation block is
  integrated with classic RK4, which is stable for any practical step
  size at \\\lambda^{F\emptyset} = 10^6\\.

- A Gillespie step is applied in full, so the reported `time_sim` can
  exceed the requested `T`, by a wide margin when the total rate is low.

## References

Saravia, L.A. and de la Fuente, R. Dynamics of invasion in forest-fire
models (manuscript). Analysis code:
<https://github.com/EcoComplex/InvasionForestFires>.

Gillespie, D.T. (1977) Exact stochastic simulation of coupled chemical
reactions. Journal of Physical Chemistry 81, 2340-2361.

Hiebeler, D. (2000) Populations on fragmented landscapes with spatially
structured heterogeneities: landscape generation and local dispersal.
Ecology 81, 1629-1641.

## See also

[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
for the baseline parameters,
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
for the probability-to-rate conversion,
[FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)
for the state codes.
