# ForestFireR: Spatial and Mean-Field Forest-Fire Invasion Model

Simulation of a two-species vegetation system in which fire mediates the
competition between a native species and an invader, on an explicit
lattice or in its well-mixed limit. The engine is written in C++ and
exposed to R through Rcpp.

## Details

See
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)
for the definition of the model: states, transitions, parameters and the
mean-field equations.

To run it:
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
builds an initial spatial configuration,
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
runs the lattice Gillespie simulation,
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
integrates the deterministic well-mixed limit and
[`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)
simulates the well-mixed stochastic process at a finite population size.
[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
gives the baseline parameter set and
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
converts the dimensionless spread/regrowth probabilities into rate
constants.
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)
and
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)
cover initial conditions that
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
cannot express, and
[`simulate_spatial_movie()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_movie.md)
renders a run as a video.

Worked scripts (an invasion trajectory, a landscape-clustering
comparison, a fire-spread/regrowth phase diagram) are installed with the
package:

    list.files(system.file("examples", package = "ForestFireR"))

## See also

Useful links:

- <https://github.com/EcoComplex/ForestFireR>

- <https://EcoComplex.github.io/ForestFireR/>

- Report bugs at <https://github.com/EcoComplex/ForestFireR/issues>

## Author

**Maintainer**: Leonardo A. Saravia <lasaravia@untdf.edu.ar>

Authors:

- Leonardo A. Saravia <lasaravia@untdf.edu.ar>

- Rebeca de la Fuente
