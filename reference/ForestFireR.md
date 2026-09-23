# ForestFireR: Spatial and Mean-Field Forest-Fire Invasion Model

R/Rcpp interface to the spatial stochastic and mean-field forest-fire
invasion model from "Dynamics of invasion in forest-fire models".

## Details

Start with
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
to build an initial spatial configuration,
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
to run the full spatially-explicit Gillespie simulation,
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
for the deterministic well-mixed limit, and
[`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)
for its stochastic (well-mixed, no spatial structure) counterpart.
[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
and
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
help set up rate parameters using a compact dimensionless notation
(baseline values, and the \\\xi\\/\\\eta\\ dimensionless spread/regrowth
probabilities).

See the `examples/` directory shipped with the package source for worked
scripts illustrating the model's main behaviors (an invasion trajectory,
a landscape-heterogeneity sweep, a fire-spread/regrowth phase diagram).
If you use this package, please cite the paper it implements – see the
README's Citation section.
