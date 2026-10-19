# Package index

## The model

- [`ForestFireR-model`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)
  [`model`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)
  : The fire-driven invasion model simulated by ForestFireR
- [`ForestFireR`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-package.md)
  [`ForestFireR-package`](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-package.md)
  : ForestFireR: Spatial and Mean-Field Forest-Fire Invasion Model

## Parameters

- [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
  : Baseline rate constants
- [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
  [`eta2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
  : Convert a spread or regrowth probability to a rate constant

## Initial landscapes

- [`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
  : Generate an initial spatial landscape
- [`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)
  : Build a landscape from several placement layers
- [`FF_STATE`](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)
  : Named grid-state codes

## Simulation

- [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)
  : Simulate the spatially-explicit stochastic model
- [`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)
  : Simulate the spatially-explicit model from a supplied landscape
- [`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
  : Simulate the deterministic mean-field model
- [`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)
  : Simulate the well-mixed stochastic (non-spatial) model
- [`simulate_spatial_movie()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_movie.md)
  : Render an animation of a spatial simulation
