# Generate an initial spatial landscape

Places the invader on a grid of natives at a given initial density and
with a given degree of clustering `p`. The invader is grown cell by cell
from a random seed site while the neighbour pair statistics are steered
towards the targets set by `p`: writing \\\rho\\ for the invader density
reached so far, the targeted probabilities for a randomly chosen pair of
adjacent sites are \\\rho p\\ for two invader sites, \\2\rho(1-p)\\ for
a mixed pair and \\1 + \rho(p-2)\\ for two native sites. `p = 0` thus
asks for no invader-invader pairs at all (maximally dispersed) and
`p = 1` for every neighbour of an invader site to be an invader (one
compact cluster).

## Usage

``` r
generate_landscape(L = 100, density2 = 0.1, p = 1, seed = NULL)
```

## Arguments

- L:

  Grid side length (grid has `L*L` sites).

- density2:

  Initial invader density, in `(0, 1)`.

- p:

  Clustering of the invader, in `[0, 1]`: `0` scatters it over the grid,
  `1` aggregates it into a single compact cluster (see Details).

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

## Value

An `L`x`L` integer matrix: `1` = native vegetation, `2` = invader.

## Details

The realised clustering increases monotonically with `p` but lags the
target, because the growth process stops once the requested density is
reached. Measured on a 100x100 grid at `density2 = 0.15`, the fraction
of an invader site's four neighbours that are also invader is about 0.00
at `p = 0`, 0.02 at `p = 0.5`, 0.50 at `p = 0.75` and 0.95 at `p = 1`;
values above `1` give the same configuration as `p = 1`. Compare
landscapes at several `p` rather than reading `p` itself as the
clustering.

## See also

[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md),
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)

## Examples

``` r
g <- generate_landscape(L = 60, density2 = 0.1, p = 1, seed = 1)
table(g)                       # 1 = native, 2 = invader
#> g
#>    1    2 
#> 3240  360 

# the same density, scattered instead of clustered
g0 <- generate_landscape(L = 60, density2 = 0.1, p = 0, seed = 1)

op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
image(g0, col = c("forestgreen", "firebrick"), axes = FALSE, main = "p = 0")
image(g,  col = c("forestgreen", "firebrick"), axes = FALSE, main = "p = 1")

par(op)
```
