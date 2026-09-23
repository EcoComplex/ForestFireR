# Generate an initial spatial landscape

Places the invader species on a grid of natives according to the initial
invader density and the spatial heterogeneity parameter `p` (called
\\p\\ in the paper, `q22` internally), following the pair-correlation
landscape-generation algorithm described in the paper's supplementary
material. `p = 0` disperses the invader uniformly at random; `p = 1`
aggregates it into a single compact cluster; intermediate values
interpolate between the two.

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

  Spatial heterogeneity parameter, in `[0, 2]` (0 = random, 1 =
  neutral/no preference in the pairwise statistic, 2 = maximally
  segregated/clustered – see the paper's supplementary material for the
  full derivation).

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

## Value

An `L`x`L` integer matrix: `1` = native vegetation, `2` = invader.

## Examples

``` r
if (FALSE) { # \dontrun{
g <- generate_landscape(L = 100, density2 = 0.1, p = 1, seed = 1)
image(g)
} # }
```
