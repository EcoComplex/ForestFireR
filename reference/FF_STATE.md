# Named grid-state codes

The integer codes used in the `L`x`L` grids returned/accepted by
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md),
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)
and
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md):
`1` = native vegetation (paper's \\V^{nat}\\), `2` = invader vegetation
(\\V^{inv}\\), `3` = active fire (\\F\\), `4` = unoccupied,
never-yet-burned space (\\\emptyset\\), `5` = unoccupied, post-fire
space (\\\emptyset^F\\, the only one of the two empty states with a
long-range regrowth channel, via `eta_nat`/`eta_inv`; the paper's own
parametrization sets `eta_nat = 0`, leaving only the invader with a
regrowth advantage there). State `0` never appears – it's an unused
index in the underlying engine's state numbering.

## Usage

``` r
FF_STATE
```

## Format

An object of class `integer` of length 5.
