# Named grid-state codes

The integer codes used in the `L`x`L` grids returned/accepted by
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md),
[`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)
and
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md):
`1` = native vegetation (\\V^{nat}\\), `2` = invader vegetation
(\\V^{inv}\\), `3` = active fire (\\F\\), `4` = unoccupied,
never-yet-burned space (\\\emptyset\\), `5` = unoccupied, post-fire
space (\\\emptyset^F\\, the only one of the two unoccupied states with a
long-range regrowth channel, via `eta_nat`/`eta_inv`; the baseline
parametrization sets `eta_nat = 0`, leaving only the invader with a
regrowth advantage there). State `0` never appears – it's an unused
index in the underlying engine's state numbering.

## Usage

``` r
FF_STATE
```

## See also

[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)

## Examples

``` r
FF_STATE
#>         native        invader           fire          empty empty_postfire 
#>              1              2              3              4              5 
g <- generate_landscape(L = 40, density2 = 0.1, p = 1, seed = 1)
sum(g == FF_STATE["invader"]) / length(g)
#> [1] 0.1
```
