# Set the underlying C RNG seed directly

The simulation engine uses the C `rand()`/`srand()` generator (as the
original model code did), which is independent of R's own RNG. Every
`simulate_*()`/[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
call below takes its own `seed` argument (preferred, for reproducible
individual runs); call `set_seed()` directly only if you need to seed
once and then make several *unseeded* calls in a reproducible sequence.

## Usage

``` r
set_seed(seed)
```

## Arguments

- seed:

  Integer seed.
