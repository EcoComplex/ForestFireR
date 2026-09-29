# Set the underlying C RNG seed directly

The simulation engine uses its own random number generator, independent
of R's: by default the 64-bit Mersenne Twister (`std::mt19937_64`),
whose random sequence for a given seed is the same on every platform
(runs themselves can still differ across platforms, because the model
computes in `long double`, whose precision is platform-dependent). Set
`options(ForestFireR.rng = "legacy")` to use the C library
`rand()`/`srand()` generator of versions up to 0.7.x instead (only to
reproduce earlier runs; its sequence differs between platforms). Every
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
