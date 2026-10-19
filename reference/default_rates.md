# Baseline rate constants

Returns the model's baseline parameter set as a named list, ready to be
passed to
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
or
[`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)
via [`do.call()`](https://rdrr.io/r/base/do.call.html), or used as a
starting point whose individual entries can be overridden. These are the
values every simulate function uses by default. Time is in years.

## Usage

``` r
default_rates()
```

## Value

A named list with the vegetation turnover rates (`L_01`, `L_02`, `L_12`,
`L_21`), the vegetation death rates (`L_10`, `L_20`, both 0), the fire
burnout rate (`L_30`, i.e. \\\lambda^{F\emptyset}\\), the two species'
spontaneous ignition rates (`Lig_13`, `Lig_23`), and the dimensionless
`xi_nat`, `xi_inv`, `eta_nat`, `eta_inv` (see
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).

## Details

Under this parametrization the native species replaces the invader twice
as fast as the reverse, so without fire the invader cannot spread; its
advantages are spontaneous ignition, a higher fire-spread probability
and post-fire regrowth. See
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)
for the mechanism and
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)
for the dimensionless probabilities.

## See also

[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md),
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)

## Examples

``` r
str(default_rates())
#> List of 13
#>  $ L_01   : num 0.03
#>  $ L_02   : num 0.03
#>  $ L_10   : num 0
#>  $ L_20   : num 0
#>  $ L_12   : num 0.005
#>  $ L_21   : num 0.01
#>  $ L_30   : num 1e+06
#>  $ Lig_13 : num 0
#>  $ Lig_23 : num 1e-04
#>  $ xi_nat : num 0.5
#>  $ xi_inv : num 0.6
#>  $ eta_nat: num 0
#>  $ eta_inv: num 0.6

# Run with the baseline set, raising only the invader's regrowth advantage
pars <- default_rates()
pars$eta_inv <- 0.8
pars[c("L_10", "L_20")] <- NULL          # not arguments of the simulate functions
mf <- do.call(simulate_mean_field, c(list(T = 200), pars))
tail(mf, 1)
#>       time        n0         n1        n2           n3           n4        n5
#> 201 200.01 0.3267421 0.02260069 0.6506568 4.601405e-08 0.0003568795 0.3263852
```
