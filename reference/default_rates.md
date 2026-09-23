# Table 1 default rate constants

Returns the baseline transition rates used throughout the paper (Table
1), as a named list ready to be passed to
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
[`simulate_mean_field()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field.md)
or
[`simulate_mean_field_stochastic()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_mean_field_stochastic.md)
via [`do.call()`](https://rdrr.io/r/base/do.call.html), or used as a
starting point that individual entries can be overridden from.

## Usage

``` r
default_rates()
```

## Value

A named list with the vegetation turnover rates (`L_01`, `L_02`, `L_12`,
`L_21`), the fire persistence rate (`L_30`, i.e.
\\\lambda^{F\emptyset}\\), the invader's spontaneous ignition rate
(`Lig_23`), the native's spontaneous ignition rate (`Lig_13`), and
default dimensionless `xi_nat`, `xi_inv`, `eta_nat`, `eta_inv` (see
[`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).
