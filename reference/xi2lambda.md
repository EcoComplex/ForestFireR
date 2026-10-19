# Convert a spread or regrowth probability to a rate constant

The model parametrizes fire spread and post-fire regrowth with two
dimensionless probabilities: \\\xi\\, the probability that a burning
neighbour ignites a given vegetation type before the fire itself burns
out, and \\\eta\\, the probability that an unoccupied post-fire site
regrows into a given species before it is recolonised some other way.
Both convert to a rate constant through the same relation, obtained by
setting \\\xi = \lambda^{sp}/(\lambda^{sp}+\lambda^{F\emptyset})\\ (and
analogously for \\\eta\\) and solving for the rate: \$\$\lambda =
\frac{\xi}{1-\xi}\\\lambda\_{ref}\$\$

## Usage

``` r
xi2lambda(xi, lambda_ref)

eta2lambda(xi, lambda_ref)
```

## Arguments

- xi:

  Probability in `[0, 1)` (named `xi` for fire spread, `eta` for
  regrowth – both use this same formula).

- lambda_ref:

  The reference rate the probability is defined against (`L_30`, i.e.
  \\\lambda^{F\emptyset}\\, for fire spread; `L_01`, i.e.
  \\\lambda\_{\emptyset nat}\\, for regrowth – see Details).

## Value

The corresponding rate constant (same units as `lambda_ref`).

## Details

Fire spread: `xi2lambda(xi_nat, L_30)` gives `Lsp_13` (native catches
fire from a burning neighbour), `xi2lambda(xi_inv, L_30)` gives `Lsp_23`
(invader catches fire from a burning neighbour). \\\xi = 0.5\\ is the
point where spread and burnout are equally fast: below it a fire usually
dies where it started, above it a fire usually propagates.

Regrowth: both species' post-fire regrowth probabilities convert the
same way, using `L_01` (not `L_02`) as the reference rate for both:
`eta2lambda(eta_nat, L_01)` gives `Lrg_01` (native long-range regrowth),
`eta2lambda(eta_inv, L_01)` gives `Lrg_02` (invader long-range
regrowth). Setting `L_01 = 0` therefore switches off both regrowth
channels as well as ordinary native colonisation.

The baseline parameter set leaves the native species without a regrowth
advantage (`eta_nat = 0`); `eta_nat` exists so the two species can be
parametrized symmetrically.

## See also

[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md),
[ForestFireR-model](https://EcoComplex.github.io/ForestFireR/reference/ForestFireR-model.md)

## Examples

``` r
# A fire-spread probability of 0.6 against a burnout rate of 1e6/year
xi2lambda(0.6, 1e6)
#> [1] 1500000

# xi = 0.5 returns exactly the reference rate (spread as fast as burnout)
xi2lambda(0.5, 1e6)
#> [1] 1e+06

# Regrowth uses L_01 as its reference rate
eta2lambda(0.6, 0.03)
#> [1] 0.045

# The conversion is what the simulate functions apply internally
r <- default_rates()
xi2lambda(r$xi_inv, r$L_30)
#> [1] 1500000
```
