# Convert a spread/regrowth probability to a rate constant

The paper parametrizes fire spread and post-fire regrowth with two
dimensionless probabilities, \\\xi\\ (probability that a burning
neighbour ignites a given vegetation type before the fire itself burns
out) and \\\eta\\ (probability that an empty-after-fire site regrows
into the invader before it is recolonised some other way). Both convert
to a rate constant through the same relation, derived from setting \\\xi
= \lambda^{sp}/(\lambda^{sp}+\lambda^{F\emptyset})\\ (and analogously
for \\\eta\\) and solving for the rate: \$\$\lambda =
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
  \\\lambda\_{\emptyset V1}\\, for regrowth – see Details).

## Value

The corresponding rate constant (same units as `lambda_ref`).

## Details

Fire spread: `xi2lambda(xi_nat, L_30)` gives `Lsp_13` (native catches
fire from a burning neighbour), `xi2lambda(xi_inv, L_30)` gives `Lsp_23`
(invader catches fire from a burning neighbour).

Regrowth: both species' post-fire regrowth probabilities convert the
same way, using `L_01` (not `L_02`) as the reference rate for both:
`eta2lambda(eta_nat, L_01)` gives `Lrg_01` (native long-range regrowth),
`eta2lambda(eta_inv, L_01)` gives `Lrg_02` (invader long-range
regrowth). The paper's own parametrization sets `eta_nat = 0` (native
has no post-fire regrowth advantage); `eta_nat` exists so the model can
be used symmetrically for other systems.
