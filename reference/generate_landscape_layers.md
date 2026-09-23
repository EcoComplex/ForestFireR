# Build a landscape from several placement layers

Generalizes
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)
to construct initial conditions it can't express on its own: more than
one pattern placed independently over a shared background, or a small
pattern placed over an already-uniform domain. The whole grid is first
filled with `fill_state`; then each layer in `layers`, in order, grows
its `pattern` state to `density * L*L` cells (a fraction of the WHOLE
grid, same convention as
[`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)'s
`density2`) with heterogeneity `p` (same meaning as there), drawing
candidate cells only from whatever `background` state is still left
after the layers before it. If a layer's target can't be reached because
its background ran out (e.g. the density budget across layers doesn't
leave enough), the engine prints a one-time notice and stops that layer
early rather than hanging – check the printed cell counts against what
you asked for if you see this.

## Usage

``` r
generate_landscape_layers(L = 100, fill_state, layers, seed = NULL)
```

## Arguments

- L:

  Grid side length.

- fill_state:

  The state
  ([FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md))
  the whole grid starts as, before any layer is placed.

- layers:

  A list of layers, each itself a list with elements `background`,
  `pattern` (states from
  [FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)),
  `density` (fraction of `L*L`) and `p` (heterogeneity, as in
  [`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)).

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

## Value

An `L`x`L` integer matrix using
[FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)'s
codes.

## Details

Two scenarios this is built for: (1) two vegetation species placed
independently over a shared background (e.g. a small amount of native
and invader vegetation scattered or clustered over an otherwise empty,
post-fire domain), and (2) a small pattern placed over an
already-uniform domain (e.g. a few active-fire cells ignited directly
over an all-invader domain), for studying post-fire regrowth or fire
spread in relative isolation from the rest of the dynamics – see the
examples below.

## Examples

``` r
if (FALSE) { # \dontrun{
# two species placed independently over a post-fire-empty domain
g <- generate_landscape_layers(L = 100, fill_state = FF_STATE["empty_postfire"],
  layers = list(
    list(background = FF_STATE["empty_postfire"], pattern = FF_STATE["native"],
         density = 0.05, p = 0.8),
    list(background = FF_STATE["empty_postfire"], pattern = FF_STATE["invader"],
         density = 0.05, p = 0.8)
  ), seed = 1)

# a small ignition patch over an all-invader domain
g2 <- generate_landscape_layers(L = 100, fill_state = FF_STATE["invader"],
  layers = list(
    list(background = FF_STATE["invader"], pattern = FF_STATE["fire"],
         density = 0.01, p = 0.8)
  ), seed = 1)
} # }
```
