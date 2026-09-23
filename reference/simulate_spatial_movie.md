# Render an animation of a spatial simulation

Runs the same Gillespie dynamics as
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md),
but instead of returning only a final state (or the throttled
fire-activity snapshots from `capture_fire_snapshots`), captures the
full `L`x`L` grid at `n_frames` evenly spaced points in SIMULATED TIME
across `[0, T]` and encodes them into a video (default) or an animated
GIF, returning the path to the rendered file.

## Usage

``` r
simulate_spatial_movie(
  T,
  L = 100,
  density2 = 0.1,
  p = 1,
  xi_nat = 0.5,
  xi_inv = 0.6,
  eta_nat = 0,
  eta_inv = 0.6,
  L_01 = 0.03,
  L_02 = 0.03,
  L_12 = 0.005,
  L_21 = 0.01,
  L_30 = 1e+06,
  Lig_13 = 0,
  Lig_23 = 1e-04,
  periodic = TRUE,
  seed = NULL,
  initial_grid = NULL,
  n_frames = 150,
  fps = 15,
  file = NULL,
  colors = NULL,
  check_extinction = TRUE,
  width = 600,
  height = 600,
  verbose = TRUE
)
```

## Arguments

- T:

  Simulation time horizon (years).

- L:

  Grid side length (grid has `L*L` sites). Runtime scales roughly with
  `L^2 * T`; for interactive exploration keep `L` around 50-100, reserve
  `L = 200` (paper resolution) for longer background runs (a single
  `L = 200`, `T = 100` run with default rates takes on the order of a
  couple of minutes in a typical environment).

- density2:

  Initial invader density.

- p:

  Spatial heterogeneity of the initial landscape (see
  [`generate_landscape()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape.md)).

- xi_nat, xi_inv:

  Fire-spread probabilities for native/invader, in `[0, 1)` (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).
  `xi = 0.5` is the spread/no-spread threshold (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)'s
  Details); the defaults here (`xi_nat = 0.5`, `xi_inv = 0.6`) put the
  native species right at that threshold and the invader moderately
  above it.

- eta_nat, eta_inv:

  Native/invader post-fire regrowth probabilities (see
  [`xi2lambda()`](https://EcoComplex.github.io/ForestFireR/reference/xi2lambda.md)).
  Default `eta_nat = 0` matches the paper's own parametrization (no
  native post-fire regrowth advantage); default `eta_inv = 0.6`
  represents a moderate post-fire regrowth advantage for the invader.

- L_01, L_02, L_12, L_21, L_30, Lig_13, Lig_23:

  Base rate constants; see
  [`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)
  for Table 1 defaults. Override individual entries to explore parameter
  space, e.g. `L_21 = 0.02`.

- periodic:

  Use periodic boundary conditions (paper default `TRUE`).

- seed:

  Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.

- initial_grid:

  An `L`x`L` integer matrix using
  [FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)'s
  codes (typically from
  [`generate_landscape_layers()`](https://EcoComplex.github.io/ForestFireR/reference/generate_landscape_layers.md)),
  used instead of building one internally from `density2`/`p` – e.g. to
  seed an active fire directly instead of waiting on spontaneous
  ignition (`Lig_23`). When supplied, `L`, `density2` and `p` are
  ignored. Remember to also set `check_extinction = FALSE` if this
  landscape has no native vegetation by design (see
  [`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)'s
  own docs for this same recipe).

- n_frames:

  Number of checkpoints to capture across `[0, T]` (the grid is also
  captured at `t = 0`, so the rendered animation has `n_frames + 1`
  frames). See Details for what happens when the model's actual dynamics
  are slower or faster than one reaction per `T / n_frames`.

- fps:

  Frames per second in the rendered output.

- file:

  Output path. The extension picks the renderer: `.gif` uses the
  `gifski` package, anything else (`.mp4`, `.mov`, `.webm`, ...) uses
  the `av` package. Defaults to a temporary `.mp4` file if `NULL`.

- colors:

  Named character vector of fill colors, keyed by
  [FF_STATE](https://EcoComplex.github.io/ForestFireR/reference/FF_STATE.md)
  names (`native`, `invader`, `fire`, `empty`, `empty_postfire`). Only
  the names you supply override the package default; you don't need to
  specify all five.

- check_extinction:

  Stop early once native vegetation goes extinct, matching
  [`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)'s
  own always-on behavior (see Details for what happens to the rendered
  video when this triggers). Set `FALSE` to keep simulating – and
  rendering – past extinction.

- width, height:

  Frame size in pixels.

- verbose:

  If `TRUE` (default), print progress as checkpoints are filled.

## Value

Invisibly, the path to the rendered video/GIF file (`file`, resolved to
an absolute path), with attributes `n_engine_calls` (how many
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)
calls it actually took to cover `[0, T]` – compare to `n_frames` per
Details), `stopped_early` (`TRUE` if `check_extinction` triggered before
`T`) and `t_elapsed` (the actual simulated time reached, `== T` unless
`stopped_early`).

## Details

Internally, this chains calls to
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md),
each one asking to advance exactly to the next unfilled checkpoint time
`k * T / n_frames`, starting from the previous call's `final_grid` – the
Gillespie process is Markovian, so this produces a trajectory
statistically identical to one continuous run.

A single Gillespie step cannot be split:
[`simulate_spatial_from_grid()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial_from_grid.md)'s
engine always applies whatever reaction it draws next in full, even when
that reaction's waiting time overshoots the requested horizon (there is
no way to ask it "how much time until the next event, without applying
it" – that would need a change to the underlying C++, which this
function deliberately avoids). So a single call can overshoot its
requested checkpoint by a lot, especially early on or with slow rate
constants (no reaction is possible faster than the model's own smallest
active rate allows) – e.g. requesting `T = 10` sliced into 30 frames
asks for one reaction roughly every `0.33` time units, but if the
landscape's next actual reaction is `12` time units away (very plausible
with, say,
[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)'s
regrowth/replacement rates, which the paper itself only shows producing
visible change over `T` on the order of `100`-`300`), that one call
alone jumps past `12/0.33 ~ 36` requested checkpoints at once. This
function handles that correctly: whenever a call's actual elapsed time
(`time_sim`) reaches or passes one or more pending checkpoints, ALL of
them are filled with that call's resulting grid (correctly showing no
change happened yet), and the next call starts from the next
still-unfilled checkpoint – so the total number of engine calls actually
made is bounded by the number of *real* reactions needed to cover
`[0, T]`, which can be far fewer than `n_frames` for slow dynamics
(efficient) or up to `n_frames` for fast ones (e.g. an active fire),
never more. The practical takeaway: if you see very few distinct frames
(or a `verbose` log with very few calls before reaching `T`), that's the
model genuinely not doing much in `[0, T]` at these rates, not a bug –
pick a larger `T` (see
[`default_rates()`](https://EcoComplex.github.io/ForestFireR/reference/default_rates.md)'s
documentation and the paper for typical horizons) rather than a smaller
`n_frames`.

Each engine call gets its own derived seed (`seed + <call index>`, when
`seed` is supplied) rather than reusing one seed repeatedly, which would
otherwise restart the underlying C RNG from the same state on every
call.

If the native population goes extinct (or the process reaches an
absorbing state, `sumA == 0`) before `T`, matching
[`simulate_spatial()`](https://EcoComplex.github.io/ForestFireR/reference/simulate_spatial.md)'s
own `check_extinction` early stop, all remaining checkpoints are filled
with that final grid, so the rendered video still runs the full
requested length, ending on a frozen final state, instead of silently
being shorter than asked.

## Examples

``` r
if (FALSE) { # \dontrun{
# Needs the 'av' package: install.packages("av")
# T chosen large enough for these default rates to actually do
# something visible -- see Details.
path <- simulate_spatial_movie(T = 150, L = 80, density2 = 0.05, p = 1,
                                n_frames = 150, fps = 15,
                                file = "invasion.mp4", seed = 1)

# Needs the 'gifski' package: install.packages("gifski")
simulate_spatial_movie(T = 150, L = 60, n_frames = 80,
                        file = "invasion.gif", seed = 1)

# Fire spread, seeded directly (default L_30/Lig_23 make fire resolve
# in ~1e-4 years -- invisible to any video; see the fire-spread-only
# reduced-dynamics recipe in simulate_spatial_from_grid()'s docs, and
# this function's Details). L_30 = 5500
# here is grounded in real fire-residence-time estimates for a 30x30 m
# cell (this package's own cell size), not a fitted package value --
# two independent estimates converge on it: (1) understory fire
# rate-of-spread in tropical forest, ~1-3 m/min, takes ~10-30 min to
# cross 30 m; (2) Saravia et al. 2025 (Oikos, doi:10.1111/oik.10764)
# fit a contact-spread fire model (same mechanism as this package's
# fire rule) on a 460 m/pixel Amazon grid with a fixed 1-day-per-site
# fire duration; rescaling that linearly by cell-size ratio
# (460/30 ~ 15x) gives ~1.6 h per 30 m cell. Both give a per-cell fire
# duration on the order of 10 min-2 h, i.e. L_30 ~ 5e3-5e4/year -- far
# from the package's L_30 = 1e6 default (effectively instantaneous)
# but also far from the earlier ad hoc guess of L_30 = 75 (~1 week per
# cell, too slow by ~2 orders of magnitude). T is rescaled down from
# years to ~15 days accordingly, so the fire front is still visibly
# spreading across the grid rather than resolving in a single frame or
# taking years. xi_inv still sets the spread-vs-burnout probability
# exactly as documented in xi2lambda() regardless of this rescaling.
ignite <- generate_landscape_layers(
  L = 80, fill_state = FF_STATE["invader"],
  layers = list(list(background = FF_STATE["invader"],
                      pattern = FF_STATE["fire"], density = 0.01, p = 1)),
  seed = 1
)
simulate_spatial_movie(T = 0.04, initial_grid = ignite,
                        xi_nat = 0, xi_inv = 0.8, eta_inv = 0,
                        L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0,
                        Lig_13 = 0, Lig_23 = 0, L_30 = 5500,
                        check_extinction = FALSE,
                        n_frames = 150, fps = 20,
                        file = "fire_spread.mp4", seed = 1)
} # }
```
