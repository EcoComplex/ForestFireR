# ForestFireR

R/Rcpp package wrapping the forest-fire invasion model from *"Dynamics of
invasion in forest-fire models"*, built from the original
`forest_fire_model.cpp` engine.

## What this is (and isn't)

This package makes the existing simulation engine callable from R with a
sensible, documented interface, using the paper's own notation. It is
**not** a reimplementation of the model: every reaction, rate and the
Gillespie/RK4 algorithms are the ones already in `forest_fire_model.cpp`.
The changes made here are (a) fixes to real bugs that prevented the code
from compiling/running correctly as a library, (b) additions needed to get
data back into R (trajectories, grid snapshots, reproducible seeding), and
(c) one numerical-stability fix in the mean-field integrator, described
below.

**Important caveat:** I do not have the original scripts/parameter grids
used to generate the paper's published Figures 2-13 -- only the model
engine. The `examples/` scripts reproduce the *type* of analysis behind
those figures, using parameter values taken directly from the manuscript
(Table 1, and the Fig. 2 caption), and the resulting dynamics are
qualitatively consistent with what the paper describes (see "Validation"
below) -- but they are not guaranteed to be pixel-identical reproductions
of any specific published panel.

## Bugs fixed relative to the shipped `ForestFireR_0.1.0.tar.gz`

1. **`GridMap`/`ReverseGridMap` were used throughout the code but never
   defined** (their declarations in the header comments were just
   comments, not real declarations). The package didn't compile. Added the
   row-major mapping `k = ix*L + jx` that all call sites are consistent
   with.
2. **`run_spatial_model`'s Rcpp export took `double` arguments but the
   internal functions require `long double&`** -- another compile error.
   Fixed by having every export copy its `double` arguments into the
   model's `long double` globals before calling the internal functions.
3. **The grid (`S`) and rate-matrix (`A`) arrays were only ever allocated
   inside `main()`**, which is never called when the code is built as a
   library -- so nothing worked even once it compiled. Added
   `AllocateGrids()`, called at the start of every simulation/landscape
   function, which also lets grid size `L` be set per-call from R instead
   of being hardcoded to 200.
4. **The RNG was only ever seeded inside `main()`** (`srand(time(NULL))`,
   with 1-second resolution) -- called from a library, `rand()` would
   produce the same sequence run after run, and repeated calls within the
   same second would collide. Every simulation function now takes an
   explicit `seed` argument (reproducible) or derives one from a
   high-resolution clock plus a call counter when `seed = NULL`.
5. **The deterministic mean-field integrator (`MeanField`/`zRK4`) diverges
   to `NaN`** under the paper's own default `lambda_F0 = 10^6`: it's an
   explicit fixed-step (h=0.01) RK4 applied to a stiff ODE (the fire
   compartment relaxes on a timescale of `1/lambda_F0`, seven orders of
   magnitude faster than the vegetation dynamics it's coupled to). This is
   not something introduced here -- the original `MeanField()` has the same
   fixed h=0.01 step and would do the same. Fixed with an IMEX
   (implicit-explicit) splitting: the fire compartment's own linear ODE is
   integrated with its closed-form (exponential) solution each step
   (unconditionally stable for any step size), and the comparatively slow
   vegetation compartments are integrated explicitly with the same RK4 as
   before, holding the fire density at that just-updated value. The
   vegetation/fire equations themselves are untouched. As a defensive
   measure (relevant only at parameter corners well outside Table 1's
   calibrated ranges), all state variables are also clamped to `[0, 1]` and
   checked for `NaN`/`Inf` after every step.
6. **`generate_landscape()` (the Appendix B pair-correlation landscape
   generator, `InitialConditionNonHomogeneous` in the C++ source) could hang
   indefinitely** for common `(density2, p)` combinations -- e.g.
   `density2 = 0.9` with almost any `p`, or `density2 >= 0.4` with
   `p = 0.5` (exactly the regime Appendix B's own Fig. A1 uses). Two
   compounding causes:
   - The inner search for a next candidate site (an "isolated" state-1 site,
     or a state-1 site adjacent to the growing cluster) is unbounded
     rejection sampling; at high target densities or clustered `p`, sites
     satisfying the strict condition become rare-to-nonexistent, so the
     search never terminates.
   - Once bounded (see below), the two searches can end up returning the
     *same* site for both candidates. The original code compared them with
     two independent `if (Dp1 < Dp2)` / `if (Dp2 < Dp1)` statements, so an
     exact tie placed nothing and the outer growth loop made zero progress
     that iteration -- confirmed by instrumentation to spin forever with
     `cluster.size()` frozen at a fixed value.

   Fixed with an attempt-limit-plus-relaxation strategy (this changes
   behavior at the specific parameter combinations that used to hang; it
   does not affect any case that previously terminated): each inner search
   is capped at 2000 attempts, after which it falls back to a relaxed
   selection (any random state-1 site, or a direct scan for a
   cluster-adjacent one) and prints a one-time `Rcpp::Rcout` warning
   explaining that the result may deviate slightly from a strict Appendix B
   run for that landscape. Independently, the two-candidate comparison was
   changed from two independent `if`s to `if (Dp1 < Dp2) {...} else {...}`,
   guaranteeing exactly one site is placed per outer iteration regardless of
   ties. Verified against the full previously-hanging matrix (`L` in
   `{20,30,50,70,100,200}`, `density2` in `{0.1,...,0.9}`, `p` in
   `{0.5,0.6,0.9,1}`): all combinations now complete in well under 2
   seconds.

## New capabilities

- Optional trajectory recording (`record_dt` argument) for both the
  spatial and mean-field models, instead of only final densities.
- Optional initial/final grid snapshots (`record_grid = TRUE`) for the
  spatial model.
- A well-mixed *stochastic* counterpart (`simulate_mean_field_stochastic()`)
  to the deterministic mean-field model, for separating the effect of
  spatial structure from demographic stochasticity.
- `generate_landscape()` to build/inspect an initial spatial configuration
  on its own.
- R-facing functions use the paper's own notation (`xi_nat`, `xi_inv`,
  `eta_inv`, `p`) with `xi2lambda()`/`eta2lambda()` handling the conversion
  to the underlying rate constants, and `default_rates()` gives Table 1's
  baseline values.

## Validation

With Fig. 2's stated parameters (`xi_inv = eta_inv = 0.6`, `xi_nat = 0.5`,
10% initial invader density, heterogeneous distribution, T up to 200
years), both `simulate_mean_field()` and `simulate_spatial()` reproduce the
qualitative pattern described in the caption: native density declines
while invader density grows over the 200-year window (mean-field: native
0.90 -> 0.08, invader 0.10 -> 0.61). Sweeping `xi_inv`/`eta_inv` (see
`examples/03_phase_diagram_xi_eta_sweep.R`) reproduces a sharp
continuous-looking transition to native extinction, consistent with the
paper's description of the model's main qualitative result. The
heterogeneity example (`examples/02_landscape_heterogeneity.R`) reproduces
the paper's finding that a segregated/clustered invader (`p = 1`) drives
faster native decline than a randomly dispersed one (`p = 0`).

Mass conservation (`n1+n2+n3+n4+n5 = 1`) was checked to hold to floating
precision for both the spatial and mean-field engines.

## Performance

The spatial model's runtime is dominated by the number of Gillespie events,
which scales with grid size and, once a fire is actively burning, with
`lambda_F0` itself (very short waiting times between fire-related events).
Rough measurements in a modest sandboxed environment:

| L   | T   | wall time |
|-----|-----|-----------|
| 50  | 20  | ~0.1 s    |
| 100 | 10  | ~0.75 s   |
| 200 | 20  | ~19 s     |
| 200 | 100 | ~2+ min (longer if/when fire is actively spreading) |

For interactive exploration, keep `L` around 50-100. Reserve `L = 200`
(paper resolution) for longer, non-interactive/background runs, and use
`record_dt = -1` (final state only) for parameter sweeps where you don't
need the full trajectory.

## Layout

- `src/forest_fire_model.cpp` -- the original model engine (fixed as above)
  plus the new recording-capable functions and the Rcpp export layer.
- `R/rates.R` -- `xi2lambda()`, `eta2lambda()`, `default_rates()`.
- `R/simulate.R` -- `generate_landscape()`, `simulate_spatial()`,
  `simulate_mean_field()`, `simulate_mean_field_stochastic()`, `set_seed()`.
- `examples/` -- runnable scripts illustrating trajectories, landscape
  heterogeneity, and a parameter-sweep phase diagram.
