.resolve_seed <- function(seed) {
  if (is.null(seed)) return(-1L)
  as.integer(seed)
}

#' Set the underlying C RNG seed directly
#'
#' The simulation engine uses the C `rand()`/`srand()` generator (as the
#' original model code did), which is independent of R's own RNG. Every
#' `simulate_*()`/`generate_landscape()` call below takes its own `seed`
#' argument (preferred, for reproducible individual runs); call
#' `set_seed()` directly only if you need to seed once and then make several
#' *unseeded* calls in a reproducible sequence.
#'
#' @param seed Integer seed.
#' @export
set_seed <- function(seed) {
  set_seed_cpp(as.integer(seed))
}

#' Generate an initial spatial landscape
#'
#' Places the invader species on a grid of natives according to the
#' initial invader density and the spatial heterogeneity parameter `p`
#' (called \eqn{p} in the paper, `q22` internally), following the
#' pair-correlation landscape-generation algorithm described in Appendix B.
#' `p = 0` disperses the invader uniformly at random; `p = 1` aggregates it
#' into a single compact cluster; intermediate values interpolate between
#' the two.
#'
#' @param L Grid side length (grid has `L*L` sites).
#' @param density2 Initial invader density, in `(0, 1)`.
#' @param p Spatial heterogeneity parameter, in `[0, 2]` (0 = random,
#'   1 = neutral/no preference in the pairwise statistic, 2 = maximally
#'   segregated/clustered -- see Appendix B).
#' @param seed Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.
#' @return An `L`x`L` integer matrix: `1` = native vegetation, `2` = invader.
#' @examples
#' \dontrun{
#' g <- generate_landscape(L = 100, density2 = 0.1, p = 1, seed = 1)
#' image(g)
#' }
#' @export
generate_landscape <- function(L = 100, density2 = 0.1, p = 1, seed = NULL) {
  generate_landscape_cpp(as.integer(L), density2, p, .resolve_seed(seed))
}

#' Simulate the spatially-explicit stochastic model
#'
#' Runs the model's Gillespie stochastic simulation algorithm on an `L`x`L`
#' periodic (or bounded) grid, starting from an invader landscape generated
#' the same way as [generate_landscape()]. This is the direct R-friendly
#' interface to the paper's spatial model (Section 2, Eq. \[the CTMC /
#' master equation\]).
#'
#' @param T Simulation time horizon (years).
#' @param L Grid side length (grid has `L*L` sites). Runtime scales roughly
#'   with `L^2 * T`; for interactive exploration keep `L` around 50-100,
#'   reserve `L = 200` (paper resolution) for longer background runs (a
#'   single `L = 200`, `T = 100` run with default rates takes on the order
#'   of a couple of minutes in a typical environment).
#' @param density2 Initial invader density.
#' @param p Spatial heterogeneity of the initial landscape (see
#'   [generate_landscape()]).
#' @param xi_nat,xi_inv Fire-spread probabilities for native/invader
#'   (see [xi2lambda()]). Table 1 fixes `xi_inv = 0.6` throughout the main
#'   results and sweeps `xi_nat` roughly across `[0.1, 0.9]`; `xi_nat = 0.5`
#'   (the default here) is the value used in the paper's Fig. 2.
#' @param eta_inv Invader post-fire regrowth probability (see
#'   [xi2lambda()]); default `0.6` (paper's Fig. 2 value), also swept
#'   roughly across `[0.1, 0.9]` for other figures.
#' @param L_01,L_02,L_12,L_21,L_30,Lig_13,Lig_23 Base rate constants; see
#'   [default_rates()] for Table 1 defaults. Override individual entries to
#'   explore parameter space, e.g. `L_21 = 0.02`.
#' @param periodic Use periodic boundary conditions (paper default `TRUE`).
#' @param seed Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.
#' @param record_dt Trajectory sampling interval, in the same time units as
#'   `T`. `-1` (default) records only the final state (fastest, use this
#'   for parameter sweeps). `0` records every single event (dense, memory
#'   heavy, only for small/short runs). A positive value samples the
#'   trajectory roughly every `record_dt` time units (use this to reproduce
#'   trajectory figures, e.g. `record_dt = T/200`).
#' @param record_grid If `TRUE`, also return the initial and final grid
#'   snapshots (`L`x`L` integer matrices) -- needed for the spatial-pattern
#'   figures.
#' @param capture_fire_snapshots If `TRUE`, also return every grid snapshot
#'   captured while the fire compartment was actively non-empty (`n3 > 0`),
#'   throttled to at most one capture per `fire_snapshot_min_gap` time
#'   units. Use this to find an illustrative "fire in progress" grid near a
#'   target time (e.g. Fig. 2's panels at `t=76`, `t=120`): `record_grid`'s
#'   start/end snapshots essentially never show fire, because Table 1's
#'   fire-spread/fire-extinction rates (~1e6/year, hours-scale) resolve any
#'   given outbreak in a sliver of simulated time invisible to a snapshot
#'   at an arbitrary fixed time -- an outbreak has to be searched for
#'   explicitly, not hoped for. After the run, pick the returned
#'   `fire_snapshot_times` entry closest to your target time and use the
#'   matching `fire_snapshot_grids` entry; there is no guarantee one exists
#'   near any particular time (ignition is itself a rare Poisson process --
#'   `lambda^ig = 1e-4`/year per invader-occupied site, Table 1), so check
#'   `length(fire_snapshot_times)` and how close the nearest one actually
#'   is before trusting it as illustrative of that moment.
#' @param fire_snapshot_min_gap Minimum simulated-time gap (years) between
#'   two fire snapshots; keeps a single outbreak from filling the returned
#'   list with near-duplicate frames. Default `0.001` (~9 hours) is well
#'   under a typical outbreak's total duration but well above the time
#'   between individual Gillespie events during one.
#' @return A list with the final `time_sim` and densities `n0`..`n5`
#'   (`n0` = total empty, `n1` = native, `n2` = invader, `n3` = fire,
#'   `n4` = empty-after-vegetation-death, `n5` = empty-after-fire), plus
#'   `trajectory` (a data frame, if `record_dt >= 0`), `initial_grid` /
#'   `final_grid` (if `record_grid = TRUE`), and `fire_snapshot_times` (a
#'   numeric vector) / `fire_snapshot_grids` (a list of `L`x`L` integer
#'   matrices, same length and order) if `capture_fire_snapshots = TRUE`.
#' @examples
#' \dontrun{
#' r <- simulate_spatial(T = 20, L = 60, density2 = 0.05, p = 1,
#'                        xi_inv = 0.6, eta_inv = 0.6, record_dt = 0.5)
#' plot(r$trajectory$time, r$trajectory$n1, type = "l")
#' }
#' @export
simulate_spatial <- function(T, L = 100, density2 = 0.1, p = 1,
                              xi_nat = 0.5, xi_inv = 0.6, eta_inv = 0.6,
                              L_01 = 0.03, L_02 = 0.03, L_12 = 0.005, L_21 = 0.01,
                              L_30 = 1e6, Lig_13 = 0, Lig_23 = 1e-4,
                              periodic = TRUE, seed = NULL,
                              record_dt = -1, record_grid = FALSE,
                              capture_fire_snapshots = FALSE, fire_snapshot_min_gap = 0.001) {
  Lsp_13 <- xi2lambda(xi_nat, L_30)
  Lsp_23 <- xi2lambda(xi_inv, L_30)
  Lrg_02 <- eta2lambda(eta_inv, L_01)

  simulate_spatial_cpp(
    T = T, Lgrid = as.integer(L), density2 = density2, p = p,
    L_01_ = L_01, L_02_ = L_02, L_10_ = 0, L_20_ = 0,
    L_12_ = L_12, L_21_ = L_21,
    L_30_ = L_30, Lig_13_ = Lig_13, Lig_23_ = Lig_23,
    Lsp_13_ = Lsp_13, Lsp_23_ = Lsp_23,
    Lrg_01_ = 0, Lrg_02_ = Lrg_02,
    Lr_01_ = 0, Lr_02_ = 0, Lr_12_ = 0, Lr_21_ = 0,
    periodic = periodic, seed = .resolve_seed(seed),
    record_dt = record_dt, record_grid = record_grid,
    capture_fire_snapshots = capture_fire_snapshots,
    fire_snapshot_min_gap = fire_snapshot_min_gap
  )
}

#' Simulate the deterministic mean-field model
#'
#' Integrates the mean-field ODE system (paper Eq. 12) describing the
#' well-mixed limit of the model, and returns the full trajectory as a
#' data frame. Internally uses a semi-implicit (IMEX) step for the fire
#' compartment so that it stays numerically stable even under the paper's
#' default `L_30 = 1e6` (an explicit fixed-step RK4 on this term alone
#' diverges for any practical step size) -- see the package's implementation
#' notes (`vignette("forest-fire-model")` or `NEWS.md`) for details. The
#' underlying vegetation/fire equations are unchanged from the original
#' model.
#'
#' @param T Simulation time horizon.
#' @param n1_0,n2_0,n3_0,n4_0,n5_0 Initial densities for native, invader,
#'   fire, empty-after-vegetation-death and empty-after-fire respectively
#'   (should sum to 1).
#' @inheritParams simulate_spatial
#' @param record_dt Trajectory sampling interval (must be `> 0`; default
#'   `T/200`).
#' @return A data frame with columns `time`, `n0`..`n5`.
#' @examples
#' \dontrun{
#' mf <- simulate_mean_field(T = 200, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05)
#' plot(mf$time, mf$n1, type = "l", ylim = c(0,1))
#' lines(mf$time, mf$n2, col = "red")
#' }
#' @export
simulate_mean_field <- function(T, n1_0 = 0.9, n2_0 = 0.05, n3_0 = 0, n4_0 = 0.05, n5_0 = 0,
                                 xi_nat = 0.5, xi_inv = 0.6, eta_inv = 0.6,
                                 L_01 = 0.03, L_02 = 0.03, L_12 = 0.005, L_21 = 0.01,
                                 L_30 = 1e6, Lig_13 = 0, Lig_23 = 1e-4,
                                 record_dt = NULL) {
  if (is.null(record_dt)) record_dt <- T / 200
  Lsp_13 <- xi2lambda(xi_nat, L_30)
  Lsp_23 <- xi2lambda(xi_inv, L_30)
  Lrg_02 <- eta2lambda(eta_inv, L_01)

  simulate_mean_field_cpp(
    T = T, n1_0 = n1_0, n2_0 = n2_0, n3_0 = n3_0, n4_0 = n4_0, n5_0 = n5_0,
    L_01_ = L_01, L_02_ = L_02, L_10_ = 0, L_20_ = 0,
    L_12_ = L_12, L_21_ = L_21, L_30_ = L_30,
    Lig_13_ = Lig_13, Lig_23_ = Lig_23, Lsp_13_ = Lsp_13, Lsp_23_ = Lsp_23,
    Lrg_01_ = 0, Lrg_02_ = Lrg_02,
    Lr_01_ = 0, Lr_02_ = 0, Lr_12_ = 0, Lr_21_ = 0,
    record_dt = record_dt
  )
}

#' Simulate the well-mixed stochastic (non-spatial) model
#'
#' Runs the well-mixed Gillespie simulation (mass-action, no spatial
#' structure) at a given population size `N`, the stochastic counterpart
#' to [simulate_mean_field()]. Useful for isolating the effect of spatial
#' structure from the effect of demographic stochasticity (paper Section
#' 3.3/3.4).
#'
#' @param T Simulation time horizon.
#' @param N Well-mixed population size (analogous to `L*L` in the spatial
#'   model, but with no spatial structure).
#' @param n1_0,n2_0,n3_0,n4_0,n5_0 Initial densities (should sum to 1).
#' @inheritParams simulate_spatial
#' @param record_dt Trajectory sampling interval. `-1` = final state only,
#'   `0` = every event, `>0` = sample interval. Default `T/200`.
#' @return A data frame with columns `time`, `n0`..`n5`. Note: the CTMC can
#'   reach a true absorbing state (e.g. all-native with no empty sites and
#'   no residual ignition source) before `T`, in which case the trajectory
#'   simply stops there.
#' @export
simulate_mean_field_stochastic <- function(T, N = 10000, n1_0 = 0.9, n2_0 = 0.05, n3_0 = 0,
                                            n4_0 = 0.05, n5_0 = 0,
                                            xi_nat = 0.5, xi_inv = 0.6, eta_inv = 0.6,
                                            L_01 = 0.03, L_02 = 0.03, L_12 = 0.005, L_21 = 0.01,
                                            L_30 = 1e6, Lig_13 = 0, Lig_23 = 1e-4,
                                            seed = NULL, record_dt = NULL) {
  if (is.null(record_dt)) record_dt <- T / 200
  Lsp_13 <- xi2lambda(xi_nat, L_30)
  Lsp_23 <- xi2lambda(xi_inv, L_30)
  Lrg_02 <- eta2lambda(eta_inv, L_01)

  simulate_mean_field_stochastic_cpp(
    T = T, sysN = N, n1_0 = n1_0, n2_0 = n2_0, n3_0 = n3_0, n4_0 = n4_0, n5_0 = n5_0,
    L_01_ = L_01, L_02_ = L_02, L_10_ = 0, L_20_ = 0,
    L_12_ = L_12, L_21_ = L_21, L_30_ = L_30,
    Lig_13_ = Lig_13, Lig_23_ = Lig_23, Lsp_13_ = Lsp_13, Lsp_23_ = Lsp_23,
    Lrg_01_ = 0, Lrg_02_ = Lrg_02,
    seed = .resolve_seed(seed), record_dt = record_dt
  )
}

#' Named grid-state codes
#'
#' The integer codes used in the `L`x`L` grids returned/accepted by
#' [generate_landscape()], [generate_landscape_layers()] and
#' [simulate_spatial_from_grid()]: `1` = native vegetation (paper's
#' \eqn{V^{nat}}), `2` = invader vegetation (\eqn{V^{inv}}), `3` = active
#' fire (\eqn{F}), `4` = unoccupied, never-yet-burned space
#' (\eqn{\emptyset}), `5` = unoccupied, post-fire space (\eqn{\emptyset^F},
#' the only one of the two empty states with a regrowth-rate advantage for
#' the invader, via `eta_inv`). State `0` never appears -- it's an unused
#' index in the underlying engine's state numbering.
#' @export
FF_STATE <- c(native = 1L, invader = 2L, fire = 3L, empty = 4L, empty_postfire = 5L)

#' Build a landscape from several placement layers
#'
#' Generalizes [generate_landscape()] to construct initial conditions it
#' can't express on its own: more than one pattern placed independently
#' over a shared background, or a small pattern placed over an
#' already-uniform domain. The whole grid is first filled with
#' `fill_state`; then each layer in `layers`, in order, grows its `pattern`
#' state to `density * L*L` cells (a fraction of the WHOLE grid, same
#' convention as [generate_landscape()]'s `density2`) with heterogeneity
#' `p` (same meaning as there), drawing candidate cells only from whatever
#' `background` state is still left after the layers before it. If a
#' layer's target can't be reached because its background ran out (e.g.
#' the density budget across layers doesn't leave enough), the engine
#' prints a one-time notice and stops that layer early rather than
#' hanging -- check the printed cell counts against what you asked for if
#' you see this.
#'
#' Used to build the initial conditions for Figs. 8-9 (5% native + 5%
#' invader, placed independently, over an all-`empty_postfire` domain) and
#' Figs. 10-11 (1% fire placed over an all-`invader` domain) -- see
#' `Source/R_reproduce_figures/` in the manuscript repository for the full
#' reproduction scripts.
#'
#' @param L Grid side length.
#' @param fill_state The state ([FF_STATE]) the whole grid starts as,
#'   before any layer is placed.
#' @param layers A list of layers, each itself a list with elements
#'   `background`, `pattern` (states from [FF_STATE]), `density` (fraction
#'   of `L*L`) and `p` (heterogeneity, as in [generate_landscape()]).
#' @param seed Integer RNG seed, or `NULL` for a fresh, non-reproducible seed.
#' @return An `L`x`L` integer matrix using [FF_STATE]'s codes.
#' @examples
#' \dontrun{
#' # 5% native + 5% invader over a post-fire-empty domain (Figs. 8-9)
#' g <- generate_landscape_layers(L = 100, fill_state = FF_STATE["empty_postfire"],
#'   layers = list(
#'     list(background = FF_STATE["empty_postfire"], pattern = FF_STATE["native"],
#'          density = 0.05, p = 0.8),
#'     list(background = FF_STATE["empty_postfire"], pattern = FF_STATE["invader"],
#'          density = 0.05, p = 0.8)
#'   ), seed = 1)
#'
#' # 1% fire over an all-invader domain (Figs. 10-11)
#' g2 <- generate_landscape_layers(L = 100, fill_state = FF_STATE["invader"],
#'   layers = list(
#'     list(background = FF_STATE["invader"], pattern = FF_STATE["fire"],
#'          density = 0.01, p = 0.8)
#'   ), seed = 1)
#' }
#' @export
generate_landscape_layers <- function(L = 100, fill_state, layers, seed = NULL) {
  background_states <- vapply(layers, function(l) as.integer(l$background), integer(1))
  pattern_states <- vapply(layers, function(l) as.integer(l$pattern), integer(1))
  densities <- vapply(layers, function(l) as.numeric(l$density), numeric(1))
  ps <- vapply(layers, function(l) as.numeric(l$p), numeric(1))

  generate_landscape_layers_cpp(
    Lgrid = as.integer(L), fill_state = as.integer(fill_state),
    background_states = background_states, pattern_states = pattern_states,
    densities = densities, ps = ps, seed = .resolve_seed(seed)
  )
}

#' Simulate the spatially-explicit model from a supplied landscape
#'
#' Like [simulate_spatial()], but starts from an initial grid you supply
#' (typically from [generate_landscape_layers()]) instead of building one
#' internally from `density2`/`p`. Use this for any initial condition
#' [generate_landscape()] can't express on its own -- e.g. Figs. 8-11,
#' which also need some of the rate constants below zeroed out to match
#' the reduced reaction sets those figures use (see
#' `Source/R_reproduce_figures/` in the manuscript repository): Figs. 8-9
#' disable fire entirely (`xi_nat = 0, xi_inv = 0, Lig_23 = 0`) and vegetation
#' replacement (`L_12 = 0, L_21 = 0`), keeping only colonization/regrowth;
#' Figs. 10-11 disable colonization/regrowth (`L_01 = 0, L_02 = 0, eta_inv = 0`)
#' and native-side reactions (irrelevant since there's no native vegetation
#' in that initial condition), keeping only fire spread/burnout.
#'
#' @param T Simulation time horizon.
#' @param initial_grid An `L`x`L` integer matrix using [FF_STATE]'s codes
#'   (e.g. from [generate_landscape_layers()]).
#' @param check_extinction The engine's usual early stop -- treat the
#'   native species as extinct and end the run as soon as its density
#'   drops below `0.0001` -- assumes `initial_grid` actually contains
#'   native vegetation to begin with. Set `FALSE` for an `initial_grid`
#'   that has none by design (e.g. Figs. 10-11's fire-over-invader
#'   domain), or that early stop triggers on the very first step. Leave at
#'   the default `TRUE` whenever native vegetation is present (matches
#'   [simulate_spatial()], which always uses the check).
#' @inheritParams simulate_spatial
#' @export
simulate_spatial_from_grid <- function(T, initial_grid,
                                        xi_nat = 0.5, xi_inv = 0.6, eta_inv = 0.6,
                                        L_01 = 0.03, L_02 = 0.03, L_12 = 0.005, L_21 = 0.01,
                                        L_30 = 1e6, Lig_13 = 0, Lig_23 = 1e-4,
                                        periodic = TRUE, seed = NULL,
                                        record_dt = -1, record_grid = FALSE,
                                        check_extinction = TRUE,
                                        capture_fire_snapshots = FALSE, fire_snapshot_min_gap = 0.001) {
  Lsp_13 <- xi2lambda(xi_nat, L_30)
  Lsp_23 <- xi2lambda(xi_inv, L_30)
  Lrg_02 <- eta2lambda(eta_inv, L_01)

  storage.mode(initial_grid) <- "integer"

  simulate_spatial_from_grid_cpp(
    T = T, initial_grid = initial_grid,
    L_01_ = L_01, L_02_ = L_02, L_10_ = 0, L_20_ = 0,
    L_12_ = L_12, L_21_ = L_21,
    L_30_ = L_30, Lig_13_ = Lig_13, Lig_23_ = Lig_23,
    Lsp_13_ = Lsp_13, Lsp_23_ = Lsp_23,
    Lrg_01_ = 0, Lrg_02_ = Lrg_02,
    Lr_01_ = 0, Lr_02_ = 0, Lr_12_ = 0, Lr_21_ = 0,
    periodic = periodic, seed = .resolve_seed(seed),
    record_dt = record_dt, record_grid = record_grid,
    check_extinction = check_extinction,
    capture_fire_snapshots = capture_fire_snapshots,
    fire_snapshot_min_gap = fire_snapshot_min_gap
  )
}
