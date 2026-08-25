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
#' @return A list with the final `time_sim` and densities `n0`..`n5`
#'   (`n0` = total empty, `n1` = native, `n2` = invader, `n3` = fire,
#'   `n4` = empty-after-vegetation-death, `n5` = empty-after-fire), plus
#'   `trajectory` (a data frame, if `record_dt >= 0`) and `initial_grid` /
#'   `final_grid` (if `record_grid = TRUE`).
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
                              record_dt = -1, record_grid = FALSE) {
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
    record_dt = record_dt, record_grid = record_grid
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
