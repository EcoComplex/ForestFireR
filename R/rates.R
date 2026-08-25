#' Convert a spread/regrowth probability to a rate constant
#'
#' The paper parametrizes fire spread and post-fire regrowth with two
#' dimensionless probabilities, \eqn{\xi} (probability that a burning
#' neighbour ignites a given vegetation type before the fire itself burns
#' out) and \eqn{\eta} (probability that an empty-after-fire site regrows
#' into the invader before it is recolonised some other way). Both convert
#' to a rate constant through the same relation, derived from setting
#' \eqn{\xi = \lambda^{sp}/(\lambda^{sp}+\lambda^{F\emptyset})} (and
#' analogously for \eqn{\eta}) and solving for the rate:
#' \deqn{\lambda = \frac{\xi}{1-\xi}\,\lambda_{ref}}
#'
#' @param xi Probability in `[0, 1)` (named `xi` for fire spread, `eta` for
#'   regrowth -- both use this same formula).
#' @param lambda_ref The reference rate the probability is defined against
#'   (`L_30`, i.e. \eqn{\lambda^{F\emptyset}}, for fire spread; `L_01`,
#'   i.e. \eqn{\lambda_{\emptyset V1}}, for regrowth -- see Details).
#' @return The corresponding rate constant (same units as `lambda_ref`).
#'
#' @details Fire spread: `xi2lambda(xi_nat, L_30)` gives `Lsp_13` (native
#'   catches fire from a burning neighbour), `xi2lambda(xi_inv, L_30)` gives
#'   `Lsp_23` (invader catches fire from a burning neighbour).
#'
#'   Regrowth: only the invader's post-fire regrowth advantage `eta_inv` is
#'   defined in the paper (this is intentional -- see the manuscript's
#'   discussion of why an analogous `eta_nat` is not defined). It converts
#'   to `Lrg_02` using `L_01` (not `L_02`) as the reference rate:
#'   `eta2lambda(eta_inv, L_01)` gives `Lrg_02`.
#' @export
xi2lambda <- function(xi, lambda_ref) {
  if (any(xi < 0 | xi >= 1)) stop("xi/eta must be in [0, 1)")
  (xi * lambda_ref) / (1 - xi)
}

#' @rdname xi2lambda
#' @export
eta2lambda <- function(xi, lambda_ref) xi2lambda(xi, lambda_ref)

#' Table 1 default rate constants
#'
#' Returns the baseline transition rates used throughout the paper (Table 1),
#' as a named list ready to be passed to [simulate_spatial()],
#' [simulate_mean_field()] or [simulate_mean_field_stochastic()] via
#' `do.call()`, or used as a starting point that individual entries can be
#' overridden from.
#'
#' @return A named list with the vegetation turnover rates
#'   (`L_01`, `L_02`, `L_12`, `L_21`), the fire persistence rate (`L_30`,
#'   i.e. \eqn{\lambda^{F\emptyset}}), the invader's spontaneous ignition
#'   rate (`Lig_23`), and default dimensionless `xi_nat`, `xi_inv`,
#'   `eta_inv` (see [xi2lambda()]).
#' @export
default_rates <- function() {
  list(
    L_01 = 0.03,     # empty_V -> native   (paper's lambda_{0,V1})
    L_02 = 0.03,     # empty_V -> invader  (lambda_{0,V2})
    L_10 = 0,        # native  -> empty_V  (not used in the paper's model)
    L_20 = 0,        # invader -> empty_V  (not used in the paper's model)
    L_12 = 0.005,    # native  -> invader, competitive replacement (lambda_{12})
    L_21 = 0.01,     # invader -> native,  competitive replacement (lambda_{21})
    L_30 = 1e6,      # fire -> empty_F, i.e. lambda^{F\\emptyset} (fire persistence)
    Lig_13 = 0,      # native  spontaneous ignition (paper: 0, native is not fire-adapted)
    Lig_23 = 1e-4,   # invader spontaneous ignition, lambda^{ig}
    xi_nat = 0.5,    # native  fire-spread probability (Fig. 2 value; paper sweeps this in [~0.1, 0.9])
    xi_inv = 0.6,    # invader fire-spread probability (fixed across the paper's main results)
    eta_inv = 0.6    # invader post-fire regrowth probability (Fig. 2 value; paper sweeps this in [~0.1, 0.9])
  )
}
