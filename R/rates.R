#' Convert a spread or regrowth probability to a rate constant
#'
#' The model parametrizes fire spread and post-fire regrowth with two
#' dimensionless probabilities: \eqn{\xi}, the probability that a burning
#' neighbour ignites a given vegetation type before the fire itself burns
#' out, and \eqn{\eta}, the probability that an unoccupied post-fire site
#' regrows into a given species before it is recolonised some other way.
#' Both convert to a rate constant through the same relation, obtained by
#' setting \eqn{\xi = \lambda^{sp}/(\lambda^{sp}+\lambda^{F\emptyset})} (and
#' analogously for \eqn{\eta}) and solving for the rate:
#' \deqn{\lambda = \frac{\xi}{1-\xi}\,\lambda_{ref}}
#'
#' @param xi Probability in `[0, 1)` (named `xi` for fire spread, `eta` for
#'   regrowth -- both use this same formula).
#' @param lambda_ref The reference rate the probability is defined against
#'   (`L_30`, i.e. \eqn{\lambda^{F\emptyset}}, for fire spread; `L_01`,
#'   i.e. \eqn{\lambda_{\emptyset nat}}, for regrowth -- see Details).
#' @return The corresponding rate constant (same units as `lambda_ref`).
#'
#' @details Fire spread: `xi2lambda(xi_nat, L_30)` gives `Lsp_13` (native
#'   catches fire from a burning neighbour), `xi2lambda(xi_inv, L_30)` gives
#'   `Lsp_23` (invader catches fire from a burning neighbour).
#'   \eqn{\xi = 0.5} is the point where spread and burnout are equally fast:
#'   below it a fire usually dies where it started, above it a fire usually
#'   propagates.
#'
#'   Regrowth: both species' post-fire regrowth probabilities convert the
#'   same way, using `L_01` (not `L_02`) as the reference rate for both:
#'   `eta2lambda(eta_nat, L_01)` gives `Lrg_01` (native long-range
#'   regrowth), `eta2lambda(eta_inv, L_01)` gives `Lrg_02` (invader
#'   long-range regrowth). Setting `L_01 = 0` therefore switches off both
#'   regrowth channels as well as ordinary native colonisation.
#'
#'   The baseline parameter set leaves the native species without a regrowth
#'   advantage (`eta_nat = 0`); `eta_nat` exists so the two species can be
#'   parametrized symmetrically.
#' @examples
#' # A fire-spread probability of 0.6 against a burnout rate of 1e6/year
#' xi2lambda(0.6, 1e6)
#'
#' # xi = 0.5 returns exactly the reference rate (spread as fast as burnout)
#' xi2lambda(0.5, 1e6)
#'
#' # Regrowth uses L_01 as its reference rate
#' eta2lambda(0.6, 0.03)
#'
#' # The conversion is what the simulate functions apply internally
#' r <- default_rates()
#' xi2lambda(r$xi_inv, r$L_30)
#' @seealso [default_rates()], [ForestFireR-model]
#' @export
xi2lambda <- function(xi, lambda_ref) {
  if (any(xi < 0 | xi >= 1)) stop("xi/eta must be in [0, 1)")
  (xi * lambda_ref) / (1 - xi)
}

#' @rdname xi2lambda
#' @export
eta2lambda <- function(xi, lambda_ref) xi2lambda(xi, lambda_ref)

#' Baseline rate constants
#'
#' Returns the model's baseline parameter set as a named list, ready to be
#' passed to [simulate_spatial()], [simulate_mean_field()] or
#' [simulate_mean_field_stochastic()] via `do.call()`, or used as a starting
#' point whose individual entries can be overridden. These are the values
#' every simulate function uses by default. Time is in years.
#'
#' Under this parametrization the native species replaces the invader twice
#' as fast as the reverse, so without fire the invader cannot spread; its
#' advantages are spontaneous ignition, a higher fire-spread probability and
#' post-fire regrowth. See [ForestFireR-model] for the mechanism and
#' [xi2lambda()] for the dimensionless probabilities.
#'
#' @return A named list with the vegetation turnover rates
#'   (`L_01`, `L_02`, `L_12`, `L_21`), the vegetation death rates
#'   (`L_10`, `L_20`, both 0), the fire burnout rate (`L_30`,
#'   i.e. \eqn{\lambda^{F\emptyset}}), the two species' spontaneous ignition
#'   rates (`Lig_13`, `Lig_23`), and the dimensionless `xi_nat`, `xi_inv`,
#'   `eta_nat`, `eta_inv` (see [xi2lambda()]).
#' @examples
#' str(default_rates())
#'
#' # Run with the baseline set, raising only the invader's regrowth advantage
#' pars <- default_rates()
#' pars$eta_inv <- 0.8
#' pars[c("L_10", "L_20")] <- NULL          # not arguments of the simulate functions
#' mf <- do.call(simulate_mean_field, c(list(T = 200), pars))
#' tail(mf, 1)
#' @seealso [xi2lambda()], [ForestFireR-model]
#' @export
default_rates <- function() {
  list(
    L_01 = 0.03,     # empty -> native,  per native neighbour (lambda_{0,nat})
    L_02 = 0.03,     # empty -> invader, per invader neighbour (lambda_{0,inv})
    L_10 = 0,        # native  -> empty  (no background mortality in this model)
    L_20 = 0,        # invader -> empty  (idem)
    L_12 = 0.005,    # native  -> invader, competitive replacement (lambda_{nat,inv})
    L_21 = 0.01,     # invader -> native,  competitive replacement (lambda_{inv,nat})
    L_30 = 1e6,      # fire -> empty_postfire, i.e. lambda^{F,empty} (fire burnout)
    Lig_13 = 0,      # native  spontaneous ignition (the native is not fire-adapted)
    Lig_23 = 1e-4,   # invader spontaneous ignition, lambda^{ig}
    xi_nat = 0.5,    # native  fire-spread probability (the spread/burnout threshold)
    xi_inv = 0.6,    # invader fire-spread probability (moderately above the threshold)
    eta_nat = 0,     # native  post-fire regrowth probability (no native advantage)
    eta_inv = 0.6    # invader post-fire regrowth probability (a moderate advantage)
  )
}
