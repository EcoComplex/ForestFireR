#' ForestFireR: Spatial and Mean-Field Forest-Fire Invasion Model
#'
#' R/Rcpp interface to the spatial stochastic and mean-field forest-fire
#' invasion model from "Dynamics of invasion in forest-fire models".
#'
#' Start with [generate_landscape()] to build an initial spatial
#' configuration, [simulate_spatial()] to run the full spatially-explicit
#' Gillespie simulation, [simulate_mean_field()] for the deterministic
#' well-mixed limit, and [simulate_mean_field_stochastic()] for its
#' stochastic (well-mixed, no spatial structure) counterpart.
#' [default_rates()] and [xi2lambda()] help set up rate parameters using
#' a compact dimensionless notation (baseline values, and the
#' \eqn{\xi}/\eqn{\eta} dimensionless spread/regrowth probabilities).
#'
#' See the `examples/` directory shipped with the package source for
#' worked scripts illustrating the model's main behaviors (an invasion
#' trajectory, a landscape-heterogeneity sweep, a fire-spread/regrowth
#' phase diagram). If you use this package, please cite the paper it
#' implements -- see the README's Citation section.
#'
#' @docType package
#' @name ForestFireR
#' @useDynLib ForestFireR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
