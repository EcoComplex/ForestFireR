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
#' the paper's own notation (Table 1, and the \eqn{\xi}/\eqn{\eta}
#' dimensionless spread/regrowth probabilities).
#'
#' See the `examples/` directory shipped with the package source for
#' scripts that reproduce the paper's main figures.
#'
#' @docType package
#' @name ForestFireR
#' @useDynLib ForestFireR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
