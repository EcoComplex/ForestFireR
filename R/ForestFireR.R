#' ForestFireR: Spatial and Mean-Field Forest-Fire Invasion Model
#'
#' Simulation of a two-species vegetation system in which fire mediates the
#' competition between a native species and an invader, on an explicit
#' lattice or in its well-mixed limit. The engine is written in C++ and
#' exposed to R through Rcpp.
#'
#' See [ForestFireR-model] for the definition of the model: states,
#' transitions, parameters and the mean-field equations.
#'
#' To run it: [generate_landscape()] builds an initial spatial
#' configuration, [simulate_spatial()] runs the lattice Gillespie
#' simulation, [simulate_mean_field()] integrates the deterministic
#' well-mixed limit and [simulate_mean_field_stochastic()] simulates the
#' well-mixed stochastic process at a finite population size.
#' [default_rates()] gives the baseline parameter set and [xi2lambda()]
#' converts the dimensionless spread/regrowth probabilities into rate
#' constants. [generate_landscape_layers()] and
#' [simulate_spatial_from_grid()] cover initial conditions that
#' [generate_landscape()] cannot express, and [simulate_spatial_movie()]
#' renders a run as a video.
#'
#' Worked scripts (an invasion trajectory, a landscape-clustering
#' comparison, a fire-spread/regrowth phase diagram) are installed with the
#' package:
#'
#' ```r
#' list.files(system.file("examples", package = "ForestFireR"))
#' ```
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib ForestFireR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
