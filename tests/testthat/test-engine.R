# The "tree" and "linear" engines select events from the same cumulative
# rate intervals, but they add the rates in a different order (the tree
# reads the total from its root; the linear scan keeps a running total).
# Rounding therefore differs slightly, and sooner or later a random draw
# falls on an interval boundary and the two engines pick different events.
# From then on the runs diverge, as any stochastic simulation does. When
# that happens depends on the precision of `long double` (80 bits on x86,
# 64 bits on Apple Silicon), so an exact match is only required over a
# short horizon, and equivalence over longer runs is tested in
# distribution across seeds.

engine_test_grid <- function() {
  set.seed(1)
  matrix(sample(1:2, 40 * 40, TRUE, c(0.3, 0.7)), 40)
}

engine_run <- function(engine, T, seed, g, record_dt = 1) {
  simulate_spatial_from_grid(T = T, initial_grid = g, xi_nat = 0.1, xi_inv = 0.9,
                             eta_nat = 0.0066, eta_inv = 0.99, Lig_23 = 1 / 1600,
                             seed = seed, record_dt = record_dt,
                             check_extinction = FALSE, engine = engine)
}

cols <- c("n1", "n2", "n3", "n4", "n5")

test_that("tree and linear engines agree over a short horizon for a given seed", {
  g <- engine_test_grid()
  a <- engine_run("linear", T = 5, seed = 7, g = g)
  b <- engine_run("tree",   T = 5, seed = 7, g = g)
  expect_identical(nrow(a$trajectory), nrow(b$trajectory))
  # Allow at most a few sites of difference (1 site = 1/1600).
  max_diff <- max(abs(as.matrix(a$trajectory[, cols]) - as.matrix(b$trajectory[, cols])))
  expect_lte(max_diff, 3 / 1600)
  expect_equal(a$trajectory$time, b$trajectory$time, tolerance = 1e-6)
})

test_that("simulate_spatial accepts both engines and defaults to tree", {
  expect_identical(eval(formals(simulate_spatial)$engine), c("tree", "linear"))
  r1 <- simulate_spatial(T = 5, L = 30, seed = 3, engine = "tree")
  r2 <- simulate_spatial(T = 5, L = 30, seed = 3, engine = "linear")
  for (r in list(r1, r2)) {
    expect_true(all(is.finite(c(r$n1, r$n2))))
    expect_true(all(c(r$n1, r$n2) >= 0 & c(r$n1, r$n2) <= 1))
  }
  expect_error(simulate_spatial(T = 5, L = 30, seed = 3, engine = "fast"))
})

test_that("tree and linear engines give the same dynamics in distribution", {
  skip_on_cran()
  g <- engine_test_grid()
  seeds <- 1:30
  # Per run: time-averaged density of each state over t in [10, 30].
  summarise_run <- function(engine, seed) {
    tr <- engine_run(engine, T = 30, seed = seed, g = g)$trajectory
    colMeans(tr[tr$time >= 10, cols, drop = FALSE])
  }
  lin  <- t(sapply(seeds, function(s) summarise_run("linear", s)))
  tree <- t(sapply(seeds, function(s) summarise_run("tree", s)))
  for (cl in cols) {
    se  <- sqrt(var(lin[, cl]) / length(seeds) + var(tree[, cl]) / length(seeds))
    tol <- max(4 * se, 0.01)
    expect_lte(abs(mean(lin[, cl]) - mean(tree[, cl])), tol,
               label = paste0("|mean difference| for ", cl))
  }
})
