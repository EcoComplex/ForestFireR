# The five compartments (native, invader, fire, empty-after-death,
# empty-after-fire) partition every site, so n1+n2+n3+n4+n5 must equal 1
# (within floating-point tolerance) at every recorded time, for every
# engine (spatial, mean-field ODE, well-mixed stochastic) and at every
# point of the trajectory, not just the endpoint. This is a known a
# priori invariant of the model, independent of parameters or seed.

expect_conserved <- function(df, tol = 1e-6) {
  total <- df$n1 + df$n2 + df$n3 + df$n4 + df$n5
  expect_true(all(abs(total - 1) < tol))
}

test_that("spatial model conserves total density (endpoint)", {
  r <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 1, record_dt = -1)
  total <- r$n1 + r$n2 + r$n3 + r$n4 + r$n5
  expect_equal(total, 1, tolerance = 1e-6)
  # n0 is reported as the sum of the two empty sub-states
  expect_equal(r$n0, r$n4 + r$n5, tolerance = 1e-6)
})

test_that("spatial model conserves total density along the whole trajectory", {
  r <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 1, record_dt = 0.1)
  expect_conserved(r$trajectory)
})

test_that("mean-field ODE model conserves total density along the trajectory", {
  mf <- simulate_mean_field(T = 50, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05, record_dt = 1)
  expect_conserved(mf)
})

test_that("well-mixed stochastic model conserves total density along the trajectory", {
  ms <- simulate_mean_field_stochastic(T = 50, N = 2000, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                                        seed = 1, record_dt = 1)
  expect_conserved(ms)
})
