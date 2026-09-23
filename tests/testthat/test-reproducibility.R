# The engine seeds its own C RNG explicitly per call (see ?set_seed);
# the same seed and parameters must reproduce the exact same run, not
# just similar aggregate densities. This guards against any accidental
# reliance on ambient/global RNG state or timing-dependent seeding.

test_that("simulate_spatial is exactly reproducible for a fixed seed", {
  rA <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 42, record_dt = -1)
  rB <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 42, record_dt = -1)
  expect_identical(rA$n1, rB$n1)
  expect_identical(rA$n2, rB$n2)
  expect_identical(rA$n3, rB$n3)
  expect_identical(rA$time_sim, rB$time_sim)
})

test_that("simulate_spatial's initial grid is exactly reproducible for a fixed seed", {
  rA <- simulate_spatial(T = 0.001, L = 25, density2 = 0.2, p = 1, seed = 7,
                          record_dt = -1, record_grid = TRUE)
  rB <- simulate_spatial(T = 0.001, L = 25, density2 = 0.2, p = 1, seed = 7,
                          record_dt = -1, record_grid = TRUE)
  expect_identical(rA$initial_grid, rB$initial_grid)
})

test_that("different seeds are (almost always) not identical", {
  # a sanity check in the other direction: two different seeds should not
  # coincidentally produce the exact same trajectory endpoint
  rA <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 1, record_dt = -1)
  rB <- simulate_spatial(T = 5, L = 30, density2 = 0.1, p = 1, seed = 2, record_dt = -1)
  expect_false(identical(rA$n1, rB$n1) && identical(rA$time_sim, rB$time_sim))
})

test_that("simulate_mean_field_stochastic is exactly reproducible for a fixed seed", {
  rA <- simulate_mean_field_stochastic(T = 20, N = 2000, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                                        seed = 11, record_dt = -1)
  rB <- simulate_mean_field_stochastic(T = 20, N = 2000, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                                        seed = 11, record_dt = -1)
  expect_identical(rA$n1[nrow(rA)], rB$n1[nrow(rB)])
  expect_identical(rA$time[nrow(rA)], rB$time[nrow(rB)])
})
