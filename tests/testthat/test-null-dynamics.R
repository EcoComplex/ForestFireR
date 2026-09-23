# With every reaction rate set to zero, nothing in the model can change:
# no colonization, no competitive replacement, no ignition, no fire
# spread, no regrowth. The grid must therefore come out of the engine
# byte-for-byte identical to how it went in, and simulated time should
# not advance (the Gillespie total propensity is exactly zero, so the
# engine has no event to draw and stops immediately). This is a known,
# exact a priori result, not a statistical one.

test_that("all-zero rates freeze the grid exactly (spatial)", {
  g0 <- generate_landscape(L = 15, density2 = 0.1, p = 1, seed = 1)
  r <- simulate_spatial_from_grid(
    T = 5, initial_grid = g0,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,
    L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0, Lig_13 = 0, Lig_23 = 0,
    seed = 1, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_identical(r$initial_grid, r$final_grid)
  expect_identical(r$initial_grid, g0)
  expect_equal(r$n3, 0)
  expect_equal(r$time_sim, 0)
})

test_that("all-zero rates freeze the grid exactly regardless of initial composition", {
  # a mixed grid (native, invader, both empty sub-states, no fire) should
  # freeze exactly the same way
  g <- matrix(c(1L, 2L, 4L, 5L), nrow = 10, ncol = 10)
  r <- simulate_spatial_from_grid(
    T = 100, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,
    L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0, Lig_13 = 0, Lig_23 = 0,
    seed = 2, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_identical(r$final_grid, g)
  expect_equal(r$time_sim, 0)
})
