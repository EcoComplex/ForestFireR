# Fire can only start via the ignition reaction (rate Lig_23 for the
# invader, Lig_13 for the native -- currently always 0 in the paper's own
# parametrization, condition iv). With both ignition rates at exactly
# zero, fire (state 3) and post-fire empty ground (state 5) must never
# appear anywhere in the trajectory, no matter how large the fire-spread
# rates are set (spread only matters once something is already burning).

test_that("fire and post-fire empty ground never appear when ignition is off", {
  g <- generate_landscape(L = 20, density2 = 0.5, p = 1, seed = 3)
  r <- simulate_spatial_from_grid(
    T = 50, initial_grid = g,
    xi_nat = 0.5, xi_inv = 0.9, eta_inv = 0.6,     # fire spread rates set high...
    L_01 = 0.03, L_02 = 0.03, L_12 = 0.01, L_21 = 0.02,
    Lig_13 = 0, Lig_23 = 0,                        # ...but ignition is off
    seed = 3, record_dt = 1, check_extinction = FALSE
  )
  expect_true(all(r$trajectory$n3 == 0))
  expect_true(all(r$trajectory$n5 == 0))
  expect_equal(r$n3, 0)
  expect_equal(r$n5, 0)
})

test_that("Table 1's default native ignition rate (Lig_13 = 0) alone keeps fire native-free", {
  # even with the invader able to ignite, native-only ignition must stay off:
  # start from an ALL-NATIVE grid (no invader anywhere) with Lig_23 active --
  # ignition should never occur, because Lig_13 (native) is 0 and there is no
  # invader present to trigger Lig_23.
  g <- matrix(1L, nrow = 15, ncol = 15)
  r <- simulate_spatial_from_grid(
    T = 50, initial_grid = g,
    xi_nat = 0.9, xi_inv = 0.9, eta_inv = 0.6,
    L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0,
    Lig_13 = 0, Lig_23 = 1e-4,
    seed = 4, record_dt = 5, check_extinction = FALSE
  )
  expect_true(all(r$trajectory$n3 == 0))
  expect_equal(r$n3, 0)
})
