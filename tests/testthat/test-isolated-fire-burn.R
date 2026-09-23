# With fire spread disabled (xi_nat = xi_inv = 0) and no new ignition
# (Lig_13 = Lig_23 = 0), a fire that is already burning cannot spread to
# any neighbor and cannot be replaced by a new outbreak. Starting from a
# grid with exactly one site on fire, exactly that one site must burn
# out to post-fire empty ground (state 5) and nothing else may change:
# an exact, deterministic, a priori known outcome.

test_that("a single isolated fire burns out to exactly one empty-postfire site", {
  g <- matrix(1L, nrow = 10, ncol = 10) # all native
  g[5, 5] <- 3L                          # one site already on fire
  r <- simulate_spatial_from_grid(
    T = 10, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,   # no spread
    L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0,
    Lig_13 = 0, Lig_23 = 0,                # no new ignition
    L_30 = 1e6,
    seed = 4, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_equal(sum(r$final_grid == 5), 1)   # exactly the one burned site
  expect_equal(sum(r$final_grid == 3), 0)   # fire itself is gone
  expect_equal(sum(r$final_grid == 1), 99)  # everything else untouched
  expect_equal(r$final_grid[5, 5], 5L)      # it's the SAME site that burned
})

test_that("several isolated fires each burn out independently with no spread", {
  g <- matrix(1L, nrow = 12, ncol = 12)
  fire_sites <- rbind(c(2, 2), c(2, 10), c(10, 2), c(10, 10), c(6, 6))
  for (i in seq_len(nrow(fire_sites))) g[fire_sites[i, 1], fire_sites[i, 2]] <- 3L

  r <- simulate_spatial_from_grid(
    T = 10, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,
    L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0,
    Lig_13 = 0, Lig_23 = 0,
    L_30 = 1e6,
    seed = 5, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_equal(sum(r$final_grid == 5), nrow(fire_sites))
  expect_equal(sum(r$final_grid == 3), 0)
  for (i in seq_len(nrow(fire_sites))) {
    expect_equal(r$final_grid[fire_sites[i, 1], fire_sites[i, 2]], 5L)
  }
})
