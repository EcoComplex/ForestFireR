# With only one species present and no fire, ordinary colonization
# (reaction 1) is the only active process, and it can only ever fill
# empty ground with that same species. With no competing sink term the
# system has a single deterministic absorbing state: every empty site
# eventually becomes native. This isolates reaction 1 from everything
# else in the model (competition, fire, regrowth all set to zero).

test_that("native-only colonization deterministically fills the grid", {
  layers <- list(list(background = 1, pattern = 4, density = 0.3, p = 1))
  g <- generate_landscape_layers(L = 15, fill_state = 1, layers = layers, seed = 2)
  expect_true(any(g == 4)) # sanity: the initial grid does have empty sites to fill

  r <- simulate_spatial_from_grid(
    T = 1000, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,          # no fire at all
    L_01 = 0.03, L_02 = 0.03, L_12 = 0, L_21 = 0, Lig_13 = 0, Lig_23 = 0,
    seed = 2, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_equal(r$n1, 1, tolerance = 1e-9)   # fully colonized by the native
  expect_equal(r$n0, 0, tolerance = 1e-9)   # no empty ground left
  expect_true(all(r$final_grid == 1))       # every site, not just the aggregate density
})

test_that("invader-only colonization deterministically fills the grid the same way", {
  layers <- list(list(background = 2, pattern = 4, density = 0.3, p = 1))
  g <- generate_landscape_layers(L = 15, fill_state = 2, layers = layers, seed = 3)

  r <- simulate_spatial_from_grid(
    T = 1000, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_inv = 0,
    L_01 = 0.03, L_02 = 0.03, L_12 = 0, L_21 = 0, Lig_13 = 0, Lig_23 = 0,
    seed = 3, record_dt = -1, record_grid = TRUE, check_extinction = FALSE
  )
  expect_equal(r$n2, 1, tolerance = 1e-9)
  expect_equal(r$n0, 0, tolerance = 1e-9)
  expect_true(all(r$final_grid == 2))
})
