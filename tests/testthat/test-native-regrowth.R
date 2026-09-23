# Post-fire regrowth (state 5, empty_postfire -> native or invader) is
# governed by Lrg_01 (native, from eta_nat) and Lrg_02 (invader, from
# eta_inv) -- see ?xi2lambda. Both channels convert through the exact same
# formula AND the exact same reference rate, L_01 (not L_02 -- see
# ?xi2lambda's Details), so L_01 must stay nonzero here for eta_nat/eta_inv
# to produce any regrowth rate at all; L_01/L_02 are left at their Table 1
# defaults throughout. Ordinary colonization (L_01*(neighbor density),
# L_02*(neighbor density)) is ALSO active from empty-postfire ground, but
# starting every cell at empty-postfire (no natives or invaders anywhere)
# means that channel can never bootstrap a first individual on its own --
# the only way either species can ever appear is through its own Lrg_0x, so
# isolating eta_nat/eta_inv in turn (fire and vegetation-replacement rates
# all off) gives an a priori known outcome: only the species whose eta is
# nonzero can ever appear.

test_that("eta_nat alone drives native regrowth from empty-postfire ground", {
  g <- matrix(FF_STATE["empty_postfire"], nrow = 15, ncol = 15)
  r <- simulate_spatial_from_grid(
    T = 50, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_nat = 0.6, eta_inv = 0,  # native regrowth only
    L_01 = 0.03, L_02 = 0.03, L_12 = 0, L_21 = 0,         # no replacement
    Lig_13 = 0, Lig_23 = 0,                               # no ignition (irrelevant, no fuel)
    seed = 6, record_dt = -1, check_extinction = FALSE
  )
  expect_true(r$n1 > 0)   # native regrew
  expect_equal(r$n2, 0)   # invader never bootstraps: Lrg_02 = 0 and no invader neighbors ever exist
  expect_equal(r$n1 + r$n5, 1)  # everything that isn't still empty_postfire is native
})

test_that("eta_inv alone drives invader regrowth from empty-postfire ground", {
  g <- matrix(FF_STATE["empty_postfire"], nrow = 15, ncol = 15)
  r <- simulate_spatial_from_grid(
    T = 50, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_nat = 0, eta_inv = 0.6,  # invader regrowth only
    L_01 = 0.03, L_02 = 0.03, L_12 = 0, L_21 = 0,
    Lig_13 = 0, Lig_23 = 0,
    seed = 6, record_dt = -1, check_extinction = FALSE
  )
  expect_true(r$n2 > 0)   # invader regrew
  expect_equal(r$n1, 0)   # native never bootstraps: Lrg_01 = 0 and no native neighbors ever exist
  expect_equal(r$n2 + r$n5, 1)
})

test_that("eta_nat = eta_inv = 0 leaves empty-postfire ground untouched", {
  g <- matrix(FF_STATE["empty_postfire"], nrow = 10, ncol = 10)
  r <- simulate_spatial_from_grid(
    T = 50, initial_grid = g,
    xi_nat = 0, xi_inv = 0, eta_nat = 0, eta_inv = 0,   # both regrowth channels off
    L_01 = 0.03, L_02 = 0.03, L_12 = 0, L_21 = 0,
    Lig_13 = 0, Lig_23 = 0,
    seed = 6, record_dt = -1, check_extinction = FALSE
  )
  expect_equal(r$n5, 1)
  expect_equal(r$n1, 0)
  expect_equal(r$n2, 0)
})
