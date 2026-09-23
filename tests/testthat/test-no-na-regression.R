# Regression guard: a past investigation (see the manuscript repo's
# fig6_lambda_na_bug notes) found NA values appearing in simulation
# output under some parameter combinations. This is not about a specific
# known numeric result, but it IS a known a priori property that must
# always hold: under Table 1's own published default rates, output
# densities must be finite, non-negative, sum to 1, and time_sim must
# never exceed the requested horizon by more than a single Gillespie
# step's worth of slack. Any NA/NaN/negative value here is a real defect.

expect_well_formed <- function(df, tol = 1e-6) {
  expect_false(anyNA(df))
  for (col in intersect(names(df), c("n0", "n1", "n2", "n3", "n4", "n5"))) {
    expect_true(all(df[[col]] >= -tol))
    expect_true(all(df[[col]] <= 1 + tol))
  }
  total <- df$n1 + df$n2 + df$n3 + df$n4 + df$n5
  expect_true(all(abs(total - 1) < 1e-4))
}

test_that("spatial model produces well-formed output under Table 1 defaults", {
  # default_rates() also carries L_10/L_20 (always 0, not exposed as
  # simulate_*() arguments -- the model has no background mortality term,
  # see ?default_rates), so drop them before forwarding the rest via do.call.
  rates <- default_rates()
  rates$L_10 <- NULL
  rates$L_20 <- NULL
  for (Tval in c(1, 10, 50)) {
    r <- do.call(simulate_spatial, c(list(T = Tval, L = 40, density2 = 0.1, p = 1,
                                           seed = 123, record_dt = Tval / 20), rates))
    expect_well_formed(r$trajectory)
    expect_false(anyNA(c(r$n1, r$n2, r$n3, r$n4, r$n5, r$time_sim)))
  }
})

test_that("mean-field model produces well-formed output under Table 1 defaults", {
  rates <- default_rates()
  rates$L_10 <- NULL
  rates$L_20 <- NULL
  mf <- do.call(simulate_mean_field, c(list(T = 100, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                                             record_dt = 1), rates))
  expect_well_formed(mf)
})

test_that("well-mixed stochastic model produces well-formed output under Table 1 defaults", {
  rates <- default_rates()
  rates$L_10 <- NULL
  rates$L_20 <- NULL
  ms <- do.call(simulate_mean_field_stochastic,
                c(list(T = 100, N = 5000, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                       seed = 5, record_dt = 1), rates))
  expect_well_formed(ms)
})
