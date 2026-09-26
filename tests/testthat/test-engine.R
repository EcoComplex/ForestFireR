test_that("tree and linear engines give the same run for a given seed", {
  set.seed(1)
  g <- matrix(sample(1:2, 40 * 40, TRUE, c(0.3, 0.7)), 40)
  run <- function(engine) {
    simulate_spatial_from_grid(T = 30, initial_grid = g, xi_nat = 0.1, xi_inv = 0.9,
                               eta_nat = 0.0066, eta_inv = 0.99, Lig_23 = 1 / 1600,
                               seed = 7, record_dt = 1, check_extinction = FALSE,
                               engine = engine)
  }
  a <- run("linear"); b <- run("tree")
  cols <- c("n1", "n2", "n3", "n4", "n5")
  expect_identical(a$trajectory[, cols], b$trajectory[, cols])
  expect_equal(a$trajectory$time, b$trajectory$time, tolerance = 1e-8)
})

test_that("simulate_spatial accepts both engines and defaults to tree", {
  expect_identical(eval(formals(simulate_spatial)$engine), c("tree", "linear"))
  r1 <- simulate_spatial(T = 5, L = 30, seed = 3, engine = "tree")
  r2 <- simulate_spatial(T = 5, L = 30, seed = 3, engine = "linear")
  expect_equal(c(r1$n1, r1$n2), c(r2$n1, r2$n2))
  expect_error(simulate_spatial(T = 5, L = 30, seed = 3, engine = "fast"))
})
