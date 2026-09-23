# xi2lambda() is a pure R function with an exact closed-form result, so
# every case here is checked against a value worked out by hand from the
# formula in ?xi2lambda: lambda = xi/(1-xi) * lambda_ref.

test_that("xi2lambda gives exact values at known points", {
  expect_equal(xi2lambda(0, 5), 0)                 # no spread probability -> zero rate
  expect_equal(xi2lambda(0.5, 5), 5)                # xi/(1-xi) = 1 at xi = 0.5
  expect_equal(xi2lambda(0.6, 1), 0.6 / 0.4)        # = 1.5, worked out by hand
  expect_equal(xi2lambda(0.9, 10), 90)              # xi/(1-xi) = 9 at xi = 0.9
})

test_that("xi2lambda is monotonically increasing in xi (for lambda_ref > 0)", {
  xs <- seq(0, 0.99, by = 0.01)
  lams <- xi2lambda(xs, 1)
  expect_true(all(diff(lams) > 0))
})

test_that("xi2lambda rejects xi outside [0, 1)", {
  expect_error(xi2lambda(-0.01, 5))
  expect_error(xi2lambda(1, 5))
  expect_error(xi2lambda(1.5, 5))
})

test_that("eta2lambda is an alias of xi2lambda", {
  expect_equal(eta2lambda(0.6, 1), xi2lambda(0.6, 1))
})

test_that("default_rates() gives Table 1's baseline values", {
  r <- default_rates()
  expect_equal(r$L_01, 0.03)
  expect_equal(r$L_02, 0.03)
  expect_equal(r$L_12, 0.005)
  expect_equal(r$L_21, 0.01)
  expect_equal(r$Lig_13, 0)     # native is not fire-adapted (condition iv)
  expect_equal(r$Lig_23, 1e-4)
  expect_equal(r$L_10, 0)       # no background mortality in the paper's model
  expect_equal(r$L_20, 0)
})
