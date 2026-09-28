# Random number generator (v0.8.0): std::mt19937_64 by default, the C
# library rand() with options(ForestFireR.rng = "legacy").

with_rng <- function(kind, code) {
  old <- options(ForestFireR.rng = kind); on.exit(options(old)); force(code)
}

test_that("the default generator is reproducible for a given seed", {
  a <- simulate_spatial(T = 10, L = 30, seed = 11)
  b <- simulate_spatial(T = 10, L = 30, seed = 11)
  expect_identical(c(a$n1, a$n2, a$time_sim), c(b$n1, b$n2, b$time_sim))
})

test_that("mt19937 and legacy generators give different runs for the same seed", {
  g1 <- with_rng("mt19937", generate_landscape(L = 30, density2 = 0.2, p = 0, seed = 5))
  g2 <- with_rng("legacy",  generate_landscape(L = 30, density2 = 0.2, p = 0, seed = 5))
  expect_false(identical(g1, g2))
  expect_identical(sum(g1 == 2L), sum(g2 == 2L))   # same density either way
})

test_that("mt19937 landscapes are the same on every platform", {
  # Reference values computed on Linux (glibc, x86-64); std::mt19937_64's
  # output is fixed by the C++ standard, so any platform must match.
  g <- with_rng("mt19937", generate_landscape(L = 30, density2 = 0.2, p = 0.8, seed = 42))
  expect_identical(c(sum(g == 2L), sum(which(g == 2L))), c(180L, 89021L))
  g <- with_rng("mt19937", generate_landscape(L = 30, density2 = 0.2, p = 0, seed = 42))
  expect_identical(c(sum(g == 2L), sum(which(g == 2L))), c(180L, 80562L))
})

test_that("an invalid ForestFireR.rng option is rejected", {
  expect_error(with_rng("xorshift", simulate_spatial(T = 1, L = 10, seed = 1)))
})
