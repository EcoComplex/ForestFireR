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

test_that("the mt19937 random stream is the same on every platform", {
  # std::mt19937_64's output for a given seed is fixed by the C++ standard.
  # Reference values computed on Linux (glibc, x86-64). Uniforms are
  # k / 2^53 with integer k, so k is compared exactly.
  r <- with_rng("mt19937", { ForestFireR:::.ffr_apply_rng(); ForestFireR:::rng_draws_cpp(3L, 42L) })
  expect_identical(r$u * 2^53, c(6801836353641660, 5755883094484128, 6774721691634639))
  expect_identical(r$k, c(136272683L, 903268966L, 94068311L))
})

test_that("mt19937 landscapes are reproducible for a given seed", {
  # Landscapes (and simulations) are reproducible on a given machine, but
  # not bit-identical across platforms: the model computes in `long double`,
  # which is 80-bit on x86-64 and 64-bit on Apple Silicon, so comparisons of
  # nearly equal quantities can go different ways even with the same random
  # stream.
  g1 <- with_rng("mt19937", generate_landscape(L = 30, density2 = 0.2, p = 0.8, seed = 42))
  g2 <- with_rng("mt19937", generate_landscape(L = 30, density2 = 0.2, p = 0.8, seed = 42))
  expect_identical(g1, g2)
})

test_that("an invalid ForestFireR.rng option is rejected", {
  expect_error(with_rng("xorshift", simulate_spatial(T = 1, L = 10, seed = 1)))
})
