# Parameter sweep over (xi_inv, eta_inv): native extinction phase diagram
# Reproduces the *type* of analysis behind Figs. 3-8 of the paper: how the
# final native/invader densities depend on the two positive-feedback
# parameters, fire spread (xi_inv) and post-fire regrowth advantage
# (eta_inv) of the invader.
#
# This example uses the fast deterministic mean-field model so the whole
# grid runs in well under a minute; swap simulate_mean_field() for
# simulate_spatial() (with replicates, averaged) to get the noisier but
# more faithful spatial/stochastic version -- much slower, see README for
# expected runtimes.

library(ForestFireR)

xi_seq  <- seq(0, 0.9, by = 0.15)
eta_seq <- seq(0, 0.9, by = 0.15)
T <- 300

grid <- expand.grid(xi_inv = xi_seq, eta_inv = eta_seq)
grid$n1_final <- NA_real_
grid$n2_final <- NA_real_

for (i in seq_len(nrow(grid))) {
  mf <- simulate_mean_field(T = T, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                             xi_inv = grid$xi_inv[i], eta_inv = grid$eta_inv[i],
                             record_dt = T) # only need the final point
  last <- tail(mf, 1)
  grid$n1_final[i] <- last$n1
  grid$n2_final[i] <- last$n2
}

cat("Native density at T =", T, "across the (xi_inv, eta_inv) grid:\n")
print(round(matrix(grid$n1_final, nrow = length(xi_seq),
                    dimnames = list(xi_inv = xi_seq, eta_inv = eta_seq)), 3))

if (interactive() || !is.null(getOption("device"))) {
  m <- matrix(grid$n1_final, nrow = length(xi_seq))
  image(xi_seq, eta_seq, m, xlab = "xi_inv (fire spread)", ylab = "eta_inv (regrowth)",
        main = "Native density (mean field, final)", col = hcl.colors(20, "Greens", rev = TRUE))
}

# ---- Spatial/stochastic version (slower; uncomment to run) --------------
# reps <- 5
# for (i in seq_len(nrow(grid))) {
#   vals <- replicate(reps, {
#     simulate_spatial(T = T, L = 80, density2 = 0.05, p = 1,
#                       xi_inv = grid$xi_inv[i], eta_inv = grid$eta_inv[i],
#                       seed = NULL, record_dt = -1)$n1
#   })
#   grid$n1_final[i] <- mean(vals)
# }
