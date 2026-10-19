# Parameter sweep over (xi_inv, eta_inv): native extinction phase diagram
#
# The invader's two fire-mediated advantages are its fire-spread probability
# (xi_inv) and its post-fire regrowth probability (eta_inv); together they
# set the strength of the invasion-fire feedback (?`ForestFireR-model`).
# This sweep maps the final native density over that plane.
#
# The grid is run with the deterministic mean-field model, which takes
# milliseconds per point; the commented block at the end does the same with
# the lattice model, averaged over replicates.

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
                            record_dt = T) # only the final point is needed
  last <- tail(mf, 1)
  grid$n1_final[i] <- last$n1
  grid$n2_final[i] <- last$n2
}

cat("Native density at T =", T, "across the (xi_inv, eta_inv) grid:\n")
print(round(matrix(grid$n1_final, nrow = length(xi_seq),
                   dimnames = list(xi_inv = xi_seq, eta_inv = eta_seq)), 3))

image(xi_seq, eta_seq, matrix(grid$n1_final, nrow = length(xi_seq)),
      xlab = "xi_inv (fire spread)", ylab = "eta_inv (post-fire regrowth)",
      main = "Final native density (mean field)",
      col = hcl.colors(20, "Greens", rev = TRUE))

# ---- Lattice version (slower; uncomment to run) -------------------------
# reps <- 5
# for (i in seq_len(nrow(grid))) {
#   vals <- replicate(reps, {
#     simulate_spatial(T = T, L = 80, density2 = 0.05, p = 1,
#                      xi_inv = grid$xi_inv[i], eta_inv = grid$eta_inv[i],
#                      seed = NULL, record_dt = -1)$n1
#   })
#   grid$n1_final[i] <- mean(vals)
# }
