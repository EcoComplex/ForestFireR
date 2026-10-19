# Trajectories: mean field, well-mixed stochastic and lattice
#
# The three engines of the same model side by side: the deterministic
# well-mixed limit, the well-mixed stochastic process at N = L*L, and the
# spatially explicit lattice simulation, for one (xi_inv, eta_inv) pair.
# What differs between them is spatial structure and the discreteness of
# fire -- see ?`ForestFireR-model`.
#
# Runtime: a few seconds.

library(ForestFireR)

xi_inv  <- 0.6   # invader fire-spread probability
eta_inv <- 0.7   # invader post-fire regrowth probability
T       <- 150
L       <- 80

mf <- simulate_mean_field(T = T, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                          xi_inv = xi_inv, eta_inv = eta_inv)

ws <- simulate_mean_field_stochastic(T = T, N = L * L,
                                     n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                                     xi_inv = xi_inv, eta_inv = eta_inv,
                                     seed = 1, record_dt = T / 200)

sp <- simulate_spatial(T = T, L = L, density2 = 0.05, p = 1,
                       xi_inv = xi_inv, eta_inv = eta_inv,
                       seed = 1, record_dt = T / 200)

panel <- function(time, n1, n2, n3, main) {
  plot(time, n1, type = "l", col = "forestgreen", ylim = c(0, 1),
       xlab = "time (years)", ylab = "density", main = main)
  lines(time, n2, col = "firebrick")
  lines(time, n3, col = "grey40")
}

op <- par(mfrow = c(1, 3))
panel(mf$time, mf$n1, mf$n2, mf$n3, "Mean field")
panel(ws$time, ws$n1, ws$n2, ws$n3, "Well mixed, stochastic")
panel(sp$trajectory$time, sp$trajectory$n1, sp$trajectory$n2,
      sp$trajectory$n3, "Lattice")
legend("right", c("native", "invader", "fire"),
       col = c("forestgreen", "firebrick", "grey40"), lty = 1, bty = "n")
par(op)

cat(sprintf("Mean field            final: n1=%.3f n2=%.3f\n",
            tail(mf$n1, 1), tail(mf$n2, 1)))
cat(sprintf("Well mixed stochastic final: n1=%.3f n2=%.3f (t=%.1f)\n",
            tail(ws$n1, 1), tail(ws$n2, 1), tail(ws$time, 1)))
cat(sprintf("Lattice               final: n1=%.3f n2=%.3f\n", sp$n1, sp$n2))
