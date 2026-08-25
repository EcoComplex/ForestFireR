# Trajectories: mean-field vs. spatial model
# Reproduces the *type* of comparison shown in Figs. 2 and 9 of the paper:
# native/invader/fire densities over time, deterministic mean-field vs.
# the spatially-explicit stochastic model, for a given (xi_inv, eta_inv).
#
# Runtime: a few seconds at these (reduced, interactive) settings. Scale
# L up to 200 and T up to ~200 for paper-resolution runs (background job).

library(ForestFireR)

xi_inv  <- 0.6   # invader fire-spread probability
eta_inv <- 0.7   # invader post-fire regrowth probability
T       <- 100
L       <- 80    # bump to 200 for paper resolution (slower, see README)

mf <- simulate_mean_field(T = T, n1_0 = 0.9, n2_0 = 0.05, n4_0 = 0.05,
                           xi_inv = xi_inv, eta_inv = eta_inv)

sp <- simulate_spatial(T = T, L = L, density2 = 0.05, p = 1,
                        xi_inv = xi_inv, eta_inv = eta_inv,
                        seed = 1, record_dt = T / 200)

if (interactive() || !is.null(getOption("device"))) {
  op <- par(mfrow = c(1, 2))
  plot(mf$time, mf$n1, type = "l", col = "forestgreen", ylim = c(0, 1),
       xlab = "time", ylab = "density", main = "Mean field")
  lines(mf$time, mf$n2, col = "firebrick")
  lines(mf$time, mf$n3, col = "grey40")
  legend("right", c("native", "invader", "fire"),
         col = c("forestgreen", "firebrick", "grey40"), lty = 1, bty = "n")

  plot(sp$trajectory$time, sp$trajectory$n1, type = "l", col = "forestgreen",
       ylim = c(0, 1), xlab = "time", ylab = "density", main = "Spatial (stochastic)")
  lines(sp$trajectory$time, sp$trajectory$n2, col = "firebrick")
  lines(sp$trajectory$time, sp$trajectory$n3, col = "grey40")
  par(op)
}

cat(sprintf("Mean field   final: n1=%.3f n2=%.3f n3=%.4f\n",
            tail(mf$n1,1), tail(mf$n2,1), tail(mf$n3,1)))
cat(sprintf("Spatial      final: n1=%.3f n2=%.3f n3=%.4f\n", sp$n1, sp$n2, sp$n3))
