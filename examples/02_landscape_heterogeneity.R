# Spatial heterogeneity of the initial landscape
# Reproduces the *type* of panels shown in Figs. 10-13 of the paper:
# initial invader configurations at different values of the heterogeneity
# parameter p (0 = randomly dispersed, 1 = single compact cluster).
#
# Runtime: a few seconds.

library(ForestFireR)

p_values <- c(0, 0.5, 1, 1.5)
L <- 100
density2 <- 0.15

landscapes <- lapply(p_values, function(p) {
  generate_landscape(L = L, density2 = density2, p = p, seed = 123)
})

if (interactive() || !is.null(getOption("device"))) {
  op <- par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
  for (i in seq_along(p_values)) {
    image(landscapes[[i]], col = c("forestgreen", "firebrick"),
          main = paste0("p = ", p_values[i]), axes = FALSE)
  }
  par(op)
}

# Then feed any of these landscapes into a spatial run, e.g. to compare how
# fast the native species goes extinct depending on how segregated the
# invader initially is (paper's main heterogeneity result):
for (i in seq_along(p_values)) {
  r <- simulate_spatial(T = 50, L = L, density2 = density2, p = p_values[i],
                         xi_inv = 0.6, eta_inv = 0.6, seed = 123, record_dt = -1)
  cat(sprintf("p=%.1f  ->  n1=%.3f  n2=%.3f  (after T=%d)\n",
              p_values[i], r$n1, r$n2, 50))
}
