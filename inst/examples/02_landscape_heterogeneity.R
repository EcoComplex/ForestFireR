# Clustering of the initial landscape
#
# How clustered the invader is at t = 0, set by the parameter p of
# generate_landscape(), and what it does to the invasion. p = 0 scatters the
# invader over the grid, p = 1 packs it into a single compact cluster.
#
# Runtime: a few seconds.

library(ForestFireR)

p_values <- c(0, 0.5, 0.75, 1)
L <- 100
density2 <- 0.15

landscapes <- lapply(p_values, function(p) {
  generate_landscape(L = L, density2 = density2, p = p, seed = 123)
})

op <- par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
for (i in seq_along(p_values)) {
  image(landscapes[[i]], col = c("forestgreen", "firebrick"),
        main = paste0("p = ", p_values[i]), axes = FALSE)
}
par(op)

# Realised clustering: the fraction of an invader cell's four neighbours
# that are also invader (periodic boundaries, as in the simulation).
q22 <- function(g, pattern = 2L) {
  n <- nrow(g); up <- c(n, 1:(n - 1)); dn <- c(2:n, 1)
  nb <- (g[up, ] == pattern) + (g[dn, ] == pattern) +
        (g[, up] == pattern) + (g[, dn] == pattern)
  sum(nb[g == pattern]) / (4 * sum(g == pattern))
}

# Then run each landscape forward and compare how far the invasion gets.
for (i in seq_along(p_values)) {
  r <- simulate_spatial(T = 150, L = L, density2 = density2, p = p_values[i],
                        xi_inv = 0.6, eta_inv = 0.6, seed = 123, record_dt = -1)
  cat(sprintf("p=%.2f  neighbour fraction=%.2f  ->  n1=%.3f  n2=%.3f (T=150)\n",
              p_values[i], q22(landscapes[[i]]), r$n1, r$n2))
}
