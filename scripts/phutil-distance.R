library(phutil)

# Plot one trefoil point cloud and one arch spiral point cloud
oldpar <- par(mfrow = c(1, 2))
plot(
  tdaunif::sample_trefoil(n = 1000, sd = 0.1),
  pch = 19,
  cex = 0.5,
  col = "blue",
  asp = 1,
  main = "Trefoil point cloud",
  xlab = "",
  ylab = ""
)
plot(
  tdaunif::sample_arch_spiral(n = 1000, sd = 0.1, arms = 2),
  pch = 19,
  cex = 0.5,
  col = "red",
  asp = 1,
  main = "Arch spiral point cloud",
  xlab = "",
  ylab = ""
)
par(oldpar)

ph_set <- as_persistence_set(c(trefoils, arch_spirals))
D <- wasserstein_pairwise_distances(ph_set, p = 2, dimension = 0, ncores = 12)
# visualize the distance matrix in a square image plot
image(
  as.matrix(D),
  axes = FALSE,
  main = "2-Wasserstein distance between persistence diagrams (dimension 0)",
  sub = "Image plot of the distance matrix",
  col = viridis::viridis(256),
  asp = 1
)
# visualize the distance matrix using multidimensional scaling (MDS)
P <- cmdscale(D, k = 2)
plot(
  P,
  col = rep(c("red", "blue"), each = 24),
  pch = 19,
  xlab = "MDS dimension 1",
  ylab = "MDS dimension 2",
  main = "2-Wasserstein distance between persistence diagrams (dimension 0)",
  sub = "Multidimensional scaling (MDS) projection"
)
legend(
  "topleft",
  legend = c("trefoil", "arch spiral"),
  col = c("red", "blue"),
  pch = 19
)
# Save the distance matrix as RDS file
saveRDS(D, file = "data/wasserstein_distances.rds")
