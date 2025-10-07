library(inphr)
library(ggfortify)
library(ggtda)

row_label_1 <- patchwork::wrap_elements(
  panel = grid::textGrob(expression(sigma == 0.05), rot = 270)
)
row_label_2 <- patchwork::wrap_elements(
  panel = grid::textGrob(expression(sigma == 0.10), rot = 270)
)

set.seed(1234)
n <- 1000
nrep <- 24
# Sample two collections of point clouds
low_noise <- lapply(1:nrep, function(i) {
  tdaunif::sample_arch_spiral(n, sd = 0.05, arms = 2)
})
high_noise <- lapply(1:nrep, function(i) {
  tdaunif::sample_arch_spiral(n, sd = 0.10, arms = 2)
})

# Visualize a few samples
low_noise_pts_plots <- purrr::map(low_noise[1:5], \(.x) {
  .x |>
    as_tibble() |>
    ggplot(aes(x, y)) +
    geom_point() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_fixed()
})
high_noise_pts_plots <- purrr::map(high_noise[1:5], \(.x) {
  .x |>
    as_tibble() |>
    ggplot(aes(x, y)) +
    geom_point() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_fixed()
})

patchwork::wrap_plots(
  c(
    c(low_noise_pts_plots, list(row_label_1)),
    c(high_noise_pts_plots, list(row_label_2))
  ),
  nrow = 2L
) +
  patchwork::plot_layout(widths = c(rep(1, 5), 0.2))

# Compute persistence diagrams
low_noise_dgms <- lapply(low_noise, function(x) {
  TDA::ripsDiag(x, maxdimension = 1, maxscale = 2)$diagram
})
high_noise_dgms <- lapply(high_noise, function(x) {
  TDA::ripsDiag(x, maxdimension = 1, maxscale = 2)$diagram
})
# Coerce to 'persistence' objects
low_noise_pers <- lapply(low_noise_dgms, phutil::as_persistence)
high_noise_pers <- lapply(high_noise_dgms, phutil::as_persistence)
# Combine into a 'persistence_set' object
ph_set <- phutil::as_persistence_set(c(low_noise_pers, high_noise_pers))

prox <- 1

low_noise_plts <- purrr::map(low_noise_pers[1:5], \(.x) {
  pd <- as_tibble(.x)
  pd$dimension <- factor(pd$dimension, levels = 0:2)
  max_prox <- max(pd$death)
  pd |>
    ggplot() +
    coord_fixed() +
    stat_persistence(
      aes(
        start = birth,
        end = death,
        colour = dimension,
        shape = dimension
      ),
      show.legend = TRUE
    ) +
    geom_abline(slope = 1) +
    labs(x = "Birth", y = "Death", color = "Dimension", shape = "Dimension") +
    lims(x = c(0, max_prox), y = c(0, max_prox)) +
    scale_color_discrete(drop = FALSE) +
    scale_linetype_discrete(drop = FALSE) +
    scale_shape_discrete(drop = FALSE) +
    theme_persist()
})

high_noise_plts <- purrr::map(high_noise_pers[1:5], \(.x) {
  pd <- as_tibble(.x)
  pd$dimension <- factor(pd$dimension, levels = 0:2)
  max_prox <- max(pd$death)
  pd |>
    ggplot() +
    coord_fixed() +
    stat_persistence(
      aes(
        start = birth,
        end = death,
        colour = dimension,
        shape = dimension
      ),
      show.legend = TRUE
    ) +
    geom_abline(slope = 1) +
    labs(x = "Birth", y = "Death", color = "Dimension", shape = "Dimension") +
    lims(x = c(0, max_prox), y = c(0, max_prox)) +
    scale_color_discrete(drop = FALSE) +
    scale_linetype_discrete(drop = FALSE) +
    scale_shape_discrete(drop = FALSE) +
    theme_persist()
})

patchwork::wrap_plots(
  c(
    c(low_noise_plts, list(row_label_1)),
    c(high_noise_plts, list(row_label_2))
  ),
  nrow = 2L
) +
  patchwork::plot_layout(guides = "collect", widths = c(rep(1, 5), 0.2)) &
  theme(legend.position = "bottom") # +
# labs(title = "Persistence Diagrams from Noisy Archimedean Spirals")

# Compute pairwise 2-Wasserstein distances between persistence diagrams (dimension 0)
D0 <- phutil::wasserstein_pairwise_distances(
  ph_set,
  p = 2,
  dimension = 0,
  ncores = 12
)

# Compute pairwise 2-Wasserstein distances between persistence diagrams (dimension 1)
D1 <- phutil::wasserstein_pairwise_distances(
  ph_set,
  p = 2,
  dimension = 1,
  ncores = 12
)

# Test in the space of persistence diagrams (dimension 0)
diag_test_dim0 <- two_sample_diagram_test(
  x = low_noise_pers,
  y = high_noise_pers,
  dimension = 0,
  p = 2,
  ncores = 12
)
saveRDS(diag_test_dim0, "data/diag_test_dim0.rds")

diag_test_dim1 <- two_sample_diagram_test(
  x = low_noise_pers,
  y = high_noise_pers,
  dimension = 1,
  p = 2,
  ncores = 12
)
saveRDS(diag_test_dim1, "data/diag_test_dim1.rds")

func_test_dim0 <- two_sample_functional_test(
  x = low_noise_pers,
  y = high_noise_pers,
  dimension = 0,
  representation = "betti"
)
saveRDS(func_test_dim0, "data/func_test_dim0.rds")
plot(func_test_dim0$iwt)

func_test_dim1 <- two_sample_functional_test(
  x = low_noise_pers,
  y = high_noise_pers,
  dimension = 1,
  representation = "betti"
)
saveRDS(func_test_dim1, "data/func_test_dim1.rds")
plot(func_test_dim1$iwt)
