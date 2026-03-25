my_upper <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # Pearson correlation
  r <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
  
  ggplot() +
    geom_tile(aes(0, 0, fill = r)) +
    geom_text(aes(0, 0, label = sprintf("%.2f", r)), size = 9) +
    scale_fill_gradientn(
      name = "Pearson\ncorrelation",
      colours = c(
        "#0073a8", "#66b3e6", "#99d2e7", "white", "#ffb3b3", "#ff6666", "#cc3333"
      ),
      values = scales::rescale(c(-1, -0.5, -0.1, 0, 0.1, 0.5, 1)),
      limits = c(-1, 1)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_test()
}

