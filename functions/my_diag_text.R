my_diag_text <- function(data, mapping, ...) {
  var_name <- as_label(mapping$x)
  
  ggplot() +
    geom_text(
      aes(0, 0, label = var_name),
      size = 5,
      fontface = "bold"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void()
}