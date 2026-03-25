my_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      contour = TRUE
    ) +
    scale_fill_gradientn(
      colours = c(
        "#FFFFFF", "#17BECF", "#1F77B4", "#2CA02C",
        "#FFFF00", "#FF7F0E", "#D62728", "#8B0000"
      )
    ) +
    theme_test() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 13,
        color = "black"
      ),
      axis.text.y = element_text(
        size  = 13,
        color = "black"
      )
    )
}