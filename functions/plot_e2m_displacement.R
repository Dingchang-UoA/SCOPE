plot_e2m_displacement  <- function(transform_df){
  ggplot(transform_df) +
    geom_segment(
      aes(
        y = r.e, x = z.e,  
        yend = r.m, xend = z.m, 
        color = status
      ),
      arrow = arrow(length = unit(0.085, "cm")), 
      linewidth = 0.5
    ) +
    geom_point(aes(y = r.m, x = z.m), color = "#FF7F50", size = 1.3) +  
    geom_point(aes(y = r.e, x = z.e), color = "#698b69", size = 1.3) + 
    scale_color_manual(
      name = "Orientation",  
      values = c("inside" = "#8DB6CD", "outside" = "#cd9b9b"),
      labels = c("inside" = "Outwards", "outside" = "Inwards")
    ) +
    labs(
      y = "r",
      x = "z"
    ) +
    theme_article() + 
    coord_equal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "top") +
    facet_wrap(vars(run), ncol = 6)
}