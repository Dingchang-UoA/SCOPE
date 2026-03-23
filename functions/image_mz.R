image_mz <- function(imdata = combinedimdata, mz, tol, unit = c("ppm", "Da")[1], 
                     intensity_scale = c("log10", "log2", "raw", "sqrt"), 
                     interpolate = FALSE, 
                     transform_fn = NULL,
                     mode = c("raw", "aligned")[1],
                     aligned_coord = NULL,
                     run_levels = NULL){
  
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(viridis)))
  suppressMessages(suppressWarnings(library(HiTMaP)))
  suppressMessages(suppressWarnings(library(Cardinal)))
  suppressMessages(suppressWarnings(library(egg)))
  suppressMessages(suppressWarnings(library(LOMAR)))
  suppressMessages(suppressWarnings(library(ggvoronoi)))
  
  ref <- mz(imdata)
  matched_mz <- match_mz(x = mz, ref = ref, tol = tol, unit = unit)
  
  
  if ( mode == "raw"){
    
    coord <- retrieve_coord(imdata)
    
    if (!is.null(run_levels)) {
      coord$run <- factor(coord$run, levels = run_levels)
    }
    
    
    coord$intensity <- intensity(data = imdata, mz, tol, unit)
    
    # Apply customized transformation and store in a new column
    # Transform intensity
    
    if (!is.null(transform_fn)) {
      coord$intensity_trans <- transform_fn(coord$intensity)
      fill_label <- "Transformed(Intensity)"
    } else {
      intensity_scale <- match.arg(intensity_scale)
      coord$intensity_trans <- switch(
        intensity_scale,
        raw = coord$intensity,
        log10 = log10(coord$intensity),
        log2 = log2(coord$intensity),
        sqrt = sqrt(coord$intensity)
      )
      fill_label <- switch(
        intensity_scale,
        raw = "Normalized Intensity",
        log10 = "log10(Intensity)",
        log2 = "log2(Intensity)",
        sqrt = "sqrt(Intensity)"
      )
    }
    
    fig <- ggplot(coord, aes(x, y)) +
      geom_raster(aes(fill = intensity_trans), interpolate = interpolate) +
      facet_wrap(vars(run)) +
      scale_fill_viridis(option = "D", name = fill_label) +  
      egg::theme_article() +
      ggtitle(paste0("matched m/z: ", as.character(matched_mz)))+
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      coord_equal()
    
    return(fig)
    
  } else {
    
    if (is.null(aligned_coord)){
      stop("Aligned data points are needed to plot.")
    }
    
    coord <-  aligned_coord
    coord$intensity <- intensity(data = imdata, mz = mz, tol = tol, unit = unit)
    outline<- LOMAR::standardize_coordinates(lens_m)
    outline<- as.data.frame(outline$X)
    colnames(outline) <- c("y", "x")
    outline <- outline[, c("x", "y")]
    outline$group <- 1
    
    if (!is.null(transform_fn)) {
      coord$intensity_trans <- transform_fn(coord$intensity)
      fill_label <- "Transformed(Intensity)"
    } else {
      intensity_scale <- match.arg(intensity_scale)
      coord$intensity_trans <- switch(
        intensity_scale,
        raw = coord$intensity,
        log10 = log10(coord$intensity),
        log2 = log2(coord$intensity),
        sqrt = sqrt(coord$intensity)
      )
      fill_label <- switch(
        intensity_scale,
        raw = "Normalized Intensity",
        log10 = "log10(Intensity)",
        log2 = "log2(Intensity)",
        sqrt = "sqrt(Intensity)"
      )
    }
    
    
    fig_align <- ggplot(coord) +
      geom_voronoi(aes(z.final, r.final, fill = intensity_trans), outline = outline) +
      scale_fill_viridis(option = "D", name = fill_label) +
      coord_equal() +
      theme_article() +  
      ggtitle(paste0("matched m/z: ", as.character(matched_mz)))+
      labs(x = "x", y = "y") + 
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      facet_wrap(vars(run))
    
    
    return(fig_align)
    
    
  }
  
}
