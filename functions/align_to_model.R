align_to_model  <- function (coordata = coordata, model = lens_m ) {
  
  transform_lt <- list()
  lens_e_t_lt <- list() #List of transformed data points
  
  
  for (imsrun in unique(coordata$run)){
    message("This is ", which(unique(coordata$run) == imsrun), " out of ", length(unique(coordata$run)), " lenses")
    lens_e_full <- coordata %>% filter(run== imsrun )
    lens_e_geom <- coordata %>% filter(run== imsrun, edge ==T )
    transform<- icp_m(X=model, Y= lens_e_geom[,c("y", "x")], Y_full = lens_e_full[,c("y", "x")] ,tol = 1e-11, scale = T)
    
    # "mp" refers to outline points from model; 
    # "dp" refers the outline points from observation
    mp <- as.data.frame(transform$X)
    colnames(mp) <- c("r","z")
    mp$Source <- "Model"
    mp<-mp %>% mutate(Distance = sqrt(r^2+z^2))
    
    dp <- as.data.frame(transform$Y)
    colnames(dp) <- c("r","z")
    dp$Source <- "Data"
    dp<-dp %>% mutate(Distance = sqrt(r^2+z^2))
    
    mp$ID <- seq_len(nrow(mp))
    dp$ID <- seq_len(nrow(dp))
    
    df_pair <- bind_rows(mp, dp)
    
    df_comparison <- merge(mp,dp, by = "ID", all =T)
    
    colnames(df_comparison) <- c("ID", "r.m", "z.m", "Source_m", "Distance.m", "r.e", "z.e", "Source.e", "Distance.e")
    
    #Check if there are several observed points matched to the same model point
    df_comparison <- remove_duplicates(df_comparison, r.e = "r.e", z.e = "z.e", r.m = "r.m", z.m = "z.m")
    
    
    df_comparison <- df_comparison %>%
      mutate(
        status = ifelse(Distance.m > Distance.e, "inside", "outside"),
        run = imsrun
      )
    df_comparison <-df_comparison %>% mutate(theta_geom = atan2(z.e,r.e))  
    
    df_comparison <- df_comparison %>% mutate(delta.r = r.m -r.e,
                                              delta.z = z.m - z.e)
    transform_lt[[imsrun]] <- df_comparison
    
    #Compute movement 
    lens_e_t <- as.data.frame(transform$Y_full_t) # transformed data points
    colnames(lens_e_t )<- c("r","z")
    lens_e_t <- lens_e_t %>% mutate(theta = atan2(z,r))
    
    neighbor_indices <- t(sapply(lens_e_t$theta, function(b) find_neighbors_indices(df_comparison$theta_geom, b)))
    
    
    lens_e_t<-cbind(lens_e_t, neighbor_indices)
    
    lambda_lt <- list()
    
    for (i in 1:nrow(lens_e_t)){
      Ax <- 0
      Ay <- 0 
      Bx <- df_comparison$r.e[lens_e_t$left_index][i]
      By <- df_comparison$z.e[lens_e_t$left_index][i]
      Cx <- df_comparison$r.e[lens_e_t$right_index][i]
      Cy <- df_comparison$z.e[lens_e_t$right_index][i]
      
      Px <- lens_e_t$r[i]
      Py <- lens_e_t$z[i]
      
      lambda_lt[i] <- try(find_barycentric(Ax, Ay, Bx, By, Cx, Cy, Px, Py))
    }
    
    lambdas<- do.call(rbind, lambda_lt)
    lambdas <- data.frame(sapply(lambdas, as.numeric))
    
    lens_e_t <- cbind(lens_e_t, lambdas)
    
    lens_e_t <- lens_e_t %>%  mutate(displacement.r = lambda2*(df_comparison$r.m[left_index] - df_comparison$r.e[left_index])+
                                       lambda3*(df_comparison$r.m[right_index] - df_comparison$r.e[right_index]))
    
    lens_e_t <- lens_e_t %>%  mutate(displacement.z = lambda2*(df_comparison$z.m[left_index] - df_comparison$z.e[left_index])+
                                       lambda3*(df_comparison$z.m[right_index] - df_comparison$z.e[right_index]))    
    
    lens_e_t <- lens_e_t %>%  mutate(r.final = r + displacement.r,
                                     z.final = z +displacement.z)
    lens_e_t$run <- imsrun
    
    lens_e_t_lt[[imsrun]]<- lens_e_t
    
  }
  
  transform_df <- do.call(rbind, transform_lt)
  lens_e_t_df <- do.call(rbind, lens_e_t_lt)
  
  lens_e_t_df$run <-  factor(lens_e_t_df$run, levels = unique(coordata$run))
  transform_df$run <- factor(transform_df$run, levels = unique(coordata$run))
 
  return(list(
    transform_df = transform_df,
    lens_e_t_df = lens_e_t_df
  ))
}

