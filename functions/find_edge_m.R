find_edge_m <- function (coordata = coordata, threshold = sqrt(8) ) {
  
  findedge <- function(coordata,center_type=c("central","southeast","northeast","northwest","southwest")){
    
    if (is.null(center_type)) center_type="central"
    
    
    
    if (center_type=="central"){
      
      uniquex=unique(coordata$x)
      
      uniquey=unique(coordata$y)
      
      coordata$edge=FALSE
      
      for (x in uniquex){
        
        
        
        min=min(coordata[coordata$x==x,"y"])
        
        max=max(coordata[coordata$x==x,"y"])
        
        coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
        
        coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
        
      }
      
      
      
      for (y in uniquey){
        
        min=min(coordata[coordata$y==y,"x"])
        
        max=max(coordata[coordata$y==y,"x"])
        
        coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
        
        coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
        
      }
      
      
      
    }else { 
      
      if(center_type=="southeast"){
        
        uniquex=unique(coordata$x)
        
        uniquey=unique(coordata$y)
        
        coordata$edge=FALSE
        
        for (x in uniquex){
          
          
          
          min=min(coordata[coordata$x==x,"y"])
          
          max=max(coordata[coordata$x==x,"y"])
          
          #coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
          
          coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
          
        }
        
        
        
        for (y in uniquey){
          
          min=min(coordata[coordata$y==y,"x"])
          
          max=max(coordata[coordata$y==y,"x"])
          
          #coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
          
          coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
          
        }
        
      } else if(center_type=="northeast"){
        
        uniquex=unique(coordata$x)
        
        uniquey=unique(coordata$y)
        
        coordata$edge=FALSE
        
        for (x in uniquex){
          
          
          
          min=min(coordata[coordata$x==x,"y"])
          
          max=max(coordata[coordata$x==x,"y"])
          
          coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
          
          #coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
          
        }
        
        
        
        for (y in uniquey){
          
          min=min(coordata[coordata$y==y,"x"])
          
          max=max(coordata[coordata$y==y,"x"])
          
          #coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
          
          coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
          
        }
        
      }else if(center_type=="northwest"){
        
        uniquex=unique(coordata$x)
        
        uniquey=unique(coordata$y)
        
        coordata$edge=FALSE
        
        for (x in uniquex){
          
          
          
          min=min(coordata[coordata$x==x,"y"])
          
          max=max(coordata[coordata$x==x,"y"])
          
          coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
          
          #coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
          
        }
        
        
        
        for (y in uniquey){
          
          min=min(coordata[coordata$y==y,"x"])
          
          max=max(coordata[coordata$y==y,"x"])
          
          coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
          
          #coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
          
        }
        
      }else if(center_type=="southwest"){
        
        uniquex=unique(coordata$x)
        
        uniquey=unique(coordata$y)
        
        coordata$edge=FALSE
        
        for (x in uniquex){
          
          
          
          min=min(coordata[coordata$x==x,"y"])
          
          max=max(coordata[coordata$x==x,"y"])
          
          #coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
          
          coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
          
        }
        
        
        
        for (y in uniquey){
          
          min=min(coordata[coordata$y==y,"x"])
          
          max=max(coordata[coordata$y==y,"x"])
          
          coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
          
          #coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
          
        }
        
      }}
    
    
    
    return(coordata)
    
    
    
  }
  
  #Find edge of the individual lenses
  coordata_lt <- list()
  for (imsrun in unique(coordata$run)){
    df_lens <- coordata[coordata$run == imsrun, ]
    df_lens<- findedge(df_lens, "central")
    
    threshold <- threshold
    T_points <- df_lens %>% filter(edge == TRUE) 
    
    calc_distance <- function(x1, y1, x2, y2) {
      sqrt((x1 - x2)^2 + (y1 - y2)^2)}
    
    df_lens <- df_lens %>%
      rowwise() %>%
      mutate(edge = ifelse(
        edge == TRUE,
        {
          # Calculate distances to all other T points
          distances <- sapply(1:nrow(T_points), function(i) {
            if (x == T_points$x[i] && y == T_points$y[i]) {
              # Ignore itself by setting distance to Inf    
              return(Inf) 
            }
            calc_distance(x, y, T_points$x[i], T_points$y[i])
          })
          
          sorted_distances <- sort(distances)
          min1 <- sorted_distances[1]
          min2 <- sorted_distances[2] 
          (min1 <= threshold) & (min2 <= threshold)
        },
        FALSE
      )) %>%
      ungroup()
    
    coordata_lt[[imsrun]] <- df_lens
  }
  
  coordata <- do.call(rbind, coordata_lt)
  return(coordata)
}