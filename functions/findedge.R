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
