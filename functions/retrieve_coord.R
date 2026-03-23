retrieve_coord <- function (imdata = combinedimdata){
  
  coordata<-data.frame(x = imdata@elementMetadata@listData[["x"]],
                       y = imdata@elementMetadata@listData[["y"]],
                       run = run(imdata))
  
 
  return(coordata)
}


