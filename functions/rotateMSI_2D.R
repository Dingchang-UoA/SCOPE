rotateMSI_2D <- function (imdata, rotation_degree = 0) 
{
  rotatetmp <- as.data.frame(imdata@elementMetadata@listData)
  # rotatetmp <- as.data.frame(imdata@elementMetadata@coord@listData )
  
  
  rotatenew <- HiTMaP:::affine_grid(rotatetmp[, c("x", "y")], rotate = rotation_degree, 
    grid = T)
  rotatenew$x <- rotatenew$x - min(rotatenew$x) + 1
  rotatenew$y <- rotatenew$y - min(rotatenew$y) + 1
  rotatenew$z = rotatetmp$z
  
  rotatetmp$x <- rotatenew$x
  rotatetmp$y <- rotatenew$y 
  

  
  imdata@elementMetadata@listData  <- as.list(rotatetmp)
  names(imdata@elementMetadata@listData) <- c("x", "y", "z", "run", 
                                              "3DPositionX", 
                                              "3DPositionY",
                                              "3DPositionZ")
  return(imdata)
}
