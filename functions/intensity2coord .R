intensity2coord <- function(matrix_data) {
  matrix_data <- as.matrix(matrix_data)
  # Get the row and column indices using expand.grid
  coordinates <- expand.grid(x = 1:nrow(matrix_data), y = 1:ncol(matrix_data))
  
  # Map intensity values cell by cell
  coordinates$intensity <- matrix_data[cbind(coordinates$x, coordinates$y)]
  
  # Remove rows with NA intensity
  coordinates <- coordinates[!is.na(coordinates$intensity), ]
  
  return(coordinates)
}