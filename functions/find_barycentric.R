find_barycentric <- function (Ax, Ay, Bx, By, Cx, Cy, Px, Py){
  
  # This function aims to find the barycentric coordinates of the points within 
  # the given triangle ABC. The argument are (x,y) for three points (Ax, Ay, Bx, 
  # By, Cx, Cy); Px and Py refer to the point within the shape
  
  
  # Create the matrix for the system of equations
  A_matrix <- matrix(c(Ax, Bx, Cx, Ay, By, Cy, 1, 1, 1), nrow =  3)
  
  # Create the P matrix (coordinates of point P)
  P_matrix <- c(Px, Py, 1)
  
  # Solve the system of equations to get the barycentric coordinates
  lambda_values <- solve(t(A_matrix), P_matrix)
  
  return(list(data.frame(lambda1 = lambda_values[1], lambda2 = lambda_values[2], lambda3 = lambda_values[3])))
}