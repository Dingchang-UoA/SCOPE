find_neighbors_indices <- function(a, b) {
  left_indices <- which(a <= b)               # Indices where A is no greater than B
  right_indices <- which(a > b)              # Indices where A is greater than B
  
  if (b < min(a)) {  # If B is smaller than all values in A
    left_index <- which.max(a)               # Index of the maximum value in A
    right_index <- which.min(a)              # Index of the minimum value in A
  } else if (b >= max(a)) {  # If B is not greater than all values in A, consider
    # the possible situation that the point is on the side
    left_index <- which.max(a)               # Index of the maximum value in A
    right_index <- which.min(a)              # Index of the minimum value in A
  } else {  # Otherwise, find neighbors normally
    left_index <- ifelse(length(left_indices) > 0, 
                         left_indices[which.max(a[left_indices])], 
                         NA)                 # Index of the largest value < b
    right_index <- ifelse(length(right_indices) > 0, 
                          right_indices[which.min(a[right_indices])], 
                          NA)                # Index of the smallest value > b
  }
  
  return(c(left_index = left_index, right_index = right_index))
}