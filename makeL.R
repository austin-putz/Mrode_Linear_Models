#------------------------------------------------------------------------------#
# makeL function
#------------------------------------------------------------------------------#

# Author: Austin Putz
# Created: Sun Feb 23, 2025
# Modified: Sun Feb 23, 2025
# Contact: putz [dot] austin [at] gmail.com
# License: GPLv2

# function
makeL <- function(ped) {
  
  # Ensure the Matrix package is available
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but not installed. Please install it using install.packages('Matrix').")
  }
  
  # ped: a data frame with columns: Animal, Sire, Dam.
  #      Missing parents are coded as 0.
  #      Pedigree is assumed to be sorted oldest (top) to youngest (bottom)
  #
  # This function computes the lower triangular matrix L (and inbreeding coefficients)
  # for the pedigree. Internally, a dummy index (1) is used for missing parents.
  
  # Check pedigree numbering: Animal IDs should run 1 to n_actual
  n_actual <- nrow(ped)
  if (ped$Animal[1] != 1 || tail(ped$Animal, 1) != n_actual) {
    stop("Pedigree is not numbered 1 to n")
  }
  
  # Define total number of indices including a dummy (index 1) for missing parents
  n_total <- n_actual + 1
  
  # Create an adjusted pedigree where 0's become 1 (dummy index).
  # Only adjust Sire and Dam; optionally, adjust Animal IDs.
  ped_adj <- ped
  ped_adj$Sire   <- ped_adj$Sire + 1
  ped_adj$Dam    <- ped_adj$Dam + 1
  ped_adj$Animal <- ped_adj$Animal + 1  # now actual animals are 2:(n_total)
  
  # Initialize the L matrix and a vector to hold the cumulative sums of squares.
  L <- matrix(0, n_total, n_total)
  F_vals <- numeric(n_total)
  
  # Loop over columns (j corresponds to each animal, adjusted)
  # We start at 2 because index 1 is reserved for missing parents.
  for (j in 2:n_total) {
    # For each column, process from the diagonal (i == j) downwards.
    for (i in j:n_total) {
      if (i == j) {
        # Diagonal: calculate using the parents' inbreeding values.
        L[i, j] <- sqrt(1 - 0.25 * (F_vals[ped_adj$Sire[i-1]] + F_vals[ped_adj$Dam[i-1]]))
        # Update the cumulative sum of squares for this animal
        F_vals[j] <- sum(L[i, 1:j]^2)
      } else {
        # Off-diagonal: use the corresponding sire and dam entries.
        L[i, j] <- 0.5 * L[ped_adj$Sire[i-1], j] + 0.5 * L[ped_adj$Dam[i-1], j]
      }
    }
  }
  
  # Remove the dummy row and column before returning.
  L_final <- Matrix::Matrix(L[-1, -1])
  
  # remove first element (dummy) and subtract 1 to get F values
  F_final <- F_vals[-1] - 1
  
  # return list
  return(list(L = L_final, F_vals = F_final))
}



# Example to test with
# 
# # Example usage with your pedigree:
# ped <- data.frame(
#   Animal = c(1,2,3,4,5,6,7,8),
#   Sire   = c(0,0,0,1,1,4,4,6),
#   Dam    = c(0,0,0,2,3,5,5,7)
# )
# 
# result <- makeL(ped)
# print(result$L)
# print(result$F_vals)




