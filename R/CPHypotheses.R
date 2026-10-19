#' @title CP hypotheses
#' @description Function to generate hypotheses for the Closure Principle concept using 0/1 contrast matrices
#' @param n Number of groups
#' @param treatment.names Optional vector with names of treatment groups
#' @return Contrast matrix (all 0/1 combinations)
CP.hypotheses = function(n, treatment.names = NULL) {
	combinations_list = list()

	# Generate all possible combinations of treatments
	for (subset_size in 1:n) {
		combinations_list[[subset_size]] = utils::combn(1:n, subset_size)
	}

	hypothesis_matrix_list = list()

	# Create a hypothesis matrix for each treatment
	for (treatment in 1:n) {
		# Initialize a matrix with zeros
		hypothesis_matrix = matrix(0, ncol = n, nrow = (2^(n - 1)))

		# Set column names for the matrix
		if (!is.null(treatment.names) & length(treatment.names) == n) {
			colnames(hypothesis_matrix) = treatment.names
		} else {
			colnames(hypothesis_matrix) = paste("treatment", 1:n)
		}

		row_index = 1  # Start at the first row

		# Fill the matrix based on the combinations
		for (subset_size in 1:length(combinations_list)) {
			for (col in 1:ncol(combinations_list[[subset_size]])) {
				# Check if the current treatment is in the combination
				if (any(combinations_list[[subset_size]][, col] == treatment)) {
					# Mark the corresponding columns in the matrix
					hypothesis_matrix[row_index, combinations_list[[subset_size]][, col]] = 1
					row_index = row_index + 1  # Move to the next row
				}
			}
		}

		# Store the matrix in the list
		hypothesis_matrix_list[[treatment]] = hypothesis_matrix
	}

	# Name each hypothesis matrix according to its treatment
	names(hypothesis_matrix_list) = paste("H0: mu_0 = mu_", 1:n, sep = "")

	return(hypothesis_matrix_list)
}
