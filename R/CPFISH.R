#' @title CPFISH
#' @description For quantal data (e.g. survival data, '14 out of 20 animals died') e.g. in the form of
#' a contingency table, CPFISH was proposed by Lehmann et al. (2018). Like CPCAT, CPFISH is based on
#' the Closure Principle (CP), but instead of a bootstrapping approach, a Fisher test is performed for
#' all sub-hypotheses to be analyzed. For details on the structure of the input table, please refer to
#' the dataset 'CPFISH.contingency.table' provided alongside this package.
#' @param contingency.table Matrix with observed data (e.g. survival counts, survival must be in first row)
#' @param control.name Character string with control group name (optional)
#' @param simulate.p.value Use simulated p-values in Fisher test or not
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param show.output Show/hide output
#' @return R object with results and information from CPFISH calculations
#' @references
#' Lehmann, R.; Bachmann, J.; Karaoglan, B.; Lacker, J.; Polleichtner, C.; Ratte, H.; Ratte, M. (2018): An alternative approach to overcome shortcomings with multiple testing of binary data in ecotoxicology. In: Stochastic Environmental Research and Risk Assessment, 2018, 32(1), p. 213-222, https://doi.org/10.1007/s00477-017-1392-1
#' @examples
#' CPFISH.contingency.table	# example data provided alongside the package
#'
#' # Test CPFISH
#' CPFISH(contingency.table = CPFISH.contingency.table,
#'			control.name = NULL,
#'			simulate.p.value = TRUE,
#'			use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'			show.output = TRUE)
#' @export
CPFISH = function(contingency.table,					# contingency.table is a matrix with observed data (e.g. survival counts, survival must be in first row)
				  control.name = NULL,					# character string with control group name
				  simulate.p.value = TRUE,				# use simulated p-values or not
				  use.fixed.random.seed = NULL,         # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
				  show.output = TRUE){					# show/hide output

	# setup information to be stored
	info = data.frame(matrix(nrow = 0, ncol = 1))
	info = rbind(info, "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

	# Assign new order of levels if control.name was specified
	if (!is.null(control.name)) {
		if (!is.character(control.name)) {
			stop("Specified control must be provided as a character string!")
		}
		if (!is.element(control.name, unique(colnames(contingency.table)))) {
			stop("Specified control cannot be found!")
		}

		# Put desired control in the first place
		dat.temp.1 = data.frame(contingency.table[, which(colnames(contingency.table) == control.name)])
		colnames(dat.temp.1) = control.name
		dat.temp.2 = contingency.table[, which(colnames(contingency.table) != control.name)]
		contingency.table = cbind(dat.temp.1, dat.temp.2)
	}

	treatment_names = colnames(contingency.table)
	num_treatments = ncol(contingency.table) - 1
	if (is.null(treatment_names)) {
		treatment_names = as.character(1:(num_treatments+1))
	}

	# Vector to store p-values for each treatment
	pvalues = rep(1, num_treatments)

	# Generate all possible hypotheses
	all_hypotheses = CP.hypotheses(n = length(treatment_names) - 1, treatment.names = treatment_names)
	compact_hypotheses = do.call(rbind, all_hypotheses)
	compact_hypotheses = unique(compact_hypotheses)

	# Matrix to track tested hypotheses and their p-values
	pvalue_flags = matrix(-9999, nrow = nrow(compact_hypotheses), ncol = ncol(compact_hypotheses))

	# Fix seed for random numbers if desired (enables to reproduce results)
	if (!is.null(use.fixed.random.seed)) {
		if (!is.numeric(use.fixed.random.seed)) {
			stop("use.fixed.random.seed must be a numeric value or NULL.")
		}
		set.seed(use.fixed.random.seed)
	}

	for (j in 1:num_treatments) {
		# Retrieve contrast matrix for the current treatment hypothesis
		contrasts = CP.hypotheses(n = num_treatments, treatment.names = treatment_names)[[j]]
		matching_rows = numeric()
		for (i in 1:nrow(contrasts)) {
			matching_rows = c(matching_rows, which(apply(compact_hypotheses, 1, identical, contrasts[i, ])))
		}
		already_tested = which(pvalue_flags[matching_rows, j] != -9999)

		# Remove already tested contrasts from the matrix to be tested
		if (length(already_tested) > 0) {
			contrasts = contrasts[-already_tested, ]
		}

		# Ensure contrasts matrix is still a matrix after subsetting
		if (!is.matrix(contrasts)) {
			contrasts = matrix(contrasts, nrow = 1)
		}

		not_tested = which(pvalue_flags[matching_rows, j] == -9999)
		to_be_tested = matching_rows[not_tested]

		# Calculate p-values for the current contrasts
		pvals = rep(0, nrow(contrasts))
		for (l in 1:nrow(contrasts)) {
			test_data = contingency.table[, c(1, (which(contrasts[l, ] == 1) + 1))]
			if (all(rowSums(test_data) != c(0, 0))) {
				pvals[l] = stats::fisher.test(test_data, alternative = "two.sided", simulate.p.value = simulate.p.value)[[1]]
			} else {
				pvals[l] = 1
			}
		}

		# Update p-value flags and track maximum p-value
		pvalue_flags[to_be_tested, j] = pvals
		if (j > 1) {
			pvals_combined = c(pvals, pvalue_flags[matching_rows[-not_tested], j])
		} else {
			pvals_combined = pvals
		}
		pvalues[j] = max(pvals_combined)

		# Propagate the p-values to the next step
		if (j < (length(treatment_names) - 1)) {
			pvalue_flags[, j + 1] = pvalue_flags[, j]
		}
	}

	n = length(treatment_names)
	significances = rep(NA, n - 1)

	# Assign significance levels based on p-values
	for (j in 1:(n - 1)) {
		if (pvalues[j] < 0.05) {
			if (pvalues[j] < 0.01) {
				if (pvalues[j] < 0.001) {
					significances[j] = "***"
				} else {
					significances[j] = "**"
				}
			} else {
				significances[j] = "*"
			}
		} else {
			significances[j] = "."
		}
	}

	# Get NOEC and LOEC
	NOEC = treatment_names[1]
	LOEC = treatment_names[2]
	for (j in 1:(n - 1)) {
		if (pvalues[j] < 0.05) {
			break
		}
		NOEC = treatment_names[j+1]
		if (j == (n - 1)) {
			LOEC = NA
			break
		}
		LOEC = treatment_names[j+2]
	}
	info = rbind(info, paste0("NOEC: ", NOEC, ", LOEC: ", ifelse(is.na(LOEC), "outside tested dose/concentration", LOEC),
							  ". Assuming that any effects are adverse. Otherwise, NOEC and LOEC should be reconsidered."))

	results = data.frame(Treatment = treatment_names[-1], p.values = pvalues, Signif. = significances)

	# Add information about the use of simulated p-values
	if (simulate.p.value) {
		info = rbind(info, "Simulated p-values used.")
	}

	# Set header for information object
	colnames(info) = "Information and warnings:"

	# Provide output
	if (show.output) {
		print(structure(list('Results' = results, Info=info)), row.names = F, quote = F, right = F)
	} else {
		invisible(structure(list('Results' = results, Info=info)))
	}
}
