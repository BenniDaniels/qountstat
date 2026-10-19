#' @title CPFISH power
#' @description The basic idea of CPFISH power calculations is to do parametric bootstrapping
#' for each dose/concentration group and to evaluate the proportion of results significantly different from the control.
#' @param contingency.table Matrix with observed data (e.g. survival counts, survival must be in first row)
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param bootstrap.runs Number of bootstrap runs (draw Poisson data n times)
#' @param simulate.p.value Use simulated p-values in Fisher test or not
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param show.progress Show progress for each shift of the probability
#' @param show.results Show results
#' @return Data frame with results from power analysis
#' @examples
#' CPFISH.contingency.table	# example data provided alongside the package
#'
#' # Test CPFISH power
#' CPFISH.power(contingency.table = CPFISH.contingency.table,
#'				control.name = NULL,
#'				alpha = 0.05,
#'				bootstrap.runs = 10,	# Caution: low number of bootstrap runs for testing
#'				simulate.p.value = TRUE,
#'				use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'				show.progress = TRUE,
#'				show.results = TRUE)
#' @export
CPFISH.power = function(contingency.table,				# contingency.table is a matrix with observed data (e.g. survival counts)
						control.name = NULL,			# character string with control group name
						alpha = 0.05,					# significance level
						bootstrap.runs = 200,			# number of bootstrap runs (draw Poisson data n times)
						simulate.p.value = TRUE,		# use simulated p-values or not
						use.fixed.random.seed = NULL,   # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
						show.progress = TRUE,			# show progress for each shift of the probability
						show.results = TRUE) {			# show results

	# check format of input data
	if (!is.numeric(contingency.table) | any(contingency.table < 0, na.rm = TRUE)) {
		stop("Counts must be non-negative numeric values!")
	}

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

	survivors = contingency.table[1, ]  # Survivors have to be in the first row
	introduced = contingency.table[1, ] + contingency.table[2, ]  # Total introduced individuals (survived + died)

	treatment_names = colnames(contingency.table)
	num_treatments = ncol(contingency.table) - 1
	if (is.null(treatment_names)) {
		treatment_names = as.character(1:(num_treatments+1))
	}

	# Fix seed for random numbers if desired (enables to reproduce results)
	if (!is.null(use.fixed.random.seed)) {
		if (!is.numeric(use.fixed.random.seed)) {
			stop("use.fixed.random.seed must be a numeric value or NULL.")
		}
		set.seed(use.fixed.random.seed)
	}

	power.values = c()

	warn.backup = getOption("warn")
	suppressWarnings({

	for (group.index in 2:(num_treatments+1)) {

		if (show.progress) {
			cat("... Calculating power for treatment group '", treatment_names[group.index], "' ...\n", sep = "")
		}

		# setup power
		power.value = NA

		temp.contingency.table = contingency.table

		# Draw N samples from the distribution
		significant.count = 0
		for (i in 1:bootstrap.runs) {
			# Set new data for the currently considered group (other data not changed)
			temp.contingency.table[1, group.index] = stats::rbinom(n = 1, size = introduced[group.index], prob = survivors[group.index] / introduced[group.index])
			temp.contingency.table[2, group.index] = introduced[group.index] - temp.contingency.table[1, group.index]

			# Perform test between original data and simulated data
			results = CPFISH(contingency.table = temp.contingency.table,
							 control.name = NULL,
							 simulate.p.value = T,
							 use.fixed.random.seed = NULL,		# important not to use fixed random seed here
							 show.output = F)

			# Check if the test is significant
			if (results$Results$p.values[group.index-1] < alpha) {
				significant.count = significant.count + 1
			}
		}
		power.value = significant.count / bootstrap.runs * 100
		power.values = c(power.values, power.value)
	}

	})

	# Return the final MDDs
	hypothesis = paste0("H0: ", treatment_names[1], " <-> ", treatment_names[2:length(treatment_names)])
	output = data.frame(Hypothesis = hypothesis, Percent.power = round(power.values, digits = 1))
	if (show.results) {
		cat("\n")
		print(output)
	} else {
		invisible(output)
	}
}
