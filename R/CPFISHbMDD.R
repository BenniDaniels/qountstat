#' @title CPFISH bootstrap MDD (bMDD)
#' @description The basic idea of the calculation of bootstrap MDD (bMDD) using the CPCAT approach
#' is to shift the probability of binomial distribution until there is a certain proportion of
#' results significantly different from the control.
#' @param contingency.table Matrix with observed data (e.g. survival counts, survival must be in first row)
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param shift.step Step of shift (negative as a reduction is assumed)
#' @param bootstrap.runs Number of bootstrap runs (draw Poisson data n times)
#' @param power Proportion of bootstrap.runs that return significant differences
#' @param max.iterations Max. number of iterations to not get stuck in the while loop
#' @param simulate.p.value Use simulated p-values in Fisher test or not
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param show.progress Show progress for each shift of the probability
#' @param show.results Show results
#' @return Data frame with results from bMDD analysis
#' @examples
#' CPFISH.contingency.table	# example data provided alongside the package
#'
#' # Test CPFISH bootstrap MDD
#' CPFISH.bMDD(contingency.table = CPFISH.contingency.table,
#'			   control.name = NULL,
#'			   alpha = 0.05,
#'			   shift.step = -0.1,		# Caution: big step size for testing
#'			   bootstrap.runs = 10,		# Caution: low number of bootstrap runs for testing
#'			   power = 0.8,
#'			   max.iterations = 1000,
#'			   simulate.p.value = TRUE,
#'			   use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'			   show.progress = TRUE,
#'			   show.results = TRUE)
#' @export
CPFISH.bMDD = function(contingency.table,			# contingency.table is a matrix with observed data (e.g. survival counts)
					   control.name = NULL,			# character string with control group name
					   alpha = 0.05,				# significance level
					   shift.step = -0.01,			# step of shift (negative as a reduction is assumed)
					   bootstrap.runs = 200,		# number of bootstrap runs (draw Poisson data n times)
					   power = 0.8,					# proportion of bootstrap.runs that return significant differences
					   max.iterations = 1000,		# max number of iterations to not get stuck in the while loop
					   simulate.p.value = TRUE,		# use simulated p-values or not
					   use.fixed.random.seed = NULL,# fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
					   show.progress = TRUE,		# show progress for each shift of lambda
					   show.results = TRUE) {		# show results

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

	# Estimate the initial probability from the data
	initial.prob = survivors[1] / introduced[1]

	# Fix seed for random numbers if desired (enables to reproduce results)
	if (!is.null(use.fixed.random.seed)) {
		if (!is.numeric(use.fixed.random.seed)) {
			stop("use.fixed.random.seed must be a numeric value or NULL.")
		}
		set.seed(use.fixed.random.seed)
	}


	mdds = c()
	mdds.percent = c()

	warn.backup = getOption("warn")
	suppressWarnings({

	for (shift.group.index in 2:(num_treatments+1)) {

		if (show.progress) {
			cat("Calculating bMDD for treatment group '", treatment_names[shift.group.index], "':\n", sep = "")
		}

		# Set current probability to the initial estimate
		current.prob = initial.prob

		# setup MDD and MDD%
		mdd = NA
		mdd.percent = NA

		temp.contingency.table = contingency.table

		iteration = 0
		significant.count = 0

		while (significant.count < bootstrap.runs && iteration < max.iterations) {

			iteration = iteration + 1

			# Draw 100 samples from the Poisson distribution with the shifted lambda
			significant.count = 0
			for (i in 1:bootstrap.runs) {
				# Set new data for the currently considered group (other data not changed)
				temp.contingency.table[1, shift.group.index] = stats::rbinom(n = 1, size = introduced[shift.group.index], prob = current.prob)
				temp.contingency.table[2, shift.group.index] = introduced[shift.group.index] - temp.contingency.table[1, shift.group.index]

				# Perform test between original data and simulated data
				results = CPFISH(contingency.table = temp.contingency.table,
								 control.name = NULL,
								 simulate.p.value = T,
								 use.fixed.random.seed = NULL,		# important not to use fixed random seed here
								 show.output = F)

				# Check if the test is significant
				if (results$Results$p.values[shift.group.index-1] < alpha) {
					significant.count = significant.count + 1
				}
			}

			# Set MDD and MDD%
			if (significant.count / bootstrap.runs < power) {	# reset MDD to NA if power not reached (maybe again)
				mdd = NA
				mdd.percent = NA
			} else {
				if (is.na(mdd)) {
					mdd = abs(current.prob - initial.prob)
					mdd.percent = abs((current.prob - initial.prob) / initial.prob * 100)
				}
			}

			if (iteration == max.iterations) {
				return(cat("Number of iterations exceeded max.iterations!"))
			}

			# Print iteration info
			if (show.progress) {
				cat("Step: ", sprintf("%4d", iteration), " | Prob. diff.: ", sprintf("%.2f", initial.prob-current.prob), " | Test power: ",
					sprintf("%5.1f", round(significant.count / bootstrap.runs * 100, digits=1)),
					" % (", significant.count, "/", bootstrap.runs, ")\n", sep = "")
			}

			# Shift the probability
			current.prob = current.prob + shift.step

			# Do not proceed if next lambda would be < 0 (not suitable for Poisson-distribution)
			if (current.prob < 0) {
				break
			}
		}

		if (show.results) {
			if (current.prob < 0) {
				cat("Shifted probability < 0, procedure stopped.")
			} else {
				if (!is.na(mdd) & !is.na(mdd.percent)) {
					cat(paste0("bMDD (power > ", power, ") = ", mdd, "   (bMDD% = ", sprintf("%.1f", round(mdd.percent, digits=1)), " %)\n\n"))
				}
			}
		}

		mdds = c(mdds, mdd)
		mdds.percent = c(mdds.percent, mdd.percent)
	}
	})


	# Return the final MDDs
	hypothesis = paste0("H0: ", treatment_names[1], " <-> ", treatment_names[2:length(treatment_names)])
	output = data.frame(Hypothesis = hypothesis, bMDD = mdds, bMDD.percent = mdds.percent)
	if (show.results) {
		print(output)
	} else {
		invisible(output)
	}
}
