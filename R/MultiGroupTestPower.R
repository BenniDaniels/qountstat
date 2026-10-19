#' @title Multi-group test power
#' @description Idea: Do parametric bootstrapping for each group and evaluate the proportion of significant results.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param bootstrap.runs Number of bootstrap runs
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param CPCAT.bootstrap.runs Bootstrap runs within CPCAT method
#' @param Dunnett.GLM.zero.treatment.action GLM.Dunnett method to be used for treatments only containing zeros
#' @param show.progress Show progress for each shift of lambda
#' @param show.results Show results
#' @param Dunnett.GLM.family Distribution family to be used for the GLM within the Dunnett.GLM method
#' (either "quasipoisson", "nb" or "poisson"). Only relevant if test = "GLM.Dunnett".
#' @param test Either "CPCAT" or "GLM.Dunnett"
#' @return Data frame with results from power analysis
Multi.group.test.power = function(groups,						# group vector
								  counts,						# vector with count data
								  control.name = NULL,			# character string with control group name
								  alpha = 0.05,					# significance level
								  bootstrap.runs = 200,			# number of bootstrap runs (draw Poisson data n times)
								  use.fixed.random.seed = NULL, # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
								  CPCAT.bootstrap.runs = 200,	# bootstrap runs within CPCAT method
								  Dunnett.GLM.zero.treatment.action = "log(x+1)",	# GLM.Dunnett method to be used for treatments only containing zeros
								  show.progress = TRUE,			# show progress for each shift of lambda
								  show.results = TRUE,			# show results
								  Dunnett.GLM.family = "quasipoisson",	# distribution family used for the GLM within Dunnett.GLM (alternatives: "nb", "poisson")
								  test = "CPCAT") {				# either "CPCAT" or "GLM.Dunnett"

	# check for allowed test type
	if (test != "CPCAT" & test != "GLM.Dunnett") {
		stop("Use either test=\"CPCAT\" or test=\"GLM.Dunnett\"!")
	}
	# check for allowed distribution family
	if (test == "GLM.Dunnett" & (length(Dunnett.GLM.family) != 1 | !is.character(Dunnett.GLM.family))) {
		stop("Parameter Dunnett.GLM.family must be either \"quasipoisson\", \"nb\" or \"poisson\"!")
	}
	if (test == "GLM.Dunnett" & !is.element(Dunnett.GLM.family, c("quasipoisson", "nb", "poisson"))) {
		stop("Parameter Dunnett.GLM.family must be either \"quasipoisson\", \"nb\" or \"poisson\"!")
	}
	# check if there is more than one group
	if (length(unique(groups)) < 2) {
		stop("Too few groups provided!")
	}
	# check if there is count data for each replicate (length of count and group vectors) - groups[i] is one replicate
	if (length(groups) != length(counts)) {
		stop("Lengths of groups and counts don't match!")
	}
	# check format of input data
	if (!is.numeric(counts) | any(counts < 0, na.rm = TRUE)) {		#  | !all(counts == floor(counts))
		stop("Counts must be non-negative numeric values!")
	}

	# Re-structure the input to a data frame
	dat = data.frame(Counts = counts, Groups = groups)

	# Assign new order of levels if control.name was specified
	if (!is.null(control.name)) {
		if (!is.character(control.name)) {
			stop("Specified control must be provided as a character string!")
		}
		if (!is.element(control.name, unique(dat$Groups))) {
			stop("Specified control cannot be found!")
		}

		# Put desired control in the first place
		dat.temp.1 = dat[dat$Groups == control.name,]
		dat.temp.2 = dat[dat$Groups != control.name,]
		dat = rbind(dat.temp.1, dat.temp.2)
	}

	# Convert groups column to a factor, specifying the desired order of levels
	dat$Groups = factor(dat$Groups, levels = unique(dat$Groups))

	# Use treatments vector for convenience
	treatments = levels(dat$Groups)

	# Exit if not enough data left
	if (dim(stats::na.omit(dat))[1] < 2) {
		stop("Too few valid data!")
	}
	if (dim(dat)[1] != dim(stats::na.omit(dat))[1]) {
		# note: these functions return a plain data frame and have no information object, hence the
		# note on excluded rows is issued as a message
		message(dim(dat)[1] - dim(stats::na.omit(dat))[1], " rows with NA values were excluded!")
	}
	dat = stats::na.omit(dat)

	# Get group sizes and mean values
	group.sizes = stats::aggregate(dat$Counts, by=list(dat$Groups), length)$x
	group.means = stats::aggregate(dat$Counts, by=list(dat$Groups), mean)$x

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

	for (group.index in 2:length(levels(dat$Groups))) {

		if (show.progress) {
			cat("... Calculating power for treatment group '", levels(dat$Groups)[group.index], "' ...\n", sep = "")
		}

		# setup power
		power.value = NA

		count.data = dat$Counts
		group.data = dat$Groups

		# Draw 100 samples from the Poisson distribution with the shifted lambda
		significant.count = 0
		for (i in 1:bootstrap.runs) {
			# Set new data for the currently considered group (other data not changed)
			count.data[group.data == levels(dat$Groups)[group.index]] = stats::rpois(group.sizes[group.index], lambda = group.means[group.index])

			# Perform test between original data and simulated data
			if (test == "CPCAT") {
				# do the test
				results = CPCAT(groups = group.data,
								counts = count.data,
								control.name = NULL,
								bootstrap.runs = CPCAT.bootstrap.runs,
								use.fixed.random.seed = NULL, 	# important not to use fixed random seed here
								get.contrasts.and.p.values = F,
								show.output = F)			# do not print results

				# Check if the test is significant
				if (results$Results$p.values[group.index-1] < alpha) {
					significant.count = significant.count + 1
				}
			}
			if (test == "GLM.Dunnett") {
				# do the test
				results = Dunnett.GLM(groups = group.data,
									  counts = count.data,
									  control.name = NULL,
									  zero.treatment.action = Dunnett.GLM.zero.treatment.action,
									  show.output = F,
									  family = Dunnett.GLM.family)

				# Check if the test is significant
				if (results$Results$test$pvalues[group.index-1] < alpha) {
					significant.count = significant.count + 1
				}
			}
		}
		power.value = significant.count / bootstrap.runs * 100
		power.values = c(power.values, power.value)
	}

	})

	# Return the final MDDs
	hypothesis = paste0("H0: ", levels(dat$Groups)[1], " <-> ", levels(dat$Groups)[2:length(levels(dat$Groups))])
	output = data.frame(Hypothesis = hypothesis, Percent.power = round(power.values, digits = 1))
	if (show.results) {
		cat("\n")
		print(output)
	} else {
		invisible(output)
	}
}
