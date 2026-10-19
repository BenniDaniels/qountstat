#' @title Multi-group test bMDD
#' @description Idea: shift lambda of Poisson distribution until there is a certain proportion of significant results
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param shift.step Step of shift (negative as a reduction is assumed)
#' @param bootstrap.runs Number of bootstrap runs
#' @param power Proportion of bootstrap.runs that return significant differences
#' @param max.iterations Max. number of iterations to not get stuck in the while loop
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param CPCAT.bootstrap.runs Bootstrap runs within CPCAT method
#' @param Dunnett.GLM.zero.treatment.action Dunnett.GLM method to be used for treatments only containing zeros
#' @param show.progress Show progress for each shift of lambda
#' @param show.results Show results
#' @param get.effect.and.power Return effect size (percent of control) and power for each step (only for last treatment)
#' @param Dunnett.GLM.family Distribution family to be used for the GLM within the Dunnett.GLM method
#' (either "quasipoisson", "nb" or "poisson"). Only relevant if test = "GLM.Dunnett".
#' @param test Either "CPCAT" or "GLM.Dunnett"
#' @return Data frame with results from bMDD analysis
Multi.group.test.bMDD = function(groups,						# group vector
								 counts,						# vector with count data
								 control.name = NULL,			# character string with control group name
								 alpha = 0.05,					# significance level
								 shift.step = -0.25,			# step of shift (negative as a reduction is assumed)
								 bootstrap.runs = 200,			# number of bootstrap runs (draw Poisson data n times)
								 power = 0.8,					# proportion of bootstrap.runs that return significant differences
								 max.iterations = 1000,			# max number of iterations to not get stuck in the while loop
								 use.fixed.random.seed = NULL,  # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
								 CPCAT.bootstrap.runs = 200,	# bootstrap runs within CPCAT method
								 Dunnett.GLM.zero.treatment.action = "log(x+1)",	# GLM.Dunnett method to be used for treatments only containing zeros
								 show.progress = TRUE,			# show progress for each shift of lambda
								 show.results = TRUE,			# show results
								 get.effect.and.power = FALSE,	# return effect size (% of control) and power for each step (only for last treatment)
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

	# Estimate the initial lambda parameter from the data
	initial.lambda = group.means[1]

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

	shift.group.index.start = 2
	if (get.effect.and.power) {
		shift.group.index.start = length(levels(dat$Groups))	# only consider last treatment
		effect.vec = c()
		power.vec = c()
	}

	for (shift.group.index in shift.group.index.start:length(levels(dat$Groups))) {

		if (show.progress) {
			cat("Calculating bMDD for treatment group '", levels(dat$Groups)[shift.group.index], "':\n", sep = "")
		}

		# Set current lambda to the initial estimate
		current.lambda = initial.lambda

		# setup MDD and MDD%
		mdd = NA
		mdd.percent = NA

		count.data = dat$Counts
		group.data = dat$Groups

		iteration = 0
		significant.count = 0

		while (significant.count < bootstrap.runs && iteration < max.iterations) {

			iteration = iteration + 1

			# Draw 100 samples from the Poisson distribution with the shifted lambda
			significant.count = 0
			for (i in 1:bootstrap.runs) {
				# Set new data for the currently considered group (other data not changed)
				count.data[group.data == levels(dat$Groups)[shift.group.index]] = stats::rpois(group.sizes[shift.group.index], lambda = current.lambda)

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
					if (results$Results$p.values[shift.group.index-1] < alpha) {
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
					if (results$Results$test$pvalues[shift.group.index-1] < alpha) {
						significant.count = significant.count + 1
					}
				}
			}

			# Set MDD and MDD%
			if (significant.count / bootstrap.runs < power) {	# reset MDD to NA if power not reached (maybe again)
				mdd = NA
				mdd.percent = NA
			} else {
				if (is.na(mdd)) {
					mdd = abs(current.lambda - initial.lambda)
					mdd.percent = abs((current.lambda - initial.lambda) / initial.lambda * 100)
				}
			}

			if (iteration == max.iterations) {
				return(cat("Number of iterations exceeded max.iterations!"))
			}

			# Print iteration info
			if (show.progress) {
				cat("Step: ", sprintf("%4d", iteration), " | Mean diff.: ", sprintf("%.2f", initial.lambda-current.lambda), " | Test power: ",
					sprintf("%5.1f", round(significant.count / bootstrap.runs * 100, digits=1)),
					" % (", significant.count, "/", bootstrap.runs, ")\n", sep = "")
			}

			if (get.effect.and.power) {
				effect.vec = c(effect.vec, (initial.lambda - current.lambda) / initial.lambda * 100)
				power.vec = c(power.vec, significant.count / bootstrap.runs * 100)
			}

			# Shift the lambda parameter
			current.lambda = current.lambda + shift.step

			# Do not proceed if next lambda would be < 0 (not suitable for Poisson-distribution)
			if (current.lambda < 0) {
				break
			}
		}

		if (show.results) {
			if (current.lambda < 0) {
				cat("Shifted lambda < 0, procedure stopped.")
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

	if (get.effect.and.power) {
		return(data.frame(Groups = rep.int(x = length(levels(dat$Groups)), times = length(effect.vec)),
						  Effect = effect.vec, Power = power.vec))
	}

	# Return the final MDDs
	hypothesis = paste0("H0: ", levels(dat$Groups)[1], " <-> ", levels(dat$Groups)[2:length(levels(dat$Groups))])
	output = data.frame(Hypothesis = hypothesis, bMDD = mdds, bMDD.percent = mdds.percent)
	if (show.results) {
		print(output)
	} else {
		invisible(output)
	}
}
