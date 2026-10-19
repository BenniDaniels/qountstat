#' @title CPCAT
#' @description When conducting statistical tests with multiple treatments, such as a control group and
#' increasing concentrations of a test substance, ANOVA and parametric post-hoc tests (e.g. Dunnett's test)
#' are commonly used. However, these tests require the assumptions of homogeneous variances and normally
#' distributed data. For count data (e.g. counts of animals), these assumptions are typically violated,
#' as the data are usually Poisson-distributed. Additionally, multiple testing using post-hoc tests can
#' lead to alpha-inflation. To address these issues, CPCAT was proposed by Lehmann et al. (2016). CPCAT has
#' two components. The first is the Closure Principle (CP) developed by Bretz et al. (2010), which aims
#' to eliminate alpha-inflation. CP applies a stepwise approach to identify at which concentration effects
#' begin to occur. The second part of CPCAT is the actual significance test, CAT (Computational Approach
#' Test; introduced by Chang et al., 2010), which uses a test based on the Poisson distribution rather
#' than a parametric test based on normal distribution assumptions. For details on the structure of the
#' input data, please refer to the dataset 'Daphnia.counts' provided alongside this package.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param bootstrap.runs Number of bootstrap runs
#' @param hampel.threshold Threshold for Hampel identifier (measure for over-/underdispersion). If NULL
#' (default) the Hampel identifier is not evaluated. The Hampel identifier compares the difference between
#' variance and mean of each treatment group with a fixed threshold. Since the sampling variability of
#' this difference grows proportionally to the mean while the threshold does not, the identifier becomes
#' increasingly sensitive at higher count levels and flags Poisson-distributed data as overdispersed far
#' too often (Daniels et al., 2026). It is therefore switched off by default and should only be enabled
#' for backward compatibility or at low count levels. See \code{\link{Overdispersion.test}} for a
#' scale-invariant alternative, and note that no dispersion diagnostic should be relied on alone: the
#' data should always be inspected graphically as well, for example by plotting the group variances
#' against the group means together with the identity line.
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param get.contrasts.and.p.values Get each row of the contrast matrices evaluated
#' @param show.output Show/hide output
#' @param overdispersion.test Add the result of the overdispersion pre-test to the information object.
#' Use NULL or FALSE (default) to switch the pre-test off, TRUE to run it with the default settings of
#' \code{\link{Overdispersion.test}}, or a named list with any of the arguments alpha, bootstrap.runs,
#' pretest.mode, tolerance.factor, use.fixed.random.seed, include.zero.groups and zero.group.mean to
#' control it. Enabling the pre-test does not change the CPCAT results, since the state of the random
#' number generator is restored afterwards.
#' @return R object with results and information from CPCAT calculations
#' @references
#' Daniels, B.; Kotschik, P.; Walter-Rohde, S.; Graeff, T. (2026): When Counts Deviate: Robust Overdispersion Diagnostics for Ecotoxicity Data. Poster presented at the 36th SETAC Europe Annual Meeting, 17-21 May 2026, Maastricht, The Netherlands.
#'
#' Bretz, F.; Hothorn, T.; Westfall, P. (2010): Multiple comparisons using R. 1st Edition, Chapman and Hall/CRC, New York
#'
#' Chang, C.-H.; Pal, N.; Lin, J.-J. (2010): A Note on Comparing Several Poisson Means. In: Commun. Stat. Simul. Comput., 2010, 39(8), p. 1605-1627, https://doi.org/10.1080/03610918.2010.508860
#'
#' Lehmann, R.; Bachmann, J.; Maletzki, D.; Polleichtner, C.; Ratte, H.; Ratte, M. (2016): A new approach to overcome shortcomings with multiple testing of reproduction data in ecotoxicology. In: Stochastic Environmental Research and Risk Assessment, 2016, 30(3), p. 871-882, https://doi.org/10.1007/s00477-015-1079-4
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' # Test CPCAT
#' CPCAT(groups = Daphnia.counts$Concentration,
#'		 counts = Daphnia.counts$Number_Young,
#'		 control.name = NULL,
#'		 bootstrap.runs = 10000,
#'		 use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		 get.contrasts.and.p.values = FALSE,
#'		 show.output = TRUE)
#'
#' # Test CPCAT and add the overdispersion pre-test to the information object
#' CPCAT(groups = Daphnia.counts$Concentration,
#'		 counts = Daphnia.counts$Number_Young,
#'		 control.name = NULL,
#'		 bootstrap.runs = 1000,	# Caution: low number of bootstrap runs for testing
#'		 use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		 get.contrasts.and.p.values = FALSE,
#'		 show.output = TRUE,
#'		 overdispersion.test = list(bootstrap.runs = 500, use.fixed.random.seed = 123))
#' @export
CPCAT = function(groups,									# group vector
				 counts,									# vector with count data
				 control.name = NULL,						# character string with control group name
				 bootstrap.runs = 10000,					# number of bootstrap runs
				 hampel.threshold = NULL,					# threshold for Hampel identifier, NULL switches the identifier off
				 use.fixed.random.seed = NULL,              # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
				 get.contrasts.and.p.values = FALSE,		# get each row of the contrast matrices evaluated
				 show.output = TRUE,						# show/hide output
				 overdispersion.test = NULL) {				# add the overdispersion pre-test (NULL/FALSE = off, TRUE = defaults, list = settings)

	# check if there is count data for each replicate (length of count and group vectors) - groups[i] is one replicate
	if (length(groups) != length(counts)) {
		stop("Lengths of groups and counts don't match!")
	}
	# check format of input data
	if (!is.numeric(counts) | any(counts < 0, na.rm = TRUE)) {		#  | !all(counts == floor(counts))
		stop("Counts must be non-negative numeric values!")
	}

	# setup information to be stored
	info = data.frame(matrix(nrow = 0, ncol = 1))
	info = rbind(info, "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")

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
		info = rbind(info, paste0(dim(dat)[1] - dim(stats::na.omit(dat))[1], " rows with NA values were excluded!"))
	}
	dat = stats::na.omit(dat)

	# Check for over- and under-dispersion using the Hampel identifier, if a threshold was provided.
	# The identifier is switched off by default, see the description of the hampel.threshold parameter.
	if (!is.null(hampel.threshold)) {
		if (!is.numeric(hampel.threshold) | length(hampel.threshold) != 1) {
			stop("hampel.threshold must be a single numeric value or NULL!")
		}
		mean.dat = stats::aggregate(dat$Counts, by=list(dat$Groups), mean)$x
		var.dat = stats::aggregate(dat$Counts, by=list(dat$Groups), stats::var)$x
		hampel.value = var.dat - mean.dat
		if (min(hampel.value) < -hampel.threshold) {
			info = rbind(info, paste0("There was under-dispersed data identified in treatment(s) ",
									  paste0(paste0(treatments, " (HI: ", round(hampel.value, digits=1), ")")[which(hampel.value < -hampel.threshold)], collapse = ", "),
									  ". HI = Hampel Identifier."))
		}
		if (max(hampel.value) > hampel.threshold) {
			info = rbind(info, paste0("There was over-dispersed data identified in treatment(s) ",
									  paste0(paste0(treatments, " (HI: ", round(hampel.value, digits=1), ")")[which(hampel.value > hampel.threshold)], collapse = ", "),
									  ". HI = Hampel Identifier."))
		}
	}

	# Run the overdispersion pre-test if it was requested
	overdispersion.info = run.overdispersion.diagnostic(dat, overdispersion.test)
	if (!is.null(overdispersion.info)) {
		info = rbind(info, overdispersion.info)
	}

	# All hypotheses to be tested
	n = length(levels(dat$Groups))
	allhypotheses = CP.hypotheses(n = n - 1, treatment.names = treatments)

	# Transform list to table data.frame
	allhypothesescompact = numeric()
	for (l in 1:length(allhypotheses)) {
		allhypothesescompact = rbind(allhypothesescompact, allhypotheses[[l]])
	}

	# Only unique rows are selected
	allhypothesescompact = unique(allhypothesescompact)

	results = list()
	# Flag all hypotheses which have already been tested by assigning a p-value,
	# else p-value = -9999
	flagpvalues = matrix(-9999, nrow = nrow(allhypothesescompact), ncol = ncol(allhypothesescompact))
	pvalsCPCAT = rep(1, n - 1)

	# Fix seed for random numbers if desired (enables to reproduce results)
	if (!is.null(use.fixed.random.seed)) {
		if (!is.numeric(use.fixed.random.seed)) {
			stop("use.fixed.random.seed must be a numeric value or NULL.")
		}
		set.seed(use.fixed.random.seed)
	}

	for (j in 1:(n - 1)) {
		# Identify contrasts a p-value != -9999 has been assigned to
		# These contrasts must not be tested again
		contrasts = CP.hypotheses(n = n - 1, treatment.names = treatments)[[j]]
		matchingrows = numeric()
		for (i in 1:nrow(contrasts)) {
			matchingrows = c(matchingrows, which(apply(allhypothesescompact, 1, identical, contrasts[i, ])))
		}
		alreadyflaggedindex = which(flagpvalues[matchingrows, j] != -9999)

		# Shorten contrasts to be tested by elimination of already tested contrasts
		if (length(alreadyflaggedindex) > 0) {
			contrasts = contrasts[-alreadyflaggedindex, ]
		}

		# In the last step the contrast matrix reduces to a vector
		# Make it a matrix consisting of nrow = 1
		if (is.matrix(contrasts) == FALSE) {
			contrasts = matrix(contrasts, nrow = 1)
		}
		notflaggedindex = which(flagpvalues[matchingrows, j] == -9999)
		# Flag p-values which are still -9999
		# After CPCAT corresponding p-values will be != -9999
		tobeflagged = matchingrows[notflaggedindex]

		results[[j]] = CPCAT.Poisson.test(dat = dat,
										  contrastmatrix = contrasts,
										  bootstrap.runs = bootstrap.runs)[[1]]

		if (j == 1) {
			contrasts.and.p.values = results[[j]]
		} else {
			contrasts.and.p.values = rbind(contrasts.and.p.values, results[[j]])
		}

		# Write obtained p-values into column j of flagpvalues and find max p-value
		pvalshelp = results[[j]][, ncol(results[[j]])]
		flagpvalues[tobeflagged, j] = pvalshelp

		# Put together new p-values of reduced contrast matrix and relevant p-values in flagpvalues[, j]
		if (j > 1) {  # In step j = 1 all flagpvalues equal -9999
			pvalshelp2 = c(pvalshelp, flagpvalues[matchingrows[-notflaggedindex], j])
		} else {
			pvalshelp2 = pvalshelp
		}
		pvalsCPCAT[j] = max(pvalshelp2)

		# Copy p-values obtained so far to the next column of flagpvalues
		if (j < (n - 1)) {
			flagpvalues[, j + 1] = flagpvalues[, j]
		}
	}

	# Assign significance levels based on p-values
	significances = rep(NA, n - 1)
	for (j in 1:(n - 1)) {
		if (pvalsCPCAT[j] < 0.05) {
			if (pvalsCPCAT[j] < 0.01) {
				if (pvalsCPCAT[j] < 0.001) {
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
	NOEC = treatments[1]
	LOEC = treatments[2]
	for (j in 1:(n - 1)) {
		if (pvalsCPCAT[j] < 0.05) {
			break
		}
		NOEC = treatments[j]
		if (j == (n - 1)) {
			LOEC = NA
			break
		}
		LOEC = treatments[j+1]
	}
	info = rbind(info, paste0("NOEC: ", NOEC, ", LOEC: ", ifelse(is.na(LOEC), "outside tested dose/concentration", LOEC),
							  ". Assuming that any effects are adverse. Otherwise, NOEC and LOEC should be reconsidered."))

	# Compile results into a data.frame
	results = data.frame(Hypothesis = paste0("H0: ", treatments[1], " <-> ", treatments[2:n]), p.values = pvalsCPCAT, Signif. = significances)

	# Set header for information object
	colnames(info) = "Information and warnings:"

	# Show output if desired
	if (show.output) {
		if (get.contrasts.and.p.values) {
			print(structure(list(Contrasts=data.frame(contrasts.and.p.values), Results=results, Info=info)), row.names = F, quote = F, right = F)
		} else {
			print(structure(list(Results=results, Info=info)), row.names = F, quote = F, right = F)
		}
	}

	# Provide output as object even if not shown
	if (get.contrasts.and.p.values) {
		invisible(structure(list(Contrasts=data.frame(contrasts.and.p.values), Results=results, Info=info)))
	} else {
		invisible(structure(list(Results=results, Info=info)))
	}
}
