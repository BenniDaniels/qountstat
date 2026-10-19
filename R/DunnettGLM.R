#' @title Dunnett.GLM
#' @description When conducting statistical tests with multiple treatments, such as a control group and
#' increasing concentrations of a test substance, ANOVA and parametric post-hoc tests (e.g. Dunnett's test)
#' are commonly used. However, these tests require the assumptions of homogeneous variances and normally
#' distributed data. For count data (e.g. counts of animals), these assumptions are typically violated,
#' as the data are usually Poisson-distributed. The Dunnett.GLM function is based on a GLM followed by a
#' Dunnett test to the model estimates. It was implemented to serve as an alternative approach to CPCAT
#' while using a Quasi-Poisson regression (default), a negative binomial regression or a classical
#' Poisson regression.
#' The basic approach from Hothorn and Kluxen (2020) was adjusted to overcome methodological issues (see
#' description of 'zero.treatment.action parameter').
#' For details on the structure of the input data, please refer to the dataset 'Daphnia.counts' provided
#' alongside this package.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param zero.treatment.action Method for dealing with treatments only containing zeros (use either
#' "identity.link" or "log(x+1)"). The method is only used if the data set contains dose/concentration
#' groups that exclusively contain zero values (since the basic method provides for a logarithmic
#' transformation of the data averages, it would lead to incorrect results). To deal with this
#' methodological shortcoming, two options were implemented. The 'identity.link' option: the 'identity'
#' link is used in the GLM instead of the 'log' link, i.e. the data are no longer transformed. The 'log(x+1)'
#' option: The 'log' link is retained and 1 is added to each count value at the start of the procedure so
#' that the subsequent log-transformation can be carried out without any problems. Note that both options
#' may slightly distort the results.
#' @param show.output Show/hide output
#' @param family Distribution family used for the GLM (use either "quasipoisson", "nb" or "poisson").
#' The default "quasipoisson" fits a Quasi-Poisson regression, i.e. the variance is modelled as a
#' constant multiple (the dispersion parameter) of the mean. The option "nb" fits a negative binomial
#' regression via \code{\link[MASS]{glm.nb}}, i.e. the variance is modelled as mu + mu^2/theta, with the
#' dispersion parameter theta being estimated from the data by maximum likelihood. Both of these options
#' are suited for overdispersed count data. The negative binomial model is a fully parametric likelihood
#' model and is often preferable if the overdispersion increases more than proportionally with the mean,
#' whereas the Quasi-Poisson model makes weaker distributional assumptions. The option "poisson" fits a
#' classical Poisson regression, i.e. the variance is assumed to equal the mean and the dispersion
#' parameter is fixed at 1 rather than estimated. This option is the regression counterpart of CPCAT and
#' is intended for data that are compatible with a Poisson distribution; it should only be used if that
#' assumption has been checked, for example with \code{\link{Overdispersion.test}} and by inspecting the
#' data graphically (see the note below), because standard errors are underestimated and the inference
#' becomes anticonservative if the data are overdispersed. Note that "poisson" and "nb" require
#' (approximately) integer counts, while "quasipoisson" also accepts non-integer values.
#' @param overdispersion.test Add the result of the overdispersion pre-test to the information object.
#' Use NULL or FALSE (default) to switch the pre-test off, TRUE to run it with the default settings of
#' \code{\link{Overdispersion.test}}, or a named list with any of the arguments alpha, bootstrap.runs,
#' pretest.mode, tolerance.factor, use.fixed.random.seed, include.zero.groups and zero.group.mean to
#' control it. Enabling the pre-test does not change the Dunnett.GLM results, since the state of the
#' random number generator is restored afterwards. The pre-test is primarily a decision aid for choosing
#' between a Poisson-based method such as CPCAT and the regression approach implemented here; within
#' Dunnett.GLM it serves as documentation of the dispersion actually present in the data.
#' @return R object with results and information from Dunnett.GLM calculations
#' @details
#' Note that no dispersion diagnostic should be the sole basis for this decision. Any such test reduces
#' the mean-variance relationship of a data set to a single number and, at the replication typical for
#' ecotoxicological studies, has limited power, so that a non-significant result is weak evidence for the
#' Poisson assumption rather than proof of it. The data should therefore always be inspected graphically
#' as well, for example by plotting the group variances against the group means together with the
#' identity line, or the residuals against the fitted values, and by looking at the per-group
#' variance-to-mean ratios that Overdispersion.test reports.
#'
#' The information object that is returned alongside the test results reports the distribution and the
#' link function that were used, together with the dispersion parameter. For family = "quasipoisson"
#' this is the estimated dispersion of the Quasi-Poisson model, for family = "nb" it is the estimated
#' theta of the negative binomial model including its standard error, and for family = "poisson" it is
#' noted that the dispersion is fixed at 1 by the model assumption.
#'
#' Note on zero treatments and the negative binomial model: zero.treatment.action only takes effect if
#' a treatment consists exclusively of zeros. In that situation family = "nb" cannot be combined with
#' zero.treatment.action = "identity.link", and Dunnett.GLM stops with a message that says so. The
#' reason is that an identity link forces the fitted mean of such a treatment to zero. A treatment with
#' a mean of zero also has a variance of zero and therefore carries no information on the dispersion,
#' so the negative binomial dispersion parameter is not identifiable at that point, and its maximum
#' likelihood estimation divides by the fitted mean. Whether that mean is represented as a small
#' positive number, as exactly zero or as a small negative number depends on floating point rounding
#' and differs between platforms, which would make the result machine-dependent.
#'
#' Use zero.treatment.action = "log(x+1)" together with family = "nb" instead. Adding one to every
#' count removes the degenerate mean, the log link is retained and the dispersion is estimated without
#' difficulty. If the identity link is required, use family = "quasipoisson" or family = "poisson",
#' neither of which estimates a dispersion parameter that depends on the fitted mean in this way. Note
#' that the two options also differ in interpretation: with a log link the Dunnett contrasts compare
#' the treatments with the control on a multiplicative scale, with an identity link on an additive one.
#' @references
#' Hothorn, L.; Kluxen, F. (2020): Statistical analysis of no observed effect concentrations or levels in eco-toxicological assays with overdispersed count endpoints. In: bioRxiv, 2020, https://doi.org/10.1101/2020.01.15.907881
#'
#' Venables, W. N.; Ripley, B. D. (2002): Modern Applied Statistics with S. Fourth Edition. Springer, New York.
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' # Test Dunnett.GLM with 'identity.link' option
#' Dunnett.GLM(groups = Daphnia.counts$Concentration,
#' 			   counts = Daphnia.counts$Number_Young,
#' 			   control.name = NULL,
#' 			   zero.treatment.action = "identity.link",
#' 			   show.output = TRUE)
#'
#' # Test Dunnett.GLM with 'log(x+1)' option
#' Dunnett.GLM(groups = Daphnia.counts$Concentration,
#'			   counts = Daphnia.counts$Number_Young,
#'			   control.name = NULL,
#'			   zero.treatment.action = "log(x+1)",
#'			   show.output = TRUE)
#'
#' # Test Dunnett.GLM using a negative binomial distribution instead of a Quasi-Poisson distribution.
#' # Note that the 'log(x+1)' option is used here: the example data contain a treatment consisting
#' # only of zeros, and the negative binomial model cannot be combined with the identity link in that
#' # situation (see the note on zero treatments below).
#' Dunnett.GLM(groups = Daphnia.counts$Concentration,
#'			   counts = Daphnia.counts$Number_Young,
#'			   control.name = NULL,
#'			   zero.treatment.action = "log(x+1)",
#'			   show.output = TRUE,
#'			   family = "nb")
#'
#' # Test Dunnett.GLM using a classical Poisson distribution, the regression counterpart of CPCAT
#' Dunnett.GLM(groups = Daphnia.counts$Concentration,
#'			   counts = Daphnia.counts$Number_Young,
#'			   control.name = NULL,
#'			   zero.treatment.action = "log(x+1)",
#'			   show.output = TRUE,
#'			   family = "poisson")
#'
#' # Test Dunnett.GLM and add the overdispersion pre-test to the information object
#' Dunnett.GLM(groups = Daphnia.counts$Concentration,
#'			   counts = Daphnia.counts$Number_Young,
#'			   control.name = NULL,
#'			   zero.treatment.action = "log(x+1)",
#'			   show.output = TRUE,
#'			   overdispersion.test = list(bootstrap.runs = 500, use.fixed.random.seed = 123))
#' @export
Dunnett.GLM = function(	groups,										# group vector
						counts,										# vector with count data
						control.name = NULL,						# character string with control group name
						zero.treatment.action = "identity.link",	# method for dealing with treatments only containing zeros (alternative: "log(x+1)")
						show.output = TRUE,							# show/hide output
						family = "quasipoisson",					# distribution family used for the GLM (alternatives: "nb", "poisson")
						overdispersion.test = NULL){				# add the overdispersion pre-test (NULL/FALSE = off, TRUE = defaults, list = settings)

	#library(multcomp)

	# do some error handling
	if (length(groups) != length(counts)) {
		stop("Lengths of groups and counts don't match!")
	}

	if (!is.numeric(counts) | any(counts < 0, na.rm = TRUE)) {	 # | !all(counts == floor(counts))
		stop("Counts must be non-negative numeric values!")
	}

	if (zero.treatment.action != "identity.link" & zero.treatment.action != "log(x+1)") {
		stop("Parameter zero.treatment.action must be either \"identity.link\" or \"log(x+1)\"!")
	}

	if (length(family) != 1 | !is.character(family)) {
		stop("Parameter family must be either \"quasipoisson\", \"nb\" or \"poisson\"!")
	}
	if (!is.element(family, c("quasipoisson", "nb", "poisson"))) {
		stop("Parameter family must be either \"quasipoisson\", \"nb\" or \"poisson\"!")
	}

	# The negative binomial and the Poisson model are discrete probability models and require
	# integer counts, whereas the Quasi-Poisson model only specifies the mean-variance relation
	if (is.element(family, c("nb", "poisson")) & !all(stats::na.omit(counts) == floor(stats::na.omit(counts)))) {
		stop(paste0("Counts must be integer values if family = \"", family, "\" is used!"))
	}

	# setup information to be stored
	info = data.frame(matrix(nrow = 0, ncol = 1))

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

	# Run the overdispersion pre-test if it was requested
	overdispersion.info = run.overdispersion.diagnostic(dat, overdispersion.test)
	if (!is.null(overdispersion.info)) {
		info = rbind(info, overdispersion.info)
	}

	agg = stats::aggregate(dat$Counts, by=list(dat$Groups), mean)

	# Fit the GLM using the selected distribution family. The negative binomial model ("nb") is fitted
	# via MASS::glm.nb, which estimates the dispersion parameter theta by maximum likelihood; the
	# Quasi-Poisson model ("quasipoisson") and the Poisson model ("poisson") are fitted via stats::glm.
	# The handling of treatments that exclusively contain zeros (see zero.treatment.action) is identical
	# for all families, since it is a property of the link function and not of the distribution.
	if (min(agg$x) > 0) {
		# Dunnett GLM with log link-function
		if (family == "quasipoisson") {
			mod = stats::glm(Counts~Groups, data=dat, family=stats::quasipoisson(link = "log"))
		} else if (family == "poisson") {
			mod = stats::glm(Counts~Groups, data=dat, family=stats::poisson(link = "log"))
		} else {
			mod = try(MASS::glm.nb(Counts~Groups, data=dat, link = "log"), silent = T)
			if ("try-error" %in% class(mod)){
				stop("Error in GLM calculation. Consider to change argument family to \"quasipoisson\".")
			}
		}
	} else {
		info = rbind(info, paste0("A treatment contained only zeros, hence, the zero.treatment.action \"", zero.treatment.action, "\" was applied."))

		# Use either the identity link or a data transformation if one treatment mean is zero (standard log link would fail)
		if (zero.treatment.action == "identity.link") {
			# The negative binomial model is not available here. With an identity link the fitted mean
			# of an all-zero treatment is zero, and the negative binomial dispersion parameter is not
			# identifiable at that point: such a treatment has zero variance as well and therefore
			# carries no information on the dispersion, while the maximum likelihood estimation of
			# theta divides by the fitted mean. Whether that mean comes out as a small positive number,
			# as exactly zero or as a small negative number is decided by floating point rounding and
			# differs between platforms, so the fit would succeed or fail depending on the machine.
			# The combination is therefore rejected with an explicit message instead.
			if (family == "nb") {
				stop(paste0("A treatment contains only zeros. The negative binomial model cannot be ",
							"combined with zero.treatment.action = \"identity.link\" in this case, ",
							"because its dispersion parameter is not identifiable when a fitted mean ",
							"is zero. Use zero.treatment.action = \"log(x+1)\" together with ",
							"family = \"nb\", or keep the identity link and use family = ",
							"\"quasipoisson\" or family = \"poisson\"."))
			}

			# use another link function (link="identity", i.e. E(Y) = b0 + b1 * X)
			if (family == "quasipoisson") {
				mod = try(stats::glm(Counts~Groups, data=dat, family=stats::quasipoisson(link = "identity")), silent = T)
			} else {
				mod = try(stats::glm(Counts~Groups, data=dat, family=stats::poisson(link = "identity")), silent = T)
			}
			if ("try-error" %in% class(mod)){
				stop("Error in GLM calculation. Consider to change argument zero.treatment.action to \"log(x+1)\".")
			}
		} else {
			# modify data to always get positive estimates using log link (link="log", i.e. log(E(Y)) = b0 + b1 * X)
			dat$Counts = dat$Counts + 1
			if (family == "quasipoisson") {
				mod = try(stats::glm(dat$Counts~Groups, data=dat, family=stats::quasipoisson(link = "log")), silent = T)
			} else if (family == "poisson") {
				mod = try(stats::glm(Counts~Groups, data=dat, family=stats::poisson(link = "log")), silent = T)
			} else {
				mod = try(MASS::glm.nb(Counts~Groups, data=dat, link = "log"), silent = T)
			}
			if ("try-error" %in% class(mod)){
				stop("Error in GLM calculation. Consider to change argument zero.treatment.action to \"identity.link\".")
			}
		}
	}

	results = summary(multcomp::glht(mod, linfct = multcomp::mcp(Groups = "Dunnett"), alternative="two.sided"))

	# Report the distribution family that was used together with its dispersion parameter
	if (family == "quasipoisson") {
		info = rbind(info, paste0("Distribution: Quasi-Poisson (family = \"quasipoisson\"), link = \"",
								  mod$family$link, "\", estimated dispersion parameter = ",
								  sprintf("%.4f", summary(mod)$dispersion), "."))
	} else if (family == "poisson") {
		info = rbind(info, paste0("Distribution: Poisson (family = \"poisson\"), link = \"",
								  mod$family$link, "\", dispersion parameter fixed at 1 by the model ",
								  "assumption (variance = mean)."))
	} else {
		info = rbind(info, paste0("Distribution: negative binomial (family = \"nb\"), link = \"",
								  mod$family$link, "\", estimated dispersion parameter theta = ",
								  sprintf("%.4f", mod$theta), " (SE = ", sprintf("%.4f", mod$SE.theta), ")."))
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
