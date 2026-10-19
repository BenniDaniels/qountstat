#' @title Overdispersion.test
#' @description Poisson-based methods for the determination of no/lowest observed (adverse) effect
#' concentrations from count data, such as CPCAT, rely on the assumption that the variance of the counts
#' equals their mean. Ecotoxicological count data frequently violate this assumption (overdispersion),
#' which invalidates the subsequent inference. Overdispersion.test provides a model-agnostic pre-test
#' that checks this assumption before a method is chosen. Model-agnostic means that no generalised linear
#' model has to be specified and fitted beforehand; the Poisson distribution is only used as a reference
#' framework. The test statistic is the mean of the squared within-group Poisson-standardised residuals,
#' T = mean(z^2) with z_ij = (y_ij - mu_i) / sqrt(mu_i), which is scale-invariant and therefore comparable
#' across designs with different mean levels. The decision threshold is not taken from an asymptotic
#' distribution but is calibrated by Monte Carlo simulation for the design at hand, which makes the test
#' applicable to the small sample sizes that are typical for ecotoxicological studies. The test is
#' one-sided: only unusually large values of T lead to rejection, i.e. underdispersion is not detected.
#' For details on the structure of the input data, please refer to the dataset 'Daphnia.counts' provided
#' alongside this package.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param alpha Significance level of the one-sided test
#' @param bootstrap.runs Number of Monte Carlo simulations used to calibrate the decision threshold
#' @param pretest.mode Reference distribution used for the calibration (use either "tolerant" or
#' "strict"). The "strict" option calibrates against a pure Poisson distribution, i.e. any deviation from
#' variance = mean is treated as overdispersion. The "tolerant" option calibrates against a negative
#' binomial distribution with a tolerated variance-to-mean ratio (see tolerance.factor), i.e. mild
#' deviations from the Poisson assumption are accepted. The tolerant option reflects that Poisson-based
#' methods are robust against small deviations and that a threshold is needed above which they are not.
#' @param tolerance.factor Tolerated variance-to-mean ratio used if pretest.mode = "tolerant". The
#' negative binomial size parameter is set to k = mu / (tolerance.factor - 1), which yields a variance of
#' mu * tolerance.factor, i.e. the tolerance factor is exactly the tolerated variance-to-mean ratio. Must
#' be greater than 1. Ignored if pretest.mode = "strict".
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param include.zero.groups Include treatment groups whose counts are all zero. Such groups carry no
#' information on dispersion and are excluded by default.
#' @param zero.group.mean Value that replaces a group mean of zero in the denominator of the standardised
#' residuals. Only used if include.zero.groups = TRUE.
#' @param show.output Show/hide output
#' @return R object with results, group statistics and information from the overdispersion test
#' @details
#' The observed statistic T is compared with a critical value c_alpha that is obtained as the
#' (1 - alpha) quantile of the Monte Carlo null distribution of T. Overdispersion is reported if
#' T > c_alpha. In addition, a one-sided Monte Carlo p-value is calculated as
#' (1 + number of simulated values >= T) / (1 + bootstrap.runs).
#'
#' Note on the interpretation of T: under a Poisson distribution with KNOWN group means the expected
#' value of z^2 is 1. Since the group means are estimated from the same data, the expected value of T is
#' (N - G) / N, with N being the number of observations and G the number of groups, i.e. it is smaller
#' than 1 and depends on the design (for example 0.833 for six replicates per group). This does not
#' affect the validity of the test, because the Monte Carlo calibration uses exactly the same estimator
#' and the bias therefore cancels out. To provide a directly interpretable quantity, the results
#' additionally report the bias-corrected dispersion ratio phi = T * N / (N - G), which has an expected
#' value of 1 under a Poisson distribution and can be read as a pooled variance-to-mean ratio. Since phi
#' is a strictly monotone transformation of T, the test decision is identical in both scales.
#'
#' The decision of this test should not be used as the only criterion for the choice of a method. The
#' statistic reduces the mean-variance relationship of a whole data set to a single number, and at the
#' replication typical for ecotoxicological studies the test has limited power, so that a
#' non-significant result is weak evidence for the Poisson assumption rather than proof of it. The group
#' statistics returned alongside the test report the mean, the variance and the variance-to-mean ratio of
#' every treatment, and these should be inspected as well, together with a graphical check such as a plot
#' of the group variances against the group means with the identity line added, or of the residuals
#' against the fitted values. A single treatment with a markedly deviating variance-to-mean ratio is
#' visible there but may well leave the pooled statistic non-significant.
#'
#' Groups containing fewer than two observations are excluded, since a within-group dispersion cannot be
#' evaluated for them.
#' @references
#' Daniels, B.; Kotschik, P.; Walter-Rohde, S.; Graeff, T. (2026): When Counts Deviate: Robust Overdispersion Diagnostics for Ecotoxicity Data. Poster presented at the 36th SETAC Europe Annual Meeting, 17-21 May 2026, Maastricht, The Netherlands.
#'
#' Cameron, A. C.; Trivedi, P. K. (1990): Regression-based tests for overdispersion in the Poisson model. In: Journal of Econometrics, 1990, 46(3), p. 347-364, https://doi.org/10.1016/0304-4076(90)90014-K
#'
#' Lehmann, R.; Bachmann, J.; Maletzki, D.; Polleichtner, C.; Ratte, H.; Ratte, M. (2016): A new approach to overcome shortcomings with multiple testing of reproduction data in ecotoxicology. In: Stochastic Environmental Research and Risk Assessment, 2016, 30(3), p. 871-882, https://doi.org/10.1007/s00477-015-1079-4
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' # Test for overdispersion, tolerating a variance-to-mean ratio of up to 2
#' Overdispersion.test(groups = Daphnia.counts$Concentration,
#'		 counts = Daphnia.counts$Number_Young,
#'		 alpha = 0.05,
#'		 bootstrap.runs = 500,	# Caution: low number of simulations for testing
#'		 pretest.mode = "tolerant",
#'		 tolerance.factor = 2,
#'		 use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		 show.output = TRUE)
#'
#' # Test against a pure Poisson distribution, i.e. without tolerating any deviation
#' Overdispersion.test(groups = Daphnia.counts$Concentration,
#'		 counts = Daphnia.counts$Number_Young,
#'		 alpha = 0.05,
#'		 bootstrap.runs = 500,	# Caution: low number of simulations for testing
#'		 pretest.mode = "strict",
#'		 use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		 show.output = TRUE)
#' @export
Overdispersion.test = function(groups,								# group vector
							   counts,								# vector with count data
							   alpha = 0.05,						# significance level of the one-sided test
							   bootstrap.runs = 5000,				# number of Monte Carlo simulations for the calibration
							   pretest.mode = "tolerant",			# reference distribution ("tolerant" or "strict")
							   tolerance.factor = 2,				# tolerated variance-to-mean ratio (only if "tolerant")
							   use.fixed.random.seed = NULL,		# fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
							   include.zero.groups = FALSE,			# include treatment groups containing only zeros
							   zero.group.mean = 0.5,				# replacement for a group mean of zero
							   show.output = TRUE) {				# show/hide output

	# do some error handling
	if (length(groups) != length(counts)) {
		stop("Lengths of groups and counts don't match!")
	}
	if (!is.numeric(counts) | any(counts < 0, na.rm = TRUE)) {
		stop("Counts must be non-negative numeric values!")
	}
	if (!is.numeric(alpha) | length(alpha) != 1) {
		stop("alpha must be a single numeric value between 0 and 1!")
	}
	if (alpha <= 0 | alpha >= 1) {
		stop("alpha must be a single numeric value between 0 and 1!")
	}
	if (!is.numeric(bootstrap.runs) | length(bootstrap.runs) != 1) {
		stop("bootstrap.runs must be a single numeric value of at least 100!")
	}
	if (bootstrap.runs < 100) {
		stop("bootstrap.runs must be a single numeric value of at least 100!")
	}
	if (length(pretest.mode) != 1 | !is.character(pretest.mode)) {
		stop("Parameter pretest.mode must be either \"tolerant\" or \"strict\"!")
	}
	if (pretest.mode != "tolerant" & pretest.mode != "strict") {
		stop("Parameter pretest.mode must be either \"tolerant\" or \"strict\"!")
	}
	if (pretest.mode == "tolerant") {
		if (!is.numeric(tolerance.factor) | length(tolerance.factor) != 1) {
			stop("tolerance.factor must be a single numeric value greater than 1!")
		}
		if (tolerance.factor <= 1) {
			stop("tolerance.factor must be a single numeric value greater than 1!")
		}
	}
	if (!is.numeric(zero.group.mean) | length(zero.group.mean) != 1) {
		stop("zero.group.mean must be a single positive numeric value!")
	}
	if (zero.group.mean <= 0) {
		stop("zero.group.mean must be a single positive numeric value!")
	}

	# setup information to be stored
	info = data.frame(matrix(nrow = 0, ncol = 1))

	# Re-structure the input to a data frame
	dat = data.frame(Counts = counts, Groups = groups)

	# Convert groups column to a factor, specifying the desired order of levels
	dat$Groups = factor(dat$Groups, levels = unique(dat$Groups))

	# Exit if not enough data left
	if (dim(stats::na.omit(dat))[1] < 2) {
		stop("Too few valid data!")
	}
	if (dim(dat)[1] != dim(stats::na.omit(dat))[1]) {
		info = rbind(info, paste0(dim(dat)[1] - dim(stats::na.omit(dat))[1], " rows with NA values were excluded!"))
	}
	dat = stats::na.omit(dat)
	dat$Groups = droplevels(dat$Groups)

	# Split the counts by treatment group
	group.counts = split(dat$Counts, dat$Groups)

	# Groups with less than two observations carry no information on the within-group dispersion
	too.small = names(group.counts)[sapply(group.counts, length) < 2]
	if (length(too.small) > 0) {
		info = rbind(info, paste0("Treatment(s) ", paste0(too.small, collapse = ", "),
								  " contained less than two observations and were excluded."))
		group.counts = group.counts[sapply(group.counts, length) >= 2]
	}

	# Groups containing only zeros carry no information on the dispersion either
	if (!include.zero.groups) {
		zero.groups = names(group.counts)[sapply(group.counts, mean) <= 0]
		if (length(zero.groups) > 0) {
			info = rbind(info, paste0("Treatment(s) ", paste0(zero.groups, collapse = ", "),
									  " contained only zeros and were excluded (see include.zero.groups)."))
			group.counts = group.counts[sapply(group.counts, mean) > 0]
		}
	}

	if (length(group.counts) < 1) {
		stop("No treatment groups left for the evaluation of overdispersion!")
	}

	# Group means used as the Poisson expectation, group sizes and observed statistic
	mu.vector = sapply(group.counts, mean)
	mu.vector[mu.vector <= 0] = zero.group.mean
	n.vector = sapply(group.counts, length)
	T.observed = MSSR.statistic(group.counts, zero.group.mean)

	# Number of observations and groups define the bias correction of the statistic
	n.observations = sum(n.vector)
	n.groups = length(group.counts)
	correction = n.observations / (n.observations - n.groups)

	# Fix seed for random numbers if desired (enables to reproduce results)
	if (!is.null(use.fixed.random.seed)) {
		if (!is.numeric(use.fixed.random.seed)) {
			stop("use.fixed.random.seed must be a numeric value or NULL.")
		}
		set.seed(use.fixed.random.seed)
	}

	# Monte Carlo calibration of the decision threshold for the design at hand
	T.null = calibrate.c.alpha(mu.vector = mu.vector,
							   n.vector = n.vector,
							   bootstrap.runs = bootstrap.runs,
							   pretest.mode = pretest.mode,
							   tolerance.factor = tolerance.factor,
							   zero.group.mean = zero.group.mean)

	c.alpha = unname(stats::quantile(T.null, 1 - alpha, type = 8))
	p.value = (1 + sum(T.null >= T.observed)) / (1 + bootstrap.runs)
	overdispersion.detected = T.observed > c.alpha

	# Treatment-level statistics
	group.statistics = data.frame(Group = names(group.counts),
								  n = as.vector(n.vector),
								  Mean = round(as.vector(sapply(group.counts, mean)), digits = 4),
								  Variance = round(as.vector(sapply(group.counts, stats::var)), digits = 4),
								  Var.mean.ratio = round(as.vector(sapply(group.counts, stats::var) / mu.vector), digits = 4),
								  Mean.z2 = round(as.vector(sapply(group.counts, function(y) {
								  	mu = mean(y)
								  	if (mu <= 0) mu = zero.group.mean
								  	mean(((y - mean(y)) / sqrt(mu))^2)
								  })), digits = 4))
	rownames(group.statistics) = NULL

	# Compile results into a data.frame
	results = data.frame(Parameter = c("Groups", "Observations", "Observed T", "c_alpha",
									   "Dispersion ratio (phi)", "Critical dispersion ratio",
									   "Monte Carlo p-value", "Decision"),
						 Value = c(n.groups,
								   n.observations,
								   round(T.observed, digits = 4),
								   round(c.alpha, digits = 4),
								   round(T.observed * correction, digits = 4),
								   round(c.alpha * correction, digits = 4),
								   round(p.value, digits = 4),
								   ifelse(overdispersion.detected, "Overdispersion detected", "No overdispersion detected")))

	# Describe the reference the data were compared against
	if (pretest.mode == "tolerant") {
		info = rbind(info, paste0("One-sided test at alpha = ", alpha, " against a negative binomial reference with a tolerated ",
								  "variance-to-mean ratio of ", tolerance.factor, " (pretest.mode = \"tolerant\"), ",
								  bootstrap.runs, " Monte Carlo simulations."))
	} else {
		info = rbind(info, paste0("One-sided test at alpha = ", alpha, " against a pure Poisson reference ",
								  "(pretest.mode = \"strict\"), ", bootstrap.runs, " Monte Carlo simulations."))
	}
	info = rbind(info, paste0("Expected value of T under the reference distribution is (N-G)/N = ",
							  round(1 / correction, digits = 4), ", the dispersion ratio phi = T * N/(N-G) is ",
							  "scaled to an expected value of 1 under a Poisson distribution."))
	info = rbind(info, paste0("Underdispersion is not evaluated, the test is one-sided."))

	# Set header for information object
	colnames(info) = "Information and warnings:"

	# Provide output
	output = list(Results = results, Group.statistics = group.statistics, Info = info,
				  null.distribution = T.null)
	if (show.output) {
		print(structure(output[c("Results", "Group.statistics", "Info")]), row.names = F, quote = F, right = F)
	}
	invisible(structure(output))
}


#' @title MSSR statistic
#' @description Mean of the squared within-group Poisson-standardised residuals,
#' T = 1/N * sum_i sum_j ((y_ij - mu_i) / sqrt(mu_i))^2. Internal helper of Overdispersion.test.
#' @param group.counts List with one vector of counts per treatment group
#' @param zero.group.mean Value that replaces a group mean of zero in the denominator
#' @return Numeric value of the test statistic
MSSR.statistic = function(group.counts, zero.group.mean = 0.5) {
	z2 = unlist(lapply(group.counts, function(y) {
		mu = mean(y)
		# a group mean of zero would lead to a division by zero, hence a continuity value is used
		if (mu <= 0) mu = zero.group.mean
		# the residuals are taken around the observed group mean, the Poisson expectation is only
		# used to standardise their scale
		((y - mean(y)) / sqrt(mu))^2
	}), use.names = FALSE)
	mean(z2)
}


#' @title Monte Carlo calibration of the decision threshold
#' @description Simulates the null distribution of the MSSR statistic for a given design. Internal
#' helper of Overdispersion.test.
#' @param mu.vector Vector with the expected value per treatment group
#' @param n.vector Vector with the number of observations per treatment group
#' @param bootstrap.runs Number of Monte Carlo simulations
#' @param pretest.mode Reference distribution ("tolerant" or "strict")
#' @param tolerance.factor Tolerated variance-to-mean ratio (only used if pretest.mode = "tolerant")
#' @param zero.group.mean Value that replaces a group mean of zero in the denominator
#' @return Numeric vector with the simulated null distribution of the test statistic
calibrate.c.alpha = function(mu.vector, n.vector, bootstrap.runs, pretest.mode = "tolerant",
							 tolerance.factor = 2, zero.group.mean = 0.5) {

	n.groups = length(mu.vector)

	# With Var = mu + mu^2/k a target variance-to-mean ratio of phi is obtained for k = mu/(phi-1),
	# i.e. the tolerance factor is exactly the tolerated variance-to-mean ratio
	nb.size = if (pretest.mode == "strict") rep(Inf, n.groups) else mu.vector / (tolerance.factor - 1)

	T.null = numeric(bootstrap.runs)
	for (b in seq_len(bootstrap.runs)) {
		simulated = vector("list", n.groups)
		for (g in seq_len(n.groups)) {
			if (is.infinite(nb.size[g])) {
				simulated[[g]] = stats::rpois(n.vector[g], lambda = mu.vector[g])
			} else {
				simulated[[g]] = stats::rnbinom(n.vector[g], mu = mu.vector[g], size = nb.size[g])
			}
		}
		T.null[b] = MSSR.statistic(simulated, zero.group.mean)
	}
	T.null
}


#' @title Overdispersion diagnostic within a test function
#' @description Runs Overdispersion.test from within CPCAT or Dunnett.GLM and returns a single line of
#' text for their information object. The state of the random number generator is saved and restored, so
#' that enabling the diagnostic does not change the results of the calling function. Internal helper.
#' @param dat Data frame with the columns Counts and Groups
#' @param settings Either NULL/FALSE (diagnostic switched off), TRUE (default settings) or a named list
#' with arguments that are passed on to Overdispersion.test
#' @return Character string to be added to the information object, or NULL if switched off
run.overdispersion.diagnostic = function(dat, settings) {

	# switched off
	if (is.null(settings)) {
		return(NULL)
	}
	if (is.logical(settings)) {
		if (length(settings) != 1 | is.na(settings[1])) {
			stop("Parameter overdispersion.test must be NULL, TRUE, FALSE or a named list of settings!")
		}
		if (!settings) {
			return(NULL)
		}
		settings = list()
	}
	if (!is.list(settings)) {
		stop("Parameter overdispersion.test must be NULL, TRUE, FALSE or a named list of settings!")
	}

	# only arguments of Overdispersion.test may be passed on
	allowed = c("alpha", "bootstrap.runs", "pretest.mode", "tolerance.factor",
				"use.fixed.random.seed", "include.zero.groups", "zero.group.mean")
	if (length(settings) > 0) {
		if (is.null(names(settings)) | any(names(settings) == "")) {
			stop("All settings in overdispersion.test must be named!")
		}
		if (!all(names(settings) %in% allowed)) {
			stop(paste0("Unknown setting(s) in overdispersion.test: ",
						paste0(names(settings)[!names(settings) %in% allowed], collapse = ", "),
						". Allowed are: ", paste0(allowed, collapse = ", "), "."))
		}
	}

	# The diagnostic draws random numbers. The state of the random number generator is restored
	# afterwards so that the results of the calling function are not affected by switching it on.
	seed.exists = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	if (seed.exists) {
		seed.backup = get(".Random.seed", envir = globalenv(), inherits = FALSE)
	}
	on.exit({
		if (seed.exists) {
			assign(".Random.seed", seed.backup, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)

	settings$groups = dat$Groups
	settings$counts = dat$Counts
	settings$show.output = FALSE
	result = do.call(Overdispersion.test, settings)

	values = stats::setNames(result$Results$Value, result$Results$Parameter)
	paste0("Overdispersion pre-test (Daniels et al. 2026): ", values[["Decision"]],
		   " (T = ", values[["Observed T"]],
		   ", c_alpha = ", values[["c_alpha"]],
		   ", dispersion ratio phi = ", values[["Dispersion ratio (phi)"]],
		   ", Monte Carlo p = ", values[["Monte Carlo p-value"]], ").")
}
