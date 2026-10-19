#' @title Dunnett.GLM bootstrap MDD (bMDD)
#' @description The basic idea of the calculation of bootstrap MDD (bMDD) using the Dunnett.GLM approach
#' is to shift the lambda parameter of Poisson distribution until there is a certain proportion of
#' results significantly different from the control.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param shift.step Step of shift (negative as a reduction is assumed)
#' @param bootstrap.runs Number of bootstrap runs
#' @param power Proportion of bootstrap.runs that return significant differences
#' @param max.iterations Max. number of iterations to not get stuck in the while loop
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param Dunnett.GLM.zero.treatment.action Dunnett.GLM method to be used for treatments only containing zeros
#' @param show.progress Show progress for each shift of lambda
#' @param show.results Show results
#' @param Dunnett.GLM.family Distribution family to be used for the GLM within the Dunnett.GLM method
#' (use either "quasipoisson", "nb" or "poisson", see \code{\link{Dunnett.GLM}}). Note that this argument only
#' affects the model that is fitted to the resampled data; the bootstrap samples themselves are always
#' drawn from a Poisson distribution.
#' @return Data frame with results from bMDD analysis
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' \donttest{
#' # Test Dunnett.GLM bootstrap MDD
#' Dunnett.GLM.bMDD(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		shift.step = -1,		# Caution: big step size for testing
#'		bootstrap.runs = 5,	# Caution: low number of bootstrap runs for testing
#'		power = 0.8,
#'		max.iterations = 1000,
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE)
#'
#' # Test Dunnett.GLM bootstrap MDD using a negative binomial distribution
#' Dunnett.GLM.bMDD(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		shift.step = -1,		# Caution: big step size for testing
#'		bootstrap.runs = 5,	# Caution: low number of bootstrap runs for testing
#'		power = 0.8,
#'		max.iterations = 1000,
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE,
#'		Dunnett.GLM.family = "nb")
#'
#' # Test Dunnett.GLM bootstrap MDD using a classical Poisson distribution
#' Dunnett.GLM.bMDD(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		shift.step = -1,		# Caution: big step size for testing
#'		bootstrap.runs = 5,	# Caution: low number of bootstrap runs for testing
#'		power = 0.8,
#'		max.iterations = 1000,
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE,
#'		Dunnett.GLM.family = "poisson")
#' }
#' @export
Dunnett.GLM.bMDD = function(groups,							# group vector
							counts,							# vector with count data
							control.name = NULL,			# character string with control group name
							alpha = 0.05,					# significance level
							shift.step = -0.25,				# step of shift (negative as a reduction is assumed)
							bootstrap.runs = 200,			# number of bootstrap runs (draw Poisson data n times)
							power = 0.8,					# proportion of bootstrap.runs that return significant differences
							max.iterations = 1000,			# max number of iterations to not get stuck in the while loop
							use.fixed.random.seed = NULL,   # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
							Dunnett.GLM.zero.treatment.action = "log(x+1)",	# Dunnett.GLM method to be used for treatments only containing zeros
							show.progress = TRUE,			# show progress for each shift of lambda
							show.results = TRUE,			# show results
							Dunnett.GLM.family = "quasipoisson") {	# distribution family used for the GLM within Dunnett.GLM (alternatives: "nb", "poisson")

	Multi.group.test.bMDD(groups = groups,
						  counts = counts,
						  control.name = control.name,
						  alpha = alpha,
						  shift.step = shift.step,
						  bootstrap.runs = bootstrap.runs,
						  power = power,
						  max.iterations = max.iterations,
						  use.fixed.random.seed = use.fixed.random.seed,
						  Dunnett.GLM.zero.treatment.action = Dunnett.GLM.zero.treatment.action,
						  show.progress = show.progress,
						  show.results = show.results,
						  Dunnett.GLM.family = Dunnett.GLM.family,
						  test = "GLM.Dunnett")
}
