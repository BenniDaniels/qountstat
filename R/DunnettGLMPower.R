#' @title Dunnett.GLM power
#' @description The basic idea of Dunnett.GLM power calculations is to do parametric bootstrapping
#' for each dose/concentration group and to evaluate the proportion of results significantly different from the control.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param bootstrap.runs Number of bootstrap runs
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param Dunnett.GLM.zero.treatment.action GLM.Dunnett method to be used for treatments only containing zeros
#' @param show.progress Show progress for each shift of lambda
#' @param show.results Show results
#' @param Dunnett.GLM.family Distribution family to be used for the GLM within the Dunnett.GLM method
#' (use either "quasipoisson", "nb" or "poisson", see \code{\link{Dunnett.GLM}}). Note that this argument only
#' affects the model that is fitted to the resampled data; the bootstrap samples themselves are always
#' drawn from a Poisson distribution.
#' @return Data frame with results from power analysis
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' # Test Dunnett.GLM power
#' Dunnett.GLM.power(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		bootstrap.runs = 10,	# Caution: low number of bootstrap runs for testing
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE)
#'
#' \donttest{
#' # Test Dunnett.GLM power using a negative binomial distribution
#' Dunnett.GLM.power(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		bootstrap.runs = 10,	# Caution: low number of bootstrap runs for testing
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE,
#'		Dunnett.GLM.family = "nb")
#'
#' # Test Dunnett.GLM power using a classical Poisson distribution
#' Dunnett.GLM.power(groups = Daphnia.counts$Concentration,
#'		counts = Daphnia.counts$Number_Young,
#'		control.name = NULL,
#'		alpha = 0.05,
#'		bootstrap.runs = 10,	# Caution: low number of bootstrap runs for testing
#'		use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'		Dunnett.GLM.zero.treatment.action = "log(x+1)",
#'		show.progress = TRUE,
#'		show.results = TRUE,
#'		Dunnett.GLM.family = "poisson")
#' }
#' @export
Dunnett.GLM.power = function(groups,						# group vector
							 counts,						# vector with count data
							 control.name = NULL,			# character string with control group name
							 alpha = 0.05,					# significance level
							 bootstrap.runs = 200,			# number of bootstrap runs (draw Poisson data n times)
							 use.fixed.random.seed = NULL,  # fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
							 Dunnett.GLM.zero.treatment.action = "log(x+1)",	# GLM.Dunnett method to be used for treatments only containing zeros
							 show.progress = TRUE,			# show progress for each shift of lambda
							 show.results = TRUE,			# show results
							 Dunnett.GLM.family = "quasipoisson") {	# distribution family used for the GLM within Dunnett.GLM (alternatives: "nb", "poisson")

	Multi.group.test.power(groups = groups,
						   counts = counts,
						   control.name = control.name,
						   alpha = alpha,
						   bootstrap.runs = bootstrap.runs,
						   use.fixed.random.seed = use.fixed.random.seed,
						   Dunnett.GLM.zero.treatment.action = Dunnett.GLM.zero.treatment.action,
						   show.progress = show.progress,
						   show.results = show.results,
						   Dunnett.GLM.family = Dunnett.GLM.family,
						   test = "GLM.Dunnett")
}
