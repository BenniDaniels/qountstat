#' @title CPCAT power
#' @description The basic idea of CPCAT power calculations is to do parametric bootstrapping
#' for each dose/concentration group and to evaluate the proportion of results significantly different from the control.
#' @param groups Group vector
#' @param counts Vector with count data
#' @param control.name Character string with control group name (optional)
#' @param alpha Significance level
#' @param bootstrap.runs Number of bootstrap runs
#' @param use.fixed.random.seed Use fixed seed, e.g. 123, for reproducible results. If NULL no seed is set.
#' @param CPCAT.bootstrap.runs Bootstrap runs within CPCAT method
#' @param show.progress Show progress for each shift of lambda
#' @param show.results Show results
#' @return Data frame with results from power analysis
#' @examples
#' Daphnia.counts	# example data provided alongside the package
#'
#' # Test CPCAT power
#' CPCAT.power(groups = Daphnia.counts$Concentration,
#'			   counts = Daphnia.counts$Number_Young,
#'			   control.name = NULL,
#'			   alpha = 0.05,
#'			   bootstrap.runs = 10,		# Caution: low number of bootstrap runs for testing
#'			   use.fixed.random.seed = 123,  #fixed seed for reproducible results
#'			   CPCAT.bootstrap.runs = 10,# Caution: low number of bootstrap runs for testing
#'			   show.progress = TRUE,
#'			   show.results = TRUE)
#' @export
CPCAT.power = function(groups,						# group vector
					   counts,						# vector with count data
					   control.name = NULL,			# character string with control group name
					   alpha = 0.05,				# significance level
					   bootstrap.runs = 200,		# number of bootstrap runs (draw Poisson data n times)
					   use.fixed.random.seed = NULL,# fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
					   CPCAT.bootstrap.runs = 200,	# bootstrap runs within CPCAT method
					   show.progress = TRUE,		# show progress for each shift of lambda
					   show.results = TRUE) {		# show results

	Multi.group.test.power(groups = groups,
						   counts = counts,
						   control.name = control.name,
						   alpha = alpha,
						   bootstrap.runs = bootstrap.runs,
						   use.fixed.random.seed = use.fixed.random.seed,
						   CPCAT.bootstrap.runs = CPCAT.bootstrap.runs,
						   show.progress = show.progress,
						   show.results = show.results,
						   test = "CPCAT")
}
