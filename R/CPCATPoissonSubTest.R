#' @title CPCAT Poisson sub-test
#' @description Helper function for CPCAT
#' @param dat Data frame to be evaluated
#' @param contrast Contrast matrix
#' @param bootstrap.runs Number of bootstrap runs
#' @return p-value for tested data and contrast
CPCAT.Poisson.sub.test = function(dat,							# data to be evaluated
								  contrast, 					# contrast matrix
								  bootstrap.runs = 10000) {		# number of bootstrap runs
	# Exclude the treatment groups that are not indicated by the contrast information
	datsheets = c(1, which(contrast == 1) + 1)    # Index to control and CONSIDERED treatments
	dat2 = list()                                 # Control and considered treatment data as a list

	# Populate dat2 with data from the specified treatments
	for (l in 1:length(datsheets)) {
		dat2[[l]] = dat[[datsheets[l]]]
	}

	# Initialize vectors for statistics
	musML = ns = xs = rep(0, length(dat2))

	# Calculate sample sizes, total abundance, and mean abundance for each group
	for (i in 1:length(musML)) {
		ns[i] = nrow(dat2[[i]])        # Sample size per group
		xs[i] = sum(dat2[[i]][, 1])    # Total abundance per group
		musML[i] = xs[i] / ns[i]       # Mean abundance per group
	}

	# Calculate 'total distance' etaML between control and all considered treatments
	etaML = sum((sqrt(musML[-1]) - sqrt(musML[1]))^2)
	n = sum(ns)                        # Total sample size over control and considered treatments
	x = sum(xs)                        # Total abundance over control and considered treatments
	mu0RML = x / n                     # Mean abundance over control and considered treatments

	# Return a p-value of 1 if control and all (considered) treatments are zero (i.e., no difference)
	if (mu0RML == 0) {
		return(1)
	}

	# Create artificial data (Poisson distributed)
	artificialdata = pseudomus = list()
	for (j in 1:length(dat2)) {
		artificialdata[[j]] = stats::rpois(bootstrap.runs, ns[j] * mu0RML)  # Artificial total abundance per group
		pseudomus[[j]] = artificialdata[[j]] / ns[j]    # Artificial mean abundance per group
	}

	# Calculate pseudo etas for the artificial data
	pseudoetasML = rep(0, bootstrap.runs)
	for (l in 1:bootstrap.runs) {
		pseudomushelp = vector()
		for (i in 1:length(pseudomus)) {
			pseudomushelp[i] = pseudomus[[i]][l]
		}
		pseudoetasML[l] = sum((sqrt(pseudomushelp[-1]) - sqrt(pseudomushelp[1]))^2)
	}

	# Calculate p-value from the number of artificial datasets causing a higher distance
	pvalue = length(which(pseudoetasML > etaML)) / bootstrap.runs

	return(pvalue)
}
