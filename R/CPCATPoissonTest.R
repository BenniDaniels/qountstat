#' @title CPCAT Poisson test
#' @description Helper function for CPCAT
#' @param dat Data frame to be evaluated
#' @param contrastmatrix Contrast matrix
#' @param bootstrap.runs Number of bootstrap runs
#' @return List of p-values for tested contrast matrices
CPCAT.Poisson.test = function(dat, contrastmatrix, bootstrap.runs = 10000) {
	# Re-structure the input data as a list object
	dat2 = list()
	treatments = levels(dat$Groups)

	# Organize data by treatment
	for (j in 1:length(treatments)) {
		index = which(dat$Groups == treatments[j])
		dat2[[j]] = dat[index, ]
	}

	# Initialize p-values vector
	pvalues = rep(0, nrow(contrastmatrix))

	# Perform Poisson test for each hypothesis (each row in contrastmatrix)
	for (j in 1:nrow(contrastmatrix)) {
		pvalues[j] = CPCAT.Poisson.sub.test(dat = dat2,
											contrast = contrastmatrix[j, ],
											bootstrap.runs = bootstrap.runs)
	}

	# For each treatment, provide the maximum p-value from all treatment-related hypotheses
	return(list(contrastmatrix_pvalues = cbind(contrastmatrix, pvalues), maxpvalue = max(pvalues)))
}
