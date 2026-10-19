library(multcomp)

################################################################################
##################################   DATA   ####################################
################################################################################

# Example from Hothorn et al. (2020)
daphnia =
	structure(list('Concentration' = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56,
									   1.56, 1.56, 3.12, 3.12, 3.12, 3.12, 3.12, 3.12, 3.12, 3.12, 3.12, 3.12, 6.25,
									   6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 6.25, 12.5, 12.5, 12.5, 12.5,
									   12.5, 12.5, 12.5, 12.5, 12.5, 12.5, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25),
				   Adults = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9,
				   		   10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
				   Number_Young = c(27, 30, 29, 31, 16, 15, 18, 17, 14, 27, 32, 35, 32, 26, 18, 29, 27, 16, 35, 13, 39, 30,
				   				 33, 33, 36, 33, 33, 27, 38, 44, 27, 34, 36, 34, 31, 27, 33, 21, 33, 31, 10, 13, 7, 7, 7,
				   				 10, 10, 16, 12, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
			  class = "data.frame", row.names = c(NA, 60L))

daphnia$Conc = as.factor(daphnia$Concentration)

contingency.table <- matrix(c(20, 10, 40, 10, 90, 15), nrow=2, ncol=3)
rownames(contingency.table)=c("Success", "Failure")
colnames(contingency.table)=c("Control", "T1", "T2")
contingency.table


################################################################################
##################################   CPCAT   ###################################
################################################################################

# testing CPCAT
res = CPCAT(groups=daphnia$Conc,
			counts=daphnia$Number_Young,
			control.name = NULL,
			bootstrap.runs = 10000,
			use.fixed.random.seed = 123,
			get.contrasts.and.p.values = F,
			show.output = T)

################################################################################
###############################   DUNNETT GLM   ################################
################################################################################

# test Dunnett GLM
res = Dunnett.GLM(groups=daphnia$Conc,
				  counts=daphnia$Number_Young,
				  control.name = NULL,
				  zero.treatment.action = "identity.link")

res = Dunnett.GLM(groups=daphnia$Conc,
				  counts=daphnia$Number_Young,
				  control.name = NULL,
				  zero.treatment.action = "log(x+1)")

# test Dunnett GLM using a negative binomial distribution
res = Dunnett.GLM(groups=daphnia$Conc,
				  counts=daphnia$Number_Young,
				  control.name = NULL,
				  zero.treatment.action = "log(x+1)",
				  family = "nb")

##########################   family argument checks   ##########################

testthat::test_that("the default of Dunnett.GLM still is the Quasi-Poisson distribution", {
	res.default = Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
							  zero.treatment.action = "log(x+1)", show.output = FALSE)
	res.qp		= Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
							  zero.treatment.action = "log(x+1)", show.output = FALSE,
							  family = "quasipoisson")
	# note: the adjusted p-values of multcomp::glht are obtained by numerical integration, hence
	# they are only reproducible up to the accuracy of that algorithm (a tolerance is required and
	# the attributes of the p-value vector have to be dropped)
	testthat::expect_equal(as.numeric(res.default$Results$test$pvalues),
						   as.numeric(res.qp$Results$test$pvalues), tolerance = 1e-3)
	testthat::expect_equal(as.numeric(res.default$Results$test$coefficients),
						   as.numeric(res.qp$Results$test$coefficients))
	testthat::expect_true(any(grepl("Quasi-Poisson", res.default$Info[, 1], fixed = TRUE)))
})

testthat::test_that("Dunnett.GLM works with the negative binomial distribution", {
	# note: the data contain a treatment consisting only of zeros, so the identity link cannot be
	# combined with the negative binomial model here (see ?Dunnett.GLM)
	for (zta in c("log(x+1)")) {
		res.nb = Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
							 zero.treatment.action = zta, show.output = FALSE, family = "nb")
		# the model must be a negative binomial model with an estimated theta
		testthat::expect_s3_class(res.nb$Results$model, "negbin")
		testthat::expect_true(is.numeric(res.nb$Results$model$theta))
		testthat::expect_gt(res.nb$Results$model$theta, 0)
		# one p-value per treatment group is returned and all p-values are valid
		testthat::expect_length(res.nb$Results$test$pvalues, length(levels(daphnia$Conc)) - 1)
		testthat::expect_true(all(res.nb$Results$test$pvalues >= 0 & res.nb$Results$test$pvalues <= 1))
		# the distribution that was used has to be reported
		testthat::expect_true(any(grepl("negative binomial", res.nb$Info[, 1], fixed = TRUE)))
	}
})

testthat::test_that("Dunnett.GLM rejects invalid values of the family argument", {
	testthat::expect_error(Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
									   family = "gaussian"), "family")
	testthat::expect_error(Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
									   family = NULL), "family")
	testthat::expect_error(Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young,
									   family = c("nb", "quasipoisson")), "family")
	# the negative binomial model is a discrete model and requires integer counts
	testthat::expect_error(Dunnett.GLM(groups=daphnia$Conc, counts=daphnia$Number_Young + 0.5,
									   family = "nb"), "integer")
})

testthat::test_that("the family argument is passed on by Dunnett.GLM.power and Dunnett.GLM.bMDD", {
	testthat::expect_error(Dunnett.GLM.power(groups=daphnia$Conc, counts=daphnia$Number_Young,
											 bootstrap.runs = 2, show.progress = FALSE,
											 show.results = FALSE,
											 Dunnett.GLM.family = "gaussian"), "family")
	testthat::expect_error(Dunnett.GLM.bMDD(groups=daphnia$Conc, counts=daphnia$Number_Young,
											bootstrap.runs = 2, show.progress = FALSE,
											show.results = FALSE,
											Dunnett.GLM.family = "gaussian"), "family")

	res.power = Dunnett.GLM.power(groups=daphnia$Conc, counts=daphnia$Number_Young,
								  bootstrap.runs = 5, use.fixed.random.seed = 123,
								  Dunnett.GLM.zero.treatment.action = "log(x+1)",
								  show.progress = FALSE, show.results = FALSE,
								  Dunnett.GLM.family = "nb")
	testthat::expect_s3_class(res.power, "data.frame")
	testthat::expect_equal(nrow(res.power), length(levels(daphnia$Conc)) - 1)
	testthat::expect_true(all(res.power$Percent.power >= 0 & res.power$Percent.power <= 100))
})

################################################################################
##################################   CPFISH   ##################################
################################################################################


CPFISH(contingency.table = contingency.table, control.name = NULL, simulate.p.value = T, use.fixed.random.seed = 123)


################################################################################
#################################   Test bMDD   ################################
################################################################################

###############################   CPCAT tests   ################################
# Idea: shift lambda of Poisson distribution until there is a certain proportion of significant results
CPCAT.bMDD(groups = daphnia$Conc,
		   counts = daphnia$Number_Young,
		   control.name = NULL,
		   alpha = 0.05,
		   shift.step = -1,
		   bootstrap.runs = 5,
		   power = 0.8,
		   max.iterations = 1000,
		   use.fixed.random.seed = 123,
		   CPCAT.bootstrap.runs = 10,
		   show.progress = T,
		   show.results = T)

############################   GLM.Dunnett tests   #############################

res = Dunnett.GLM.bMDD(groups = daphnia$Conc,
					   counts = daphnia$Number_Young,
					   control.name = NULL,
					   alpha = 0.05,
					   shift.step = -1,
					   bootstrap.runs = 5,
					   power = 0.8,
					   max.iterations = 1000,
					   use.fixed.random.seed = 123,
					   Dunnett.GLM.zero.treatment.action = "log(x+1)",
					   show.progress = T,
					   show.results = T)

res = Dunnett.GLM.bMDD(groups = daphnia$Conc,
					   counts = daphnia$Number_Young,
					   control.name = NULL,
					   alpha = 0.05,
					   shift.step = -1,
					   bootstrap.runs = 5,
					   power = 0.8,
					   max.iterations = 1000,
					   use.fixed.random.seed = 123,
					   Dunnett.GLM.zero.treatment.action = "log(x+1)",
					   show.progress = T,
					   show.results = T,
					   Dunnett.GLM.family = "nb")

##############################   CPFISH tests   ################################

CPFISH.bMDD(contingency.table = contingency.table, # contingency.table is a matrix with observed data (e.g. survival counts)
			control.name = NULL,		# character string with control group name
			alpha = 0.05,				# significance level
			shift.step = -0.1,		# step of shift (negative as a reduction is assumed)
			bootstrap.runs = 10,		# number of bootstrap runs (draw Poisson data n times)
			power = 0.8,			# proportion of bootstrap.runs that return significant differences
			max.iterations = 1000,		# max number of iterations to not get stuck in the while loop
			simulate.p.value = TRUE,	# use simulated p-values or not
			use.fixed.random.seed = 123,# fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
			show.progress = T,			# show progress for each shift of lambda
			show.results = T) 			# show results


################################################################################
################################   Test power   ################################
################################################################################

###############################   CPCAT tests   ################################
# Idea: shift lambda of Poisson distribution until there is a certain proportion of significant results
CPCAT.power(groups = daphnia$Conc,
			counts = daphnia$Number_Young,
			control.name = NULL,
			alpha = 0.05,
			bootstrap.runs = 10,
			use.fixed.random.seed = 123,
			CPCAT.bootstrap.runs = 10,
			show.progress = T,
			show.results = T)


############################   GLM.Dunnett tests   #############################

res = Dunnett.GLM.power(groups = daphnia$Conc,
						counts = daphnia$Number_Young,
						control.name = NULL,
						alpha = 0.05,
						bootstrap.runs = 10,
						use.fixed.random.seed = 123,
						Dunnett.GLM.zero.treatment.action = "log(x+1)",
						show.progress = T,
						show.results = T)

res = Dunnett.GLM.power(groups = daphnia$Conc,
						counts = daphnia$Number_Young,
						control.name = NULL,
						alpha = 0.05,
						bootstrap.runs = 10,
						use.fixed.random.seed = 123,
						Dunnett.GLM.zero.treatment.action = "log(x+1)",
						show.progress = T,
						show.results = T,
						Dunnett.GLM.family = "nb")

##############################   CPFISH tests   ################################

CPFISH.power(contingency.table = contingency.table, # contingency.table is a matrix with observed data (e.g. survival counts)
			 control.name = NULL,		# character string with control group name
			 alpha = 0.05,				# significance level
			 bootstrap.runs = 10,		# number of bootstrap runs (draw Poisson data n times)
			 simulate.p.value = TRUE,	# use simulated p-values or not
			 use.fixed.random.seed = 123,# fix seed, e.g. 123, for random numbers if desired (enables to reproduce results)
			 show.progress = T,			# show progress for each shift of lambda
			 show.results = T) 			# show results


################################################################################
################################################################################
###   ASSERTIONS                                                             ###
###                                                                          ###
###   The calls above only check that the functions run without an error.    ###
###   The tests below additionally check that the returned values are        ###
###   correct, reproducible and that invalid input is rejected.              ###
################################################################################
################################################################################

n.groups = length(levels(daphnia$Conc))

###############################   CP.hypotheses   ##############################

testthat::test_that("CP.hypotheses returns the closure of all group combinations", {
	h = CP.hypotheses(3)
	testthat::expect_type(h, "list")
	# one hypothesis per treatment-versus-control comparison
	testthat::expect_length(h, 3)
	testthat::expect_named(h, c("H0: mu_0 = mu_1", "H0: mu_0 = mu_2", "H0: mu_0 = mu_3"))
	# every hypothesis matrix has one column per group and contains only 0/1 indicators
	for (mat in h) {
		testthat::expect_equal(ncol(mat), 3)
		testthat::expect_true(all(mat %in% c(0, 1)))
	}
	# treatment names are passed through if provided
	h.named = CP.hypotheses(2, treatment.names = c("Control", "Dose1"))
	testthat::expect_equal(ncol(h.named[[1]]), 2)
})

##################################   CPCAT   ###################################

testthat::test_that("CPCAT returns a valid result object", {
	res = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
				bootstrap.runs = 200, use.fixed.random.seed = 123, show.output = FALSE)

	testthat::expect_named(res, c("Results", "Info"))
	testthat::expect_s3_class(res$Results, "data.frame")
	testthat::expect_named(res$Results, c("Hypothesis", "p.values", "Signif."))
	# one hypothesis per treatment group (all groups except the control)
	testthat::expect_equal(nrow(res$Results), n.groups - 1)
	testthat::expect_true(all(res$Results$p.values >= 0 & res$Results$p.values <= 1))
	# the highest concentration consists of zeros only and has to be significant
	testthat::expect_lt(res$Results$p.values[n.groups - 1], 0.05)
	# NOEC/LOEC have to be reported
	testthat::expect_true(any(grepl("NOEC", res$Info[, 1], fixed = TRUE)))
})

testthat::test_that("CPCAT is reproducible with a fixed random seed", {
	res.a = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
				  bootstrap.runs = 200, use.fixed.random.seed = 123, show.output = FALSE)
	res.b = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
				  bootstrap.runs = 200, use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_identical(res.a$Results$p.values, res.b$Results$p.values)
})

testthat::test_that("CPCAT returns the contrasts if requested", {
	res = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
				bootstrap.runs = 100, use.fixed.random.seed = 123,
				get.contrasts.and.p.values = TRUE, show.output = FALSE)
	testthat::expect_named(res, c("Contrasts", "Results", "Info"))
	testthat::expect_s3_class(res$Contrasts, "data.frame")
	testthat::expect_gt(nrow(res$Contrasts), 0)
})

##################################   CPFISH   ##################################

testthat::test_that("CPFISH returns a valid result object", {
	res = CPFISH(contingency.table = contingency.table, control.name = NULL,
				 simulate.p.value = TRUE, use.fixed.random.seed = 123, show.output = FALSE)

	testthat::expect_named(res, c("Results", "Info"))
	testthat::expect_named(res$Results, c("Treatment", "p.values", "Signif."))
	# one row per treatment, i.e. all columns of the contingency table except the control
	testthat::expect_equal(nrow(res$Results), ncol(contingency.table) - 1)
	testthat::expect_equal(as.character(res$Results$Treatment), colnames(contingency.table)[-1])
	testthat::expect_true(all(res$Results$p.values >= 0 & res$Results$p.values <= 1))
	testthat::expect_true(any(grepl("Simulated p-values used.", res$Info[, 1], fixed = TRUE)))
})

testthat::test_that("CPFISH is reproducible and rejects an unknown control", {
	res.a = CPFISH(contingency.table = contingency.table, use.fixed.random.seed = 123, show.output = FALSE)
	res.b = CPFISH(contingency.table = contingency.table, use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_identical(res.a$Results$p.values, res.b$Results$p.values)

	testthat::expect_error(CPFISH(contingency.table = contingency.table,
								  control.name = "does.not.exist", show.output = FALSE),
						   "control cannot be found")
	testthat::expect_error(CPFISH(contingency.table = contingency.table,
								  control.name = 1, show.output = FALSE),
						   "character string")
})

################################   Dunnett.GLM   ###############################

testthat::test_that("Dunnett.GLM returns a valid result object", {
	res = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
					  zero.treatment.action = "log(x+1)", show.output = FALSE)

	testthat::expect_named(res, c("Results", "Info"))
	testthat::expect_length(res$Results$test$pvalues, n.groups - 1)
	testthat::expect_true(all(res$Results$test$pvalues >= 0 & res$Results$test$pvalues <= 1))
	# the highest concentration consists of zeros only and has to be significant
	testthat::expect_lt(res$Results$test$pvalues[n.groups - 1], 0.05)
})

testthat::test_that("Dunnett.GLM respects a user defined control group", {
	res = Dunnett.GLM(groups = as.character(daphnia$Conc), counts = daphnia$Number_Young,
					  control.name = "12.5", zero.treatment.action = "log(x+1)", show.output = FALSE)
	# every comparison has to be made against the specified control
	testthat::expect_true(all(grepl("12.5$", names(res$Results$test$coefficients))))
	testthat::expect_length(res$Results$test$pvalues, n.groups - 1)

	testthat::expect_error(Dunnett.GLM(groups = as.character(daphnia$Conc), counts = daphnia$Number_Young,
									   control.name = "does.not.exist"), "control cannot be found")
})

##########################   input validation (all)   ##########################

testthat::test_that("negative counts are rejected", {
	# note: this used to be checked with min() instead of any(), which only triggered if ALL
	# values were negative, i.e. a single negative count passed the check unnoticed
	counts.neg = daphnia$Number_Young
	counts.neg[1] = -5

	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = counts.neg,
								 bootstrap.runs = 10, show.output = FALSE), "non-negative")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = counts.neg,
									   show.output = FALSE), "non-negative")
	testthat::expect_error(CPCAT.power(groups = daphnia$Conc, counts = counts.neg,
									   bootstrap.runs = 2, show.progress = FALSE,
									   show.results = FALSE), "non-negative")
	testthat::expect_error(CPCAT.bMDD(groups = daphnia$Conc, counts = counts.neg,
									  bootstrap.runs = 2, show.progress = FALSE,
									  show.results = FALSE), "non-negative")

	table.neg = contingency.table
	table.neg[1, 1] = -1
	testthat::expect_error(CPFISH.power(contingency.table = table.neg, bootstrap.runs = 2,
										show.progress = FALSE, show.results = FALSE), "non-negative")
	testthat::expect_error(CPFISH.bMDD(contingency.table = table.neg, bootstrap.runs = 2,
									   show.progress = FALSE, show.results = FALSE), "non-negative")
})

testthat::test_that("rows containing NA values are excluded instead of causing an error", {
	# note: the check for negative values used to fail with "missing value where TRUE/FALSE
	# needed" as soon as the counts contained an NA
	counts.na = daphnia$Number_Young
	counts.na[c(3, 17)] = NA

	res.cpcat = CPCAT(groups = daphnia$Conc, counts = counts.na, bootstrap.runs = 200,
					  use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_true(any(grepl("2 rows with NA values were excluded!",
									res.cpcat$Info[, 1], fixed = TRUE)))
	testthat::expect_equal(nrow(res.cpcat$Results), n.groups - 1)

	res.glm = Dunnett.GLM(groups = daphnia$Conc, counts = counts.na,
						  zero.treatment.action = "log(x+1)", show.output = FALSE)
	testthat::expect_true(any(grepl("2 rows with NA values were excluded!",
									res.glm$Info[, 1], fixed = TRUE)))

	# the power/bMDD functions return a plain data frame, hence they report this as a message
	testthat::expect_message(CPCAT.power(groups = daphnia$Conc, counts = counts.na,
										 bootstrap.runs = 2, CPCAT.bootstrap.runs = 10,
										 use.fixed.random.seed = 123, show.progress = FALSE,
										 show.results = FALSE),
							 "2 rows with NA values were excluded!")
})

testthat::test_that("inconsistent or insufficient input is rejected", {
	testthat::expect_error(CPCAT(groups = daphnia$Conc[-1], counts = daphnia$Number_Young,
								 show.output = FALSE), "don't match")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc[-1], counts = daphnia$Number_Young,
									   show.output = FALSE), "don't match")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = as.character(daphnia$Number_Young),
									   show.output = FALSE), "non-negative")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
									   zero.treatment.action = "log(x)"), "zero.treatment.action")
	testthat::expect_error(CPCAT.power(groups = rep("A", 10), counts = 1:10, bootstrap.runs = 2,
									   show.progress = FALSE, show.results = FALSE), "Too few groups")
	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
								 use.fixed.random.seed = "abc", show.output = FALSE),
						   "numeric value or NULL")
})

##############################   power and bMDD   ##############################

testthat::test_that("the power functions return valid percentages", {
	res.cpcat = CPCAT.power(groups = daphnia$Conc, counts = daphnia$Number_Young,
							bootstrap.runs = 5, CPCAT.bootstrap.runs = 20,
							use.fixed.random.seed = 123, show.progress = FALSE, show.results = FALSE)
	testthat::expect_s3_class(res.cpcat, "data.frame")
	testthat::expect_named(res.cpcat, c("Hypothesis", "Percent.power"))
	testthat::expect_equal(nrow(res.cpcat), n.groups - 1)
	testthat::expect_true(all(res.cpcat$Percent.power >= 0 & res.cpcat$Percent.power <= 100))

	res.fish = CPFISH.power(contingency.table = contingency.table, bootstrap.runs = 10,
							use.fixed.random.seed = 123, show.progress = FALSE, show.results = FALSE)
	testthat::expect_s3_class(res.fish, "data.frame")
	testthat::expect_equal(nrow(res.fish), ncol(contingency.table) - 1)
})

testthat::test_that("the bMDD functions return non-negative minimum detectable differences", {
	res.cpcat = CPCAT.bMDD(groups = daphnia$Conc, counts = daphnia$Number_Young,
						   shift.step = -2, bootstrap.runs = 5, CPCAT.bootstrap.runs = 20,
						   use.fixed.random.seed = 123, show.progress = FALSE, show.results = FALSE)
	testthat::expect_s3_class(res.cpcat, "data.frame")
	testthat::expect_named(res.cpcat, c("Hypothesis", "bMDD", "bMDD.percent"))
	testthat::expect_equal(nrow(res.cpcat), n.groups - 1)
	# a bMDD is either not reached (NA) or a non-negative difference
	valid = is.na(res.cpcat$bMDD) | res.cpcat$bMDD >= 0
	testthat::expect_true(all(valid))
	# bMDD% has to be consistent with bMDD relative to the control mean
	control.mean = mean(daphnia$Number_Young[daphnia$Conc == levels(daphnia$Conc)[1]])
	testthat::expect_equal(res.cpcat$bMDD.percent, res.cpcat$bMDD / control.mean * 100)
})

testthat::test_that("the power and bMDD functions are reproducible with a fixed random seed", {
	args = list(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 5,
				use.fixed.random.seed = 123, Dunnett.GLM.zero.treatment.action = "log(x+1)",
				show.progress = FALSE, show.results = FALSE)
	testthat::expect_identical(do.call(Dunnett.GLM.power, args)$Percent.power,
							   do.call(Dunnett.GLM.power, args)$Percent.power)
})


################################################################################
##########################   OVERDISPERSION PRE-TEST   #########################
################################################################################

testthat::test_that("Overdispersion.test returns a valid result object", {
	res = Overdispersion.test(groups = daphnia$Conc, counts = daphnia$Number_Young,
							  bootstrap.runs = 300, use.fixed.random.seed = 123, show.output = FALSE)

	testthat::expect_named(res, c("Results", "Group.statistics", "Info", "null.distribution"))
	testthat::expect_s3_class(res$Results, "data.frame")
	testthat::expect_s3_class(res$Group.statistics, "data.frame")
	testthat::expect_length(res$null.distribution, 300)

	values = stats::setNames(res$Results$Value, res$Results$Parameter)
	testthat::expect_true(as.numeric(values[["Observed T"]]) > 0)
	p.value = as.numeric(values[["Monte Carlo p-value"]])
	testthat::expect_true(p.value > 0 & p.value <= 1)
	testthat::expect_true(values[["Decision"]] %in%
						  c("Overdispersion detected", "No overdispersion detected"))

	# the all-zero treatment carries no information on dispersion and is excluded by default
	testthat::expect_false("25" %in% res$Group.statistics$Group)
	testthat::expect_equal(nrow(res$Group.statistics), length(levels(daphnia$Conc)) - 1)
})

testthat::test_that("the reported dispersion ratio is the bias-corrected statistic", {
	res = Overdispersion.test(groups = daphnia$Conc, counts = daphnia$Number_Young,
							  bootstrap.runs = 300, use.fixed.random.seed = 123, show.output = FALSE)
	values = stats::setNames(res$Results$Value, res$Results$Parameter)

	n.observations = as.numeric(values[["Observations"]])
	n.groups = as.numeric(values[["Groups"]])
	correction = n.observations / (n.observations - n.groups)

	# phi = T * N/(N-G), the same relation must hold for the critical value
	testthat::expect_equal(as.numeric(values[["Dispersion ratio (phi)"]]),
						   as.numeric(values[["Observed T"]]) * correction, tolerance = 1e-3)
	testthat::expect_equal(as.numeric(values[["Critical dispersion ratio"]]),
						   as.numeric(values[["c_alpha"]]) * correction, tolerance = 1e-3)

	# the decision must be identical in both scales, since phi is a monotone transformation of T
	decision.T = as.numeric(values[["Observed T"]]) > as.numeric(values[["c_alpha"]])
	decision.phi = as.numeric(values[["Dispersion ratio (phi)"]]) > as.numeric(values[["Critical dispersion ratio"]])
	testthat::expect_equal(decision.T, decision.phi)
})

testthat::test_that("Overdispersion.test is reproducible and reacts to the reference distribution", {
	args = list(groups = daphnia$Conc, counts = daphnia$Number_Young,
				bootstrap.runs = 300, use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_identical(do.call(Overdispersion.test, args)$null.distribution,
							   do.call(Overdispersion.test, args)$null.distribution)

	# a tolerated variance-to-mean ratio has to shift the critical value upwards
	strict = do.call(Overdispersion.test, c(args, list(pretest.mode = "strict")))
	tolerant = do.call(Overdispersion.test, c(args, list(pretest.mode = "tolerant", tolerance.factor = 2)))
	c.strict = as.numeric(stats::setNames(strict$Results$Value, strict$Results$Parameter)[["c_alpha"]])
	c.tolerant = as.numeric(stats::setNames(tolerant$Results$Value, tolerant$Results$Parameter)[["c_alpha"]])
	testthat::expect_gt(c.tolerant, c.strict)
})

testthat::test_that("Overdispersion.test detects clear overdispersion and accepts Poisson data", {
	set.seed(42)
	groups = rep(paste0("G", 1:4), each = 8)
	poisson.counts = stats::rpois(32, lambda = 20)
	# negative binomial counts with a variance-to-mean ratio of 6
	nb.counts = stats::rnbinom(32, mu = 20, size = 20 / 5)

	res.poisson = Overdispersion.test(groups, poisson.counts, bootstrap.runs = 500,
									  pretest.mode = "strict", use.fixed.random.seed = 1,
									  show.output = FALSE)
	res.nb = Overdispersion.test(groups, nb.counts, bootstrap.runs = 500,
								 pretest.mode = "strict", use.fixed.random.seed = 1,
								 show.output = FALSE)

	decision = function(x) stats::setNames(x$Results$Value, x$Results$Parameter)[["Decision"]]
	testthat::expect_equal(decision(res.poisson), "No overdispersion detected")
	testthat::expect_equal(decision(res.nb), "Overdispersion detected")
})

testthat::test_that("Overdispersion.test rejects invalid input", {
	testthat::expect_error(Overdispersion.test(daphnia$Conc, daphnia$Number_Young,
											   pretest.mode = "mild"), "pretest.mode")
	testthat::expect_error(Overdispersion.test(daphnia$Conc, daphnia$Number_Young,
											   tolerance.factor = 1), "tolerance.factor")
	testthat::expect_error(Overdispersion.test(daphnia$Conc, daphnia$Number_Young,
											   bootstrap.runs = 10), "bootstrap.runs")
	testthat::expect_error(Overdispersion.test(daphnia$Conc, daphnia$Number_Young,
											   alpha = 1), "alpha")
	testthat::expect_error(Overdispersion.test(daphnia$Conc[-1], daphnia$Number_Young), "don't match")
})

########################   integration in CPCAT / Dunnett.GLM   ################

testthat::test_that("the Hampel identifier is switched off by default", {
	res.off = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 200,
					use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_false(any(grepl("Hampel", res.off$Info[, 1], fixed = TRUE)))

	# providing a threshold enables it as before
	res.on = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 200,
				   hampel.threshold = 5, use.fixed.random.seed = 123, show.output = FALSE)
	testthat::expect_true(any(grepl("Hampel", res.on$Info[, 1], fixed = TRUE)))

	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
								 bootstrap.runs = 200, hampel.threshold = "five"), "hampel.threshold")
})

testthat::test_that("enabling the overdispersion pre-test does not change the results", {
	# the pre-test draws random numbers, but must restore the state of the random number generator
	res.without = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 500,
						use.fixed.random.seed = 123, show.output = FALSE)
	res.with = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 500,
					 use.fixed.random.seed = 123, show.output = FALSE,
					 overdispersion.test = list(bootstrap.runs = 200, use.fixed.random.seed = 99))
	testthat::expect_identical(res.without$Results$p.values, res.with$Results$p.values)
	testthat::expect_true(any(grepl("Overdispersion pre-test", res.with$Info[, 1], fixed = TRUE)))
	testthat::expect_false(any(grepl("Overdispersion pre-test", res.without$Info[, 1], fixed = TRUE)))

	res.glm.without = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
								  zero.treatment.action = "log(x+1)", show.output = FALSE)
	res.glm.with = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
							   zero.treatment.action = "log(x+1)", show.output = FALSE,
							   overdispersion.test = list(bootstrap.runs = 200, use.fixed.random.seed = 99))
	testthat::expect_equal(as.numeric(res.glm.without$Results$test$coefficients),
						   as.numeric(res.glm.with$Results$test$coefficients))
	testthat::expect_true(any(grepl("Overdispersion pre-test", res.glm.with$Info[, 1], fixed = TRUE)))
})

testthat::test_that("the overdispersion.test argument is validated", {
	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
								 bootstrap.runs = 200, overdispersion.test = "yes"), "overdispersion.test")
	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
								 bootstrap.runs = 200,
								 overdispersion.test = list(unknown.setting = 1)), "Unknown setting")
	testthat::expect_error(CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young,
								 bootstrap.runs = 200, overdispersion.test = list(5)), "named")

	# FALSE and NULL both switch the pre-test off
	res = CPCAT(groups = daphnia$Conc, counts = daphnia$Number_Young, bootstrap.runs = 200,
				use.fixed.random.seed = 123, show.output = FALSE, overdispersion.test = FALSE)
	testthat::expect_false(any(grepl("Overdispersion pre-test", res$Info[, 1], fixed = TRUE)))
})


############################   Poisson distribution   #########################

testthat::test_that("Dunnett.GLM works with a classical Poisson distribution", {
	for (zta in c("identity.link", "log(x+1)")) {
		res = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
						  zero.treatment.action = zta, show.output = FALSE, family = "poisson")
		testthat::expect_equal(res$Results$model$family$family, "poisson")
		testthat::expect_length(res$Results$test$pvalues, length(levels(daphnia$Conc)) - 1)
		testthat::expect_true(all(res$Results$test$pvalues >= 0 & res$Results$test$pvalues <= 1))
		# the dispersion is not estimated but fixed by the model assumption
		testthat::expect_true(any(grepl("dispersion parameter fixed at 1", res$Info[, 1], fixed = TRUE)))
	}
})

testthat::test_that("the Poisson model yields the Quasi-Poisson estimates with rescaled errors", {
	res.poisson = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
							  zero.treatment.action = "log(x+1)", show.output = FALSE, family = "poisson")
	res.quasi = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
							zero.treatment.action = "log(x+1)", show.output = FALSE, family = "quasipoisson")

	# both models share the same mean structure, hence identical coefficients
	testthat::expect_equal(as.numeric(res.poisson$Results$test$coefficients),
						   as.numeric(res.quasi$Results$test$coefficients))

	# the Quasi-Poisson standard errors are the Poisson ones scaled by sqrt(dispersion)
	dispersion = summary(res.quasi$Results$model)$dispersion
	testthat::expect_equal(as.numeric(res.quasi$Results$test$sigma),
						   as.numeric(res.poisson$Results$test$sigma) * sqrt(dispersion),
						   tolerance = 1e-6)
})

testthat::test_that("the Poisson option requires integer counts and a valid family name", {
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young + 0.5,
									   family = "poisson"), "integer")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
									   family = "Poisson"), "family")
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
									   family = "gaussian"), "family")
})

testthat::test_that("the Poisson option is passed on by Dunnett.GLM.power and Dunnett.GLM.bMDD", {
	res.power = Dunnett.GLM.power(groups = daphnia$Conc, counts = daphnia$Number_Young,
								  bootstrap.runs = 5, use.fixed.random.seed = 123,
								  Dunnett.GLM.zero.treatment.action = "log(x+1)",
								  show.progress = FALSE, show.results = FALSE,
								  Dunnett.GLM.family = "poisson")
	testthat::expect_s3_class(res.power, "data.frame")
	testthat::expect_equal(nrow(res.power), length(levels(daphnia$Conc)) - 1)
	testthat::expect_true(all(res.power$Percent.power >= 0 & res.power$Percent.power <= 100))

	res.bMDD = Dunnett.GLM.bMDD(groups = daphnia$Conc, counts = daphnia$Number_Young,
								shift.step = -2, bootstrap.runs = 5, use.fixed.random.seed = 123,
								Dunnett.GLM.zero.treatment.action = "log(x+1)",
								show.progress = FALSE, show.results = FALSE,
								Dunnett.GLM.family = "poisson")
	testthat::expect_s3_class(res.bMDD, "data.frame")
	testthat::expect_true(all(is.na(res.bMDD$bMDD) | res.bMDD$bMDD >= 0))

	testthat::expect_error(Dunnett.GLM.power(groups = daphnia$Conc, counts = daphnia$Number_Young,
											 bootstrap.runs = 2, show.progress = FALSE,
											 show.results = FALSE,
											 Dunnett.GLM.family = "Poisson"), "Dunnett.GLM.family")
	testthat::expect_error(Dunnett.GLM.bMDD(groups = daphnia$Conc, counts = daphnia$Number_Young,
											bootstrap.runs = 2, show.progress = FALSE,
											show.results = FALSE,
											Dunnett.GLM.family = "gaussian"), "Dunnett.GLM.family")
})

testthat::test_that("nb is rejected with an identity link if a treatment contains only zeros", {
	# the fitted mean of an all-zero treatment is zero under an identity link, where the negative
	# binomial dispersion is not identifiable; the combination must fail with a clear message
	# rather than depending on floating point rounding
	testthat::expect_error(Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
									   zero.treatment.action = "identity.link", show.output = FALSE,
									   family = "nb"),
						   "not identifiable")

	# the other families keep working with the identity link
	for (fam in c("quasipoisson", "poisson")) {
		res = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
						  zero.treatment.action = "identity.link", show.output = FALSE, family = fam)
		testthat::expect_equal(res$Results$model$family$link, "identity")
	}

	# and nb works with the recommended setting
	res.nb = Dunnett.GLM(groups = daphnia$Conc, counts = daphnia$Number_Young,
						 zero.treatment.action = "log(x+1)", show.output = FALSE, family = "nb")
	testthat::expect_s3_class(res.nb$Results$model, "negbin")

	# without an all-zero treatment the log link is used anyway, so the setting is irrelevant
	keep = daphnia$Conc != "25"
	res.sub = Dunnett.GLM(groups = droplevels(daphnia$Conc[keep]), counts = daphnia$Number_Young[keep],
						  zero.treatment.action = "identity.link", show.output = FALSE, family = "nb")
	testthat::expect_s3_class(res.sub$Results$model, "negbin")
	testthat::expect_equal(res.sub$Results$model$family$link, "log")
})
