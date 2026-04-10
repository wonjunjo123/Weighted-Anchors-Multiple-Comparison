
# HEADER BLOCK #
# Method Name: Jo's Weighted Anchors Method
# Name: Wonjun Jo
# Date: 4/10/2026
# Brief Description: This is a multiple comparison method that begins with a baseline Sidak anchor. Then we distribute alpha budget to uncertain p-values by taking away from more certain p-values. The closer a p-value is to  the threshold, the more uncertain we are whether we can reject that comparison or not.

# Bugs Fixed
# We now correct the constraint equation to use the product of (1- alpha) instead of the product of alphas as before
# We use the correct degrees of freedom in our method
# Now we report Standard Error (SE) with the FWER simulations

# Implemented Suggestions
# We reconsidered the weighting direction so that we give more budget to uncertain p-values rather than to the comparisons with large observed differences.
# We implement an empirical calibration method for eta, much like the hyperparameter tuning in a machine learning algorithm
# In order to keep the method simple, we save the "step-down" structure and incorporating correlation structure for a future reiteration of the method



#install.packages("progress")
#install.packages("knitr")
#install.packages("tidyr")
#install.packages("ggplot2")
library(progress)
library(knitr)
library(tidyr)
library(ggplot2)

# variables for oil filter demo
oil_means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J <- 9
MSE <- 0.088
global_eta <- 0.001



set.seed(123)

my_comparison <- function(means, J, MSE, alpha = 0.05, eta = 1) {
	I <- length(means)
	df <- I*(J-1)
	m <- choose(I,2)
	epsilon <- 0.0001 # variable to control for various situations where limited floating point precision causes NaN
	
	pairs <- combn(I,2) # we use combination function to make code more expressive
	results <- data.frame( # we compute differences of means
		group1 = pairs[1,],
		group2 = pairs[2,],
		diff = abs(means[pairs[1,]] - means[pairs[2,]])
	)
	
	se <- sqrt(2*MSE/J)
	results$pval <- 2*(1 - pt(results$diff / se, df))
	
	anchor <- 1 - (1 - alpha)^(1/m)
	
	# uncertain comparisons (p near 0.05) get more budget
	# already-obvious comparisons (p near 0 or 1) get less
	results$dist <- prop_dist(results$pval, anchor) # we compute how "far" pvalue is to our significance level
	results$uncertainty <- 1/results$dist # the closer a value is to signficance level, the more uncertain we are
	total_uncertainty <- sum(results$uncertainty)
	
	# not likely to happen but avoid dividing by zero
	if (total_uncertainty == 0) {
		results$weights <- rep(1/m, m)  # if somehow we get no uncertainty (basically not going to happen) fall back to equal weights
	} else {
		# the weights represent how much of the "budget" you receive
		results$weights <- results$uncertainty / total_uncertainty
	}
	
	# this is the Sidak anchor, the default threshold for controlling FWER
	
	results <- results[order(results$weights),]
	
	results$threshold <- ifelse(
		results$pval > anchor, # For p-values that we still might need to reject, but are in the top 50% of the uncertain values that are greater than the anchor
		anchor + results$weights*eta, # we allocate budget to more uncertain p-values
		anchor - min(1/results$weights, anchor/eta - epsilon)*eta # for already significant ones, we take away more by distance
			# but we need to make sure that threshold > 0, hence the minimum
			# and we also use 1/weights because we want to steal budget from the larger distance p-values which have smaller weights (it's a computational trick)
			# We need some small but finite epsilon to account for division by zero and floating point memory of computers.
	)

	results$significant <- results$pval < results$threshold
	return(results)

}

# returns how far away you are but by a proportional log scale
# we do log because we want to "tame" the exponential differences in pvalues
prop_dist <- function(x,y) {
	abs(log(y/x))
}


# You only need to calibrate_eta for the specific I and J dimension size. Once calibrated, the eta will generalize to any other dataset of the same size.
# Doing this for each scenario is not overfitting the model to the training data; We are calling calibrate_eta 5 times because there are 5 different scenarios all with different data dimensions. If each of the scenarios were of the same I and J, we only need to calibrate once.
# this is a binary search algorithm
calibrate_eta <- function(I, J, target_fwer = 0.05, n_sim = 10000, tolerance = 0.001, precision = 0.001) {
	eta_low <- 0
	eta_high <- 1
	
	while (eta_high - eta_low > precision) {
		eta_mid <- (eta_low + eta_high) / 2 # our guess
		fwer <- simulate_fwer(I, J, n_sim=n_sim, ETA=eta_mid)
		cat("FWER: ", fwer, "\n")
		if (fwer < target_fwer & fwer > target_fwer - tolerance) { # if we get good enough, stop training
			break
		}
		
		if (fwer > target_fwer) {
			eta_high <- eta_mid # Too liberal , reduce eta
		} else {
			eta_low <- eta_mid # Too conservative , increase eta
		}
	}
	
	return (eta_mid)
}

# we simulate to validate whether or not our method controls FWER
# this is code from Professor Prince-Nelson, so not much commenting
simulate_fwer <- function(I, J, sigma = 1, n_sim = 10000, alpha = 0.05, ETA=1) {
	false_rejections <- 0
	
	for (sim in 1:n_sim) {
		data <- matrix(rnorm(I*J, mean = 0, sd = sigma), nrow = I, ncol = J)
		
		means <- rowMeans(data)
		#MSE <- mean(apply(data, 1, var)) # simplified
		MSE <- sum(apply(data, 1, function(row) sum((row - mean(row))^2))) / (I * (J - 1)) # pooled
		
		results <- my_comparison(means, J, MSE, eta = ETA)
		
		if (any(results$significant)) {
			false_rejections <- false_rejections + 1
		}
	}
	return(false_rejections / n_sim)	
}

# this is code from Professor Prince-Nelson, so not much commenting
simulate_power <- function(func, d, I = 5, J = 10, sigma = 1, n_sim = 1000, alpha = 0.05) {
	correct_rejections <- 0
	
	for (sim in 1:n_sim) {
		data <- rbind(
			matrix(rnorm(3*J, mean = 0, sd = sigma), nrow = 3, ncol = J),
			matrix(rnorm(2*J, mean = d, sd = sigma), nrow = 2, ncol = J))
		
		means <- rowMeans(data)
		#MSE <- mean(apply(data, 1, var))
		MSE <- sum(apply(data, 1, function(row) sum((row - mean(row))^2))) / (I * (J - 1))
		
		results <- func(means, J, MSE)
		
		if (results[results$group1 == 1 & results$group2 == 4, ]$significant) {
			correct_rejections <- correct_rejections + 1
		}
	}
	
	return(correct_rejections / n_sim)
	
}

# this is code from Professor Prince-Nelson, so not much commenting
tukey_hsd <- function(means, J, MSE, alpha = 0.05, eta=1) {
	I <- length(means)
	df_error <- I * (J-1)
	q_crit <- qtukey(1-alpha, I, df_error)
	hsd <- q_crit * sqrt(MSE / J)
	comps <- list()
	for (i in 1:(I-1)) {
		for (j in (i+1):I) {
			key <- paste0(i, "_", j)
			comps[[key]] <- list(
				i = i, j = j, diff = abs(means[i] - means[j]))
		}
	}
	for (key in names(comps)) {
		comps[[key]]$significant <- comps[[key]]$diff > hsd
	}
	
	result <- data.frame(
		group1 = sapply(comps, function(x) x$i),
		group2 = sapply(comps, function(x) x$j),
		diff = sapply(comps, function(x) x$diff),
		significant = sapply(comps, function(x) x$significant),
		row.names = NULL
	)
	return(result)
	#return(result[order(result$diff),])
}


cat("-------------------------------------------------")
cat("\nPart 2: Demonstration of Example 11.5 Oil Filters Data\n")

cat("\nSample Means: ")
cat(oil_means, "\n")
cat("J: ")
cat(J, "\n")
cat("MSE: ")
cat(MSE, "\n")

#eta1 <- calibrate_eta(length(oil_means), J)
eta1 <- global_eta
demonstration <- my_comparison(oil_means, J, MSE, eta=eta1)
kable(demonstration[order(demonstration$diff),])


cat("-------------------------------------------------")
cat("\nPart 3: Validation via Simulation\n")

scenarios <- list(
  	A = list(name = "A",  I = 3,  J = 10),
  	B = list(name = "B", I = 5,  J = 10),
  	C = list(name = "C", I = 7,  J = 10),
  	D = list(name = "D", I = 5,  J = 5),
  	E = list(name = "E", I = 5,  J = 20)
)

pb1 <- progress_bar$new(total = length(scenarios)) # just to display progress bar while doing simulations
for (key in names(scenarios)) {
	pb1$tick()
	eta1 <- calibrate_eta(scenarios[[key]]$I, scenarios[[key]]$J)
  	scenarios[[key]]$fwer <- simulate_fwer(scenarios[[key]]$I, scenarios[[key]]$J, ETA=eta1)
}

simulation_results <- data.frame(
		Scenario = sapply(scenarios, function(x) x$name),
		I = sapply(scenarios, function(x) x$I),
		J = sapply(scenarios, function(x) x$J),
		FWER = sapply(scenarios, function(x) x$fwer),
		SE = sapply(scenarios, function(x) sqrt(x$fwer * (1 - x$fwer) / 10000)),
		row.names = NULL
	)

kable(simulation_results)


cat("-------------------------------------------------")
cat("\nPart 4: Power Analysis\n")

deltas <- list(
	A = list(d = 0.5),
	B = list(d = 1.0),
	C = list(d = 1.5),
	D = list(d = 2.0)
)

pb2 <- progress_bar$new(total = length(deltas))
for (key in names(deltas)) {
	pb2$tick()
	deltas[[key]]$powerWA <- simulate_power(my_comparison, deltas[[key]]$d)
	deltas[[key]]$powerTukey <- simulate_power(tukey_hsd, deltas[[key]]$d)
}

analysis_results <- data.frame(
	delta = sapply(deltas, function(x) x$d),
	wa_power = sapply(deltas, function(x) x$powerWA),
	tukey_power = sapply(deltas, function(x) x$powerTukey),
	row.names = NULL
)

kable(analysis_results)

analysis_results %>%
  	pivot_longer(cols = c(wa_power, tukey_power), 
  		names_to = "method", 
  		values_to = "power") %>%
	ggplot(aes(x = delta, y = power, color = method)) +
	geom_line(linewidth=1) +
  	geom_point() +
  	labs(title="Tukey vs. Weighted Anchor Power Analysis") +
  	theme_minimal() +
  	theme(legend.position="bottom", plot.title = element_text(hjust=0.5))

