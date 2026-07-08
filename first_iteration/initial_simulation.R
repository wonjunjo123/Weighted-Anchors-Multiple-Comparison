
# Need to account for relative path for user
source("/Users/wonjunjo/Library/CloudStorage/Box-Box/Academics/Fourth Year Academics/MATH310/MC Project/mc_project budget (main).R")

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
