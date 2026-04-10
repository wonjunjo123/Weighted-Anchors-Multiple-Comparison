
# Need to account for relative path of user
source("/Users/wonjunjo/Library/CloudStorage/Box-Box/Academics/Fourth Year Academics/MATH310/MC Project/mc_project budget (main).R")


data("InsectSprays")

J <- nrow(InsectSprays) / nlevels(InsectSprays$spray)  # observations per group
I <- nlevels(InsectSprays$spray)                        # number of groups

means <- tapply(InsectSprays$count, InsectSprays$spray, mean)

MSE <- tapply(InsectSprays$count, InsectSprays$spray, function(x) sum((x - mean(x))^2))
MSE <- sum(MSE) / (I * (J - 1))

cat("Group Means (mean insect count per spray):\n")
print(round(means, 3))
cat("\nJ:", J, "\n")
cat("MSE:", round(MSE, 4), "\n")

eta1 <- calibrate_eta(I = I, J = J)
cat("\nCalibrated eta:", round(eta, 4), "\n")

results <- my_comparison(means, J, MSE, eta = eta1)
cat("\nPairwise Comparison Results:\n")
kable(results[order(results$diff), ])
