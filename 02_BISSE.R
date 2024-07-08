# Load necessary libraries
library(ape)
library(diversitree)
library(evobiR)
library(plyr)
library(geiger)


tree <- read.tree("Dys_tr_names.tre")
data <- read.table("cl3_completed.txt", header = TRUE, sep = "\t")
rownames(data) <- data$Comp1

# select species with defined cheliceral morphologies
clades_of_interest <- c("A", "B", "G", "D", "C", "E", "F", "I")
data <- data[data$group %in% clades_of_interest,]

# Check for discrepancies between the tree and the data
obj <- name.check(tree, data)
tree <- drop.tip(tree, obj$tree_not_data)
obj <- name.check(tree, data)
print(obj)

# Reorder the data to match the tree structure
y <- ReorderData(tree, data, taxa.names = "row names")

# Revalue the groups to binary states for BiSSE analysis between generalist and specialist
y$group <- revalue(y$group, c("D" = 0, "E" = 0, "A" = 0, "G" = 1, "B" = 1, "F" = 1, "I" = 1, "C" = 1))

# Prepare the tree for BiSSE analysis
new_tree <- tree
data_values <- as.numeric(y[, 4])
names(data_values) <- row.names(y)
new_tree$tip.state <- data_values

# Create the BiSSE likelihood function
equation <- make.bisse(tree = new_tree, states = new_tree$tip.state)
starting_params <- starting.point.bisse(new_tree)
mle_fit <- find.mle(equation, starting_params)

# Define and fit the null BiSSE models 
bissenull_model <- constrain(equation, lambda1 ~ lambda0, mu1 ~ mu0)
bissenull_mle <- find.mle(bissenull_model, starting_params[-c(2, 4)])
print(coef(bissenull_mle))

# Perform ANOVA to compare the full and null BiSSE models
bisseAnova_GS <- anova(mle_fit, null = bissenull_mle)
print(bisseAnova_GS)


# Run MCMC for the BiSSE model and calculate posterior probabilities
sums <- matrix(ncol = 1, nrow = 1000)
for (i in seq_along(sums)) {
  prior <- make.prior.exponential(1 / (2 * 0.4))
  w <- 1
  
  mcmc_bisse <- mcmc(equation, mle_fit$par, nsteps = 1000, prior = prior, w = w, print.every = 100)
  w <- diff(sapply(mcmc_bisse[,-c(1, ncol(mcmc_bisse))], quantile, c(0.05, 0.95)))
  mcmc_bisse <- mcmc(equation, mle_fit$par, nsteps = 1000, prior = prior, w = w, print.every = 100)
  
  sums[i, 1] <- sum(mcmc_bisse$lambda0 > mcmc_bisse$lambda1) / length(mcmc_bisse$lambda0)
}


colors <- setNames(c("blue", "red"), c("Generalist", "Specialist"))
profiles.plot(mcmc_bisse[, c("lambda0", "lambda1")], col.line = colors, las = 1, bty = "n",
              xlab = expression(lambda), cex.axis = 0.7)
legend("topright", names(colors), pch = 15, col = colors, pt.cex = 1.5, bty = "n", cex = 0.7)


######  All posible binary clade combinations  #####

rm(list = ls())

# Generate all possible combinations of binary states for the clades
num_clades <- 8
combinations <- expand.grid(replicate(num_clades, 0:1, simplify = FALSE))
combinations <- combinations[-c(1, 256),]
colnames(combinations) <- c("A", "B", "C", "D", "E", "F", "G", "I")

# Initialize matrices to store results
combis <- matrix(nrow = 0, ncol = ncol(combinations))
lambdas <- matrix(nrow = 0, ncol = 2)
pvals <- matrix(nrow = 0, ncol = 1)

# Re-read the data
tree <- read.tree("Dys_tr_names.tre")
data <- read.table("cl3_completed.txt", header = TRUE, sep = "\t")
rownames(data) <- data$Comp1

# Filter the data to include only the specified clades

data <- data[data$group %in% clades_of_interest,]

# Check for discrepancies between the tree and the data
obj <- name.check(tree, data)
tree <- drop.tip(tree, obj$tree_not_data)
obj <- name.check(tree, data)
print(obj)

# Reorder the data to match the tree structure
y <- ReorderData(tree, data, taxa.names = "row names")

# Loop through each combination of binary states for the clades
for (i in 1:254) {
  y2 <- y
  y2$group <- revalue(y2$group, as.list(combinations[i, ]))
  
  new_tree <- tree
  data_values <- as.numeric(y2[, 4])
  names(data_values) <- row.names(y2)
  new_tree$tip.state <- data_values
  
  equation <- make.bisse(tree = new_tree, states = new_tree$tip.state)
  
  starting_params <- starting.point.bisse(new_tree)
  
  mle_fit <- find.mle(equation, starting_params)
  
  bissenull_model <- constrain(equation, lambda1 ~ lambda0, mu1 ~ mu0)
  bissenull_mle <- find.mle(bissenull_model, starting_params[-c(2, 4)])
  
  bisseAnova <- anova(mle_fit, null = bissenull_mle)
  print(bisseAnova)
  
  # Store the results if the specialist rate is greater than the generalist rate and the p-value is significant
  if (coef(mle_fit)[2] > coef(mle_fit)[1] && bisseAnova$`Pr(>|Chi|)`[2] < 0.05) {
    combis <- rbind(combis, combinations[i,])
    lambdas <- rbind(lambdas, coef(mle_fit)[1:2])
    pvals <- rbind(pvals, bisseAnova$`Pr(>|Chi|)`[2])
  }
}

# Combine the results into a single data frame
result <- cbind(combis, lambdas, pvals)
print(result)