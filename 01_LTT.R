
library(ape)
library(laser)
library(phytools)

my_tree_dysderids <- read.tree("Dys_tr_names.tre")

# Calculate the lineage-through-time (LTT) plot data for the tree
dysderid_ltt <- ltt(my_tree_dysderids)

# Calculate the Monte Carlo Constant Rate (MCCR) test statistics for the LTT data
dysderid_mccr <- mccr(dysderid_ltt)
print(dysderid_mccr)
plot(dysderid_mccr)


# Read multiple phylogenetic trees from a file
dysdera_tree_multi_tip <- read.tree("Dys_multi_tr_names.trees")

# Initialize a matrix to store gamma statistics for each tree
mccr_dy <- matrix(nrow = length(dysdera_tree_multi_tip), ncol = 1)

# Loop through each tree and calculate the gamma statistic from the LTT data
for (i in seq_along(dysdera_tree_multi_tip)) {
  dysderid_ltt <- ltt(dysdera_tree_multi_tip[[i]])
  mccr_dy[i, 1] <- dysderid_ltt[["gamma"]]
}
plot(density(mccr_dy), main = "Density of Gamma Statistics")

# Generate the LTT plot for all trees with logarithmic scaling
temp <- ltt(dysdera_tree_multi_tip, log = TRUE)
plot(temp)


# Loop through each tree and fit various diversification models, storing AIC and LH values
vls <- matrix(nrow = length(dysdera_tree_multi_tip), ncol = 12)
for (i in 1:1001) {
  dysderid_btimes <- branching.times(dysdera_tree_multi_tip[[i]])
  div_models_dysderids <- fitdAICrc(dysderid_btimes, modelset = c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints = 100)
  vls[i, 1:6] <- div_models_dysderids[1:6, 12]
  vls[i, 7:12] <- div_models_dysderids[1:6, 5]
}

colnames(vls) <- c("AICpureBirth", "AICbd", "AICDDX", "AICDDL", "AICyule2rate", "AICyule3rate",
                   "LHpureBirth", "LHbd", "LHDDX", "LHDDL", "LHyule2rate", "LHyule3rate")


vls_df <- as.data.frame(vls)
stacked_df <- stack(vls_df)
boxplot(stacked_df$values ~ stacked_df$ind, col = rainbow(ncol(vls_df)),
        main = "Boxplot of AIC and LH Values for Different Models",
        xlab = "Model", ylab = "Value")
