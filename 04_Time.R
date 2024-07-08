
library(phytools)  
library(geiger)    


my.tree.agamids <- read.tree("Dys_tr_names.tre")

# Create a birth-death likelihood function for the tree
bd_model <- make.bd(my.tree.agamids)

# Find maximum likelihood estimates for the birth-death model parameters
st <- c(0.03, 0.01)
bd_mle <- find.mle(func = bd_model, x.init = st)

# Create a time-dependent birth-death model with exponential speciation rate
bvar_model <- make.bd.t(my.tree.agamids, functions = c("exp.t", "constant.t"))


# Find maximum likelihood estimates for the time-dependent birth-death model parameters
st <- c(0.03, 0.01, 0.01)
bvar_mle <- find.mle(bvar_model, st)

# Compare the birth-death model and the time-dependent birth-death model 
anova(bd_mle, bvar_mle)
AIC(bd_mle, bvar_mle)

# Perform MCMC (Markov Chain Monte Carlo) to sample from the posterior distribution of the time-dependent birth-death model parameters
bvar_res_bayes <- mcmc(lik = bvar_model, x.init = st, nsteps = 1000, w = 0.01, print.every = 100)
mean(bvar_res_bayes$lambda.a)
sd(bvar_res_bayes$lambda.a)

# Remove the first 100 samples as burn-in
bvar_res_bayes <- bvar_res_bayes[-(1:100),]


postSamples <- bvar_res_bayes[, c("lambda.l", "lambda.a", "mu")]
thinnedPosterior <- postSamples[round(seq(1, nrow(postSamples), length.out = 100)),]

# Plot the time-dependent speciation and extinction rates
plot(NULL, xlim = c(0, 20), ylim = c(0, 0.2), bty = "n", xlab = "Time",
     las = 1, ylab = expression(paste("rate (", lambda, " or ", mu, ")")), cex.axis = 0.8)

legend("topright", lwd = 3, col = c("black", "lightgrey"),
       legend = c(expression(paste("Speciation (", lambda, ")")),
                  expression(paste("Extinction (", mu, ")"))),
       cex = 0.8, bty = "n")

t <- seq(0, 70, length.out = 100)

for (i in 1:nrow(thinnedPosterior)) {
  lambda.l <- thinnedPosterior[i, "lambda.l"]
  lambda.a <- thinnedPosterior[i, "lambda.a"]
  mu <- thinnedPosterior[i, "mu"]
  b <- lambda.l * exp(-lambda.a * t)
  d <- rep(mu, length(t))
  lines(t, b, col = make.transparent("black", 0.5), lwd = 1)
  lines(t, d, col = make.transparent("lightgrey", 0.5), lwd = 1)
}
