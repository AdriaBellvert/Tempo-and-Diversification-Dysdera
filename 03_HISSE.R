
library(ape)       
library(geiger)    
library(plyr)       
library(hisse)      
library(phytools) 


tree <- read.tree("Dys_tr_names.tre")
data <- read.table("cl3_completed.txt", header = TRUE, sep = "\t")
rownames(data) <- data$Comp1
colnames(data) <- c("spp", "Comp1", "Comp2", "group")

# Check for discrepancies between the tree and the data
obj <- name.check(tree, data)
tree <- drop.tip(tree, obj$tree_not_data)
new.tree <- tree
y <- data

# Revalue the groups to binary states for BiSSE analysis
y$group <- revalue(y$group, c("D" = 0, "E" = 0, "A" = 0, "G" = 1, "B" = 1, "F" = 1,
                              "I" = 1, "C" = 1, "H" = NA))
y <- na.omit(y)
y <- ReorderData(tree, y, taxa.names = "row names")

data.v <- as.numeric(y[, 4])
names(data.v) <- row.names(y)
new.tree$tip.state <- data.v
y2 <- y[, c(1, 4)]

# Create transition rate matrices for HiSSE analysis
rates.hisse <- TransMatMakerHiSSE(hidden.traits = 1)
rates.bisse <- TransMatMakerHiSSE(hidden.traits = 0)

# Fit the BiSSE model using HiSSE framework
bisse.hmle <- hisse(new.tree, y2, turnover = c(1, 2), eps = c(1, 2), hidden.states = FALSE, trans.rate = rates.bisse)

# Fit the CID model (constant rates) using HiSSE framework
cid.mle <- hisse(new.tree, y2, turnover = c(1, 1), eps = c(1, 1), hidden.states = FALSE, trans.rate = rates.bisse)

# Modify the rates matrix for CID-2 model
rates.cid2 <- rates.hisse
rates.cid2[!is.na(rates.cid2)] <- 1
print(rates.cid2)

# Fit the CID-2 model with hidden states
cid2.mle <- hisse(new.tree, y2, f = c(1, 1), turnover = c(1, 1, 2, 2), eps = c(1, 1, 2, 2), hidden.states = TRUE, trans.rate = rates.cid2)

# Perform marginal reconstruction for CID-2 model
cid2.recon <- MarginReconHiSSE(new.tree, y2, f = cid2.mle$f, pars = cid2.mle$solution, hidden.states = 2)

# Create transition rate matrix for CID-4 model
rates.cid4 <- TransMatMakerHiSSE(hidden.traits = 3)
rates.cid4[!is.na(rates.cid4)] <- 1

# Fit the CID-4 model with hidden states
cid4.mle <- hisse(new.tree, y2, f = c(1, 1), turnover = c(1, 1, 2, 2, 3, 3, 4, 4), eps = c(1, 1, 2, 2, 3, 3, 4, 4), hidden.states = TRUE, trans.rate = rates.cid4)

# Fit the HiSSE model with multiple hidden states
hisse.mle <- hisse(new.tree, y2, f = c(1, 1), hidden.states = TRUE, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), trans.rate = rates.cid4)

# Function to extract log-likelihood with degrees of freedom attribute
logLik.hisse.fit <- function(x, ...) {
  lik <- x$loglik
  attr(lik, "df") <- (x$AIC + 2 * lik) / 2
  lik
}

# Compare models using AIC and Akaike weights
model_comparison <- data.frame(
  model = c("CID", "BISSE", "HISSE CID-2", "HISSE CID-4", "HiSSE"),
  logL = sapply(list(cid.mle, bisse.hmle, cid2.mle, cid4.mle, hisse.mle), logLik),
  k = sapply(list(cid.mle, bisse.hmle, cid2.mle, cid4.mle, hisse.mle), function(x) attr(logLik(x), "df")),
  AIC = aic <- sapply(list(cid.mle, bisse.hmle, cid2.mle, cid4.mle, hisse.mle), AIC),
  Akaike.weight = unclass(aic.w(aic))
)
print(model_comparison)