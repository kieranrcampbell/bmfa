library(ggplot2)
library(rhdf5)

set.seed(123L)

#' Sigmoid function for activations
sigmoid <- function(t, phi, k, delta) {
  return( 2 * phi / (1 + exp(-k*(t - delta))))
}

C <- 100 # cells
G <- 40 # genes

branch <- rbinom(C, 1, 0.5)

gsd <- sqrt(1 / rgamma(G, 1, 1))


## We assume first G / 2 (= 20) genes are common to both branches, and the 
## final G / 2 genes exhibit branching structure. We want to build in the 
## fact that delta < 0.5 for the common genes and delta > 0.5 for the 
## branch specific genes

k <- replicate(2, runif(G, 5, 10) * sample(c(-1, 1), G, replace = TRUE))
phi <- replicate(2, runif(G, 5, 10))
delta <- replicate(2, runif(G, 0.5, 1))

inds <- 1:(G / 2)
inds2 <- (G/2 + 1):G
k[, 2] <- k[inds, 1]

k[inds2, ] <- t(apply(k[inds2, ], 1, function(r) r * sample(c(0, 1))))

phi[, 1] <- phi[, 2]
delta[inds, 2] <- delta[inds, 1] <- runif(G / 2, 0, 0.5)

## Now make it look like a branching process
for(r in inds2) {
  whichzero <- which(k[r,] == 0)
  nonzero <- which(k[r,] != 0)
  k_sign <- sign(k[r,nonzero])
  if(k_sign == 1) {
    phi[r, whichzero] <- 0
  } else {
    phi[r, whichzero] <- 2 * phi[r, nonzero]
  }
}

pst <- runif(C)

X <- sapply(seq_along(branch), function(i) {
  k_i <- k[, branch[i] + 1]
  phi_i <- phi[, branch[i] + 1]
  delta_i <- delta[, branch[i] + 1]
  mu <- sigmoid(pst[i], phi_i, k_i, delta_i)
  rnorm(length(mu), mu, gsd)
})

fname <- file.path("data", "synthetic.h5")
if(!file.exists(fname)) {
  h5createFile(fname)
}

h5createGroup(fname, "basic_branching") # we'll call this dataset basic_branching
h5write(t(X), fname, "basic_branching/X")
h5write(pst, fname, "basic_branching/pseudotime")
h5write(branch, fname, "basic_branching/branch_assignment")
h5write(k, fname, "basic_branching/k")
h5write(phi, fname, "basic_branching/phi")
h5write(delta, fname, "basic_branching/delta")
h5write(gsd, fname, "basic_branching/stdev")

d <- data.frame(t(X), pseudotime = pst, branch = as.factor(branch))

dm <- reshape2::melt(d, id.vars = c("pseudotime", "branch"),
                     variable.name = "gene", value.name = "expression")

ggplot(dm, aes(x = pseudotime, y = expression, colour = branch)) +
  geom_point() + facet_wrap(~ gene, scales = "free_y")

ggsave("figs/basic_branching.png", width = 10, height = 8)