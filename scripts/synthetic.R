library(ggplot2)


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

k[inds2, ] <- apply(k[inds2, ], 1, function(r) r * sample(c(0, 1)))

phi[, 2] <- phi[inds, 2]
delta[, 2] <- delta[inds, 1] <- runif(G / 2, 0, 0.5)

pst <- runif(C)

X <- sapply(seq_along(branch), function(i) {
  k_i <- k[, branch[i] + 1]
  phi_i <- phi[, branch[i] + 1]
  delta_i <- delta[, branch[i] + 1]
  mu <- sigmoid(pst[i], phi_i, k_i, delta_i)
  rnorm(length(mu), mu, gsd)
})

d <- data.frame(t(X), pseudotime = pst, branch = as.factor(branch))

dm <- reshape2::melt(d, id.vars = c("pseudotime", "branch"),
                     variable.name = "gene", value.name = "expression")

ggplot(dm, aes(x = pseudotime, y = expression, colour = branch)) +
  geom_point() + facet_wrap(~ gene, scales = "free_y")

