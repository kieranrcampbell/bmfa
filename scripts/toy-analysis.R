library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
pst <- h5read(fname, "basic_branching/pseudotime")

dp <- data.frame(prcomp(X)$x[,1:2], branch = as.factor(branch), pseudotime = pst)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point()
ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

source("scripts/gibbs.R")

g <- mfa_gibbs(t(X))

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(X)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point()
ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

