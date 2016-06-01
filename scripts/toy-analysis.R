library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

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




# MH - gibbs --------------------------------------------------------------

source("scripts/mh-gibbs.R")

proposals = list(kappa_t = 0.1, kappa_k = 1, kappa_delta = 0.1)
g <- mh_gibbs(X, iter = 2000, proposals = proposals, burn = 1000)

s <- to_ggmcmc(g)
ggs_traceplot(s, "lp__")
ggs_running(s, "lp__")

tmap <- posterior.mode(mcmc(g$traces$pst_trace))
plot(true_t, tmap)

sapply(g$accept, mean)

ggs_traceplot(filter(s, Parameter == "pst[1]"))
ggs_traceplot(filter(s, Parameter == "pst[2]"))

ggs_traceplot(filter(s, Parameter == "k0[2]"))
ggs_traceplot(filter(s, Parameter == "delta[1]"))

