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

proposals = list(kappa_t = 0.3, kappa_k = 1, kappa_delta = 0.1)
k_true <- h5read(fname, "basic_branching/k")
delta_true <- h5read(fname, "basic_branching/delta")
fixed <- list(k0 = k_true[,1], k1 = k_true[,2], delta = delta_true[,1])

g <- mh_gibbs(X, iter = 1000, thin = 1, proposals = proposals)#, fixed = fixed)

s <- to_ggmcmc(g)
ggs_traceplot(s, "lp__")
ggs_running(s, "lp__")

tmap <- posterior.mode(mcmc(g$traces$pst_trace))
plot(true_t, tmap)

sapply(g$accept, mean)

branch_true <- as.factor(h5read(fname, "basic_branching/branch_assignment"))
gamma_mean <- colMeans(g$traces$gamma_trace)

d <- data.frame(prcomp(y)$x[,1:2], pst = tmap, gamma = gamma_mean)
ggplot(d, aes(x = PC1, y = PC2, color = pst)) + geom_point()
ggplot(d, aes(x = PC1, y = PC2, color = gamma)) + geom_point()

qplot(branch_true, gamma_mean, geom = 'boxplot')

ggs_traceplot(filter(s, Parameter == "pst[1]"))
ggs_traceplot(filter(s, Parameter == "pst[2]"))

ggs_traceplot(filter(s, Parameter == "k0[2]"))
ggs_traceplot(filter(s, Parameter == "delta[1]"))

