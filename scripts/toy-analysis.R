library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(viridis)
library(ggmcmc)


# MH - gibbs --------------------------------------------------------------

source("scripts/mh-gibbs.R")

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

k_true <- h5read(fname, "basic_branching/k")
phi_true <- h5read(fname, "basic_branching/phi")
delta_true <- h5read(fname, "basic_branching/delta")
tau_true <- 1 / h5read(fname, "basic_branching/stdev")^2
# fixed <- list(k0 = k_true[,1], k1 = k_true[,2], delta = delta_true[,1])

pc1 <- prcomp(X)$x[,1]
pst_init <- (pc1 - min(pc1) + 0.01) / (max(pc1) - min(pc1) + 0.02)

# fixed <- list(pst = true_t, delta0 = delta_true[,1], delta1 = delta_true[,2])
fixed <- list(delta0 = delta_true[,1], delta1 = delta_true[,2])#, gamma = branch)#,

proposals = list(kappa_t = 0.1, kappa_k = 0.8, kappa_delta = 0.1); set.seed(123)
g <- mh_gibbs(X, iter = 2000, thin = 1, proposals = proposals, 
              fixed = fixed, ptype = "individual", pst_init = pst_init)

s <- to_ggmcmc(g)
ggs_traceplot(s, "lp__") + stat_smooth()
ggs_running(s, "lp__")
ggs_autocorrelation(s, "lp__")

sapply(g$accept, mean)

tmean <- colMeans(g$traces$pst_trace)
plot(true_t, tmean)
cor(true_t, tmean)

branch_true <- as.factor(h5read(fname, "basic_branching/branch_assignment"))
gamma_mean <- colMeans(g$traces$gamma_trace)


d <- data.frame(prcomp(X)$x[,1:2], pst = tmean, gamma = gamma_mean)
ggplot(d, aes(x = PC1, y = PC2, color = pst)) + geom_point() + scale_color_viridis()
ggplot(d, aes(x = PC1, y = PC2, color = gamma)) + geom_point() + scale_color_viridis()

qplot(branch_true, gamma_mean, geom = 'boxplot')

ggs_traceplot(filter(s, Parameter == "pst[1]"))
ggs_traceplot(filter(s, Parameter == "pst[2]"))

ggs_traceplot(filter(s, Parameter == "k0[1]"))
ggs_traceplot(filter(s, Parameter == "delta0[1]"))
ggs_traceplot(filter(s, Parameter == "delta1[2]"))
ggs_traceplot(filter(s, Parameter == "phi1[1]"))



k0_mean <- colMeans(g$traces$k0_trace)
k1_mean <- colMeans(g$traces$k1_trace)
qplot(k_true[,1], k0_mean)
qplot(k_true[,2], k1_mean)

phi0_mean <- colMeans(g$traces$phi0_trace)
phi1_mean <- colMeans(g$traces$phi1_trace)

qplot(phi_true[,1], phi0_mean)
qplot(phi_true[,2], phi1_mean)

taumap <- colMeans(g$traces$tau_trace)
qplot(tau_true, taumap) + scale_x_log10() + scale_y_log10()


# Play with pi ------------------------------------------------------------

k_true <- h5read(fname, "basic_branching/k")
phi_true <- h5read(fname, "basic_branching/phi")
delta_true <- h5read(fname, "basic_branching/delta")
tau_true <- 1 / h5read(fname, "basic_branching/stdev")^2

k0 <- k_true[,1]; k1 <- k_true[,2]
phi0 <- phi_true[,1]; phi1 <- phi_true[,2]
delta0 <- delta_true[,1]; delta1 <- delta_true[,2]
tau <- tau_true
pst <- true_t

pi <- sapply(seq_len(N), function(i) {
  y_i <- y[i,]
  comp0 <- sum(dnorm(y_i, mean = mu_cg(k0, phi0, delta0, pst[i]), 1 / sqrt(tau), log = TRUE))
  comp1 <- sum(dnorm(y_i, mean = mu_cg(k1, phi1, delta1, pst[i]), 1 / sqrt(tau), log = TRUE))
  pi_i <- comp0 - logSumExp(c(comp0, comp1))
  return(exp(pi_i))
})

d <- data.frame(prcomp(X)$x[,1:2], pi)
ggplot(d, aes(x = PC1, y = PC2, color = pi)) + geom_point() + scale_color_viridis()

