library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(viridis)
library(ggmcmc)

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

dp <- data.frame(prcomp(X)$x[,1:2], branch = as.factor(branch), pseudotime = true_t)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() + scale_color_viridis()
ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

source("scripts/gibbs.R")

g <- mfa_gibbs(t(X))

tmap <- posterior.mode(mcmc(g$pst_trace))
gamma_mean <- colMeans(g$gamma_trace)

dp2 <- data.frame(prcomp(X)$x[,1:2], branch = gamma_mean, pseudotime = tmap)

ggplot(dp2, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() +
  scale_color_viridis(name = "MAP\npseudotime")
ggplot(dp2, aes(x = PC1, y = PC2, colour = branch)) + geom_point() +
  scale_color_viridis(name = "MAP\ngamma")


data_frame(true_t, tmap) %>% 
  ggplot(aes(x = true_t, y = tmap)) +
  geom_point(shape = 21, fill = "grey", color = "black", size = 2) +
  ylab("MAP pseudotime") + xlab("True pseudotime") + geom_rug(alpha = 0.5)


# MH - gibbs --------------------------------------------------------------

source("scripts/mh-gibbs.R")

k_true <- h5read(fname, "basic_branching/k")
phi_true <- h5read(fname, "basic_branching/phi")
delta_true <- h5read(fname, "basic_branching/delta")
tau_true <- 1 / h5read(fname, "basic_branching/stdev")^2
# fixed <- list(k0 = k_true[,1], k1 = k_true[,2], delta = delta_true[,1])

# fixed <- list(pst = true_t, delta0 = delta_true[,1], delta1 = delta_true[,2])
fixed <- list(delta0 = delta_true[,1], delta1 = delta_true[,2])
proposals = list(kappa_t = 0.5, kappa_k = 10, kappa_delta = 0.1)
g <- mh_gibbs(X, iter = 20000, thin = 10, proposals = proposals, 
              fixed = fixed, ptype = "joint")

s <- to_ggmcmc(g)
ggs_traceplot(s, "lp__") + stat_smooth()
ggs_running(s, "lp__")
ggs_autocorrelation(s, "lp__")

sapply(g$accept, mean)

tmap <- posterior.mode(mcmc(g$traces$pst_trace))
tmean <- matrixStats::colMedians(g$traces$pst_trace)
plot(true_t, tmap)
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
ggs_histogram(filter(s, Parameter == "k0[1]"))
ggs_traceplot(filter(s, Parameter == "delta0[1]"))
ggs_traceplot(filter(s, Parameter == "delta1[2]"))

k0_mean <- colMeans(g$traces$k0_trace)
k1_mean <- colMeans(g$traces$k1_trace)
k0_map <- posterior.mode(mcmc(g$traces$k0_trace))
k1_map <- posterior.mode(mcmc(g$traces$k1_trace))
qplot(k_true[,1], k0_mean)
qplot(k_true[,2], k1_mean)

qplot(k_true[,1], k0_map)
qplot(k_true[,2], k1_map)

# Play around with Q ------------------------------------------------------

# Q <- function(y, pst, pst_p,
#               k0, k0_p, k1, k1_p,
#               phi0, phi1,
#               delta, delta_p, 
#               tau, gamma, tau_k, tau_phi, tau_delta,
#               alpha, beta, kappa_t) 

Q(X, tmap, true_t, #rep(0.5, length(true_t)),
  k_true[,1], k_true[,1], k_true[,2], k_true[,2],
  phi_true[,1], phi_true[,2],
  delta_true, delta_true, tau_true,
  branch_true, 1, 1, 1, 2, 1, proposals$kappa_t,
  include_truncation_correction = TRUE)

vary_q_t1 <- function(t1, w) {
  t <- true_t
  t[w] <- t1
  Q(X, true_t, t,
    k_true[,1], k_true[,1], k_true[,2], k_true[,2],
    phi_true[,1], phi_true[,2],
    delta_true, delta_true, tau_true,
    branch_true, 1, 1, 1, 2, 1, proposals$kappa_t)
}

w <- 7
tvals <- seq(0, 1, length.out = 100)
qv <- sapply(tvals, vary_q_t1, w)
ldf <- data_frame(tvals, qv)

ggplot(ldf, aes(x = tvals, y = qv)) + geom_line() + geom_vline(xintercept = true_t[w])
