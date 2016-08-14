library(cowplot)
library(ggplot2)
library(rhdf5)
library(coda)
library(MCMCglmm)
library(viridis)
library(ggmcmc)
library(dplyr)
library(readr)
library(scater)
library(matrixStats)
library(tidyr)
library(magrittr)
library(Rcpp)


# Synthetic ---------------------------------------------------------------

fname <- "data/synthetic.h5"
X <- h5read(fname, "basic_branching/X")
branch <- h5read(fname, "basic_branching/branch_assignment")
true_t <- h5read(fname, "basic_branching/pseudotime")

dp <- data.frame(prcomp(X)$x[,1:2], branch = as.factor(branch), pseudotime = true_t)

ggplot(dp, aes(x = PC1, y = PC2, colour = pseudotime)) + geom_point() + scale_color_viridis()
ggplot(dp, aes(x = PC1, y = PC2, colour = branch)) + geom_point()

source("houdini/R/houdini.R")
sourceCpp("houdini/src/gibbs.cpp")

g <- mfa_cpp(t(X), iter = 20000, thin = 10, collapse = TRUE)


mc <- to_ggmcmc(g)
ggs_traceplot(filter(mc, Parameter == "lp__")) + stat_smooth()


library(profvis)
profvis({
  g <- mfa_cpp(t(X), iter = 200, thin = 1, collapse = TRUE)
})

system.time( mfa_gibbs_semi_ard(t(X), iter = 10000, thin = 1, collapse = TRUE) )
#' 10000 iterations:
#'    user  system elapsed 
#'   61.987   2.623  66.210 
#' 

lamk_0 <- tau_k + tau * sum(pst[which_0]^2)
lamk_0_2 <- calculate_lamk(tau_k, tau, pst, which_0_l)


sample_k_native <- function(y, pst, c0, tau, theta, tau_k, which_0) {
  lamk_0 <- tau_k + tau * sum(pst[which_0]^2)
  
  nuk0 <- sapply(seq_len(G), function(g) {
    tau_k[g] * theta[g] + tau[g] * sum(pst[which_0] * (y[which_0, g] - c0[g]))
  })
  #if(length(which_0) > 0) nuk0 <- nuk0 
  nuk0 <- nuk0 / lamk_0
  
  # print(tau_k[17])
  
  k0_new <- rnorm(G, nuk0, 1 / sqrt(lamk_0))
  return(k0_new)
}

# sample_k(y, pst, c0, tau, theta, tau_k, which_0_l)

system.time(replicate(1000, sample_k_native(y, pst, c0, tau, theta, tau_k, which_0_l)))
system.time(replicate(1000, sample_k(y, pst, c0, tau, theta, tau_k, which_0_l)))



# C -----------------------------------------------------------------------

nuc_0 <- sapply(seq_len(G), function(g) {
  tau_c * eta[1] + tau[g] * sum(y[which_0,g] - k0_new[g] * pst[which_0])
})

nuc_0_2 <- calculate_nuc(y, pst, k0_new, tau, eta[1], tau_c, which_0_l)



sample_c_native <- function() {
  lamc_0 <- tau_c + N_0 * tau
  
  nuc_0 <- sapply(seq_len(G), function(g) {
    tau_c * eta[1] + tau[g] * sum(y[which_0,g] - k0_new[g] * pst[which_0])
  })
  
  nuc_0 <- nuc_0 / lamc_0
  c0_new <- rnorm(G, nuc_0, 1 / sqrt(lamc_0))
  return(c0_new)
}

system.time(replicate(1000, sample_c_native()))
system.time(replicate(1000, sample_c(y, pst, k0_new, tau, eta[1], tau_c, which_0_l, N_0)))


# Pseudotimes -------------------------------------------------------------

pst_updates <- t(sapply(seq_len(N), function(i) {
  k <- c <- NULL # awhhh yeah look at that R code
  if(gamma[i] == 0) {
    k <- k0_new
    c <- c0_new
  } else {
    k <- k1_new
    c <- c1_new
  }
  lam_ti <- r^2 + sum(tau * k^2)
  
  nu_ti <- sum(tau * k * (y[i,] - c))
  nu_ti <- nu_ti / lam_ti
  return(c(nu_ti, lam_ti))
}))

pst_updates_2 <- pst_update_par(y, k0_new, k1_new, c0_new, c1_new, 1, gamma, tau)



# Testing pi --------------------------------------------------------------

pi <- sapply(seq_len(N), function(i) {
  y_i <- y[i,]
  comp0 <- sum(dnorm(y_i, mean = c0_new + k0_new * pst_new[i], 1 / sqrt(tau_new), log = TRUE))
  comp1 <- sum(dnorm(y_i, mean = c1_new + k1_new * pst_new[i], 1 / sqrt(tau_new), log = TRUE))
  pi_i <- comp0 - log_sum_exp(c(comp0, comp1))
  return(exp(pi_i))
})

pi2 <- calculate_pi(y, c0_new, c1_new, k0_new, k1_new, gamma, pst_new, tau_new, FALSE)

rpi <- runif(100)
system.time(replicate(100, rbernoulli(rpi)))
system.time(replicate(100, r_bernoulli_vec(rpi)))

# Overall -----------------------------------------------------------------
system.time( mfa_gibbs_semi_ard(t(X), iter = 5000, thin = 1, collapse = FALSE) )
system.time( mfa_cpp(t(X), iter = 10000, thin = 1, collapse = FALSE) )
