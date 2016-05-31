## Gibbs sampling for mixture of factor analyzers

library(matrixStats)

rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))

mfa_gibbs <- function(y, iter = 2000, thin = 1, burn = iter / 2) {
  # iter <- 2000
  
  N <- ncol(y)
  G <- nrow(y)
  message(paste("Sampling for", N, "cells and", G, "genes"))

  ## assignments for each cell
  gamma <- sample(0:1, N, replace = TRUE)
  
  ## c & k parameters
  c0 <- c1 <- k0 <- k1 <- rep(0, G)
  
  ## precision parameters
  alpha <- 2
  beta <- 1
  tau <- rgamma(G, shape = alpha, rate = beta)
  
  
  ## percision hyperpriors
  tau_k <- tau_c <- 1
  
  ## pseudotime parameters
  r <- 1
  pst <-prcomp(t(y))$x[,1] # rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  k0_trace <- matrix(NA, nrow = G_dim[1], ncol = G_dim[2])
  k1_trace <- matrix(NA, nrow = G_dim[1], ncol = G_dim[2])
  c0_trace <- matrix(NA, nrow = G_dim[1], ncol = G_dim[2])
  c1_trace <- matrix(NA, nrow = G_dim[1], ncol = G_dim[2])
  tau_trace <- matrix(NA, nrow = G_dim[1], ncol = G_dim[2])
  gamma_trace <- matrix(NA, nrow = N_dim[1], ncol = N_dim[2])
  pst_trace <- matrix(NA, nrow = N_dim[1], ncol = N_dim[2])
  
  y <- t(y) # cell by gene!
  
  for(it in 1:iter) {
    ## update for gamma
    pi <- sapply(seq_len(N), function(i) {
      y_i <- y[i,]
      comp0 <- sum(dnorm(y_i, mean = c0 + k0 * pst[i], 1 / sqrt(tau), log = TRUE))
      comp1 <- sum(dnorm(y_i, mean = c1 + k1 * pst[i], 1 / sqrt(tau), log = TRUE))
      pi_i <- comp0 - logSumExp(c(comp0, comp1))
      return(exp(pi_i))
    })
    gamma_new <- rbernoulli(1 - pi)
    
    ## set which and N parameters
    which_0 <- which(gamma_new == 0)
    which_1 <- which(gamma_new == 1)
    N_0 <- length(which_0)
    N_1 <- length(which_1)
    
    ## update for k0
    lamk_0 <- sapply(seq_len(G), function(g) {
      tau_k + tau[g] * sum(pst[which_0]^2)
    })
    
    nuk0 <- sapply(seq_len(G), function(g) {
      tau[g] * sum(pst[which_0] * (y[which_0, g] - c0[g]))
    })
    nuk0 <- nuk0 / lamk_0
    
    k0_new <- rnorm(G, nuk0, 1 / sqrt(lamk_0))
    
    ## update for k1
    lamk_1 <- sapply(seq_len(G), function(g) {
      tau_k + tau[g] * sum(pst[which_1]^2)
    })
    
    nuk1 <- sapply(seq_len(G), function(g) {
      tau[g] * sum(pst[which_1] * (y[which_1, g] - c1[g]))
    })
    nuk1 <- nuk1 / lamk_1
    k1_new <- rnorm(G, nuk1, 1 / sqrt(lamk_1))
  
    ## update for c0
    lamc_0 <- sapply(seq_len(G), function(g) {
      tau_c + N_0 * tau[g]
    })
    
    nuc_0 <- sapply(seq_len(G), function(g) {
      tau[g] * sum(y[which_0,g] - k0_new[g] * pst[which_0])
    })
    nuc_0 <- nuc_0 / lamc_0
    c0_new <- rnorm(G, nuc_0, 1 / sqrt(lamc_0))
    
    ## update for c1
    lamc_1 <- sapply(seq_len(G), function(g) {
      tau_c + N_1 * tau[g]
    })
    
    nuc_1 <- sapply(seq_len(G), function(g) {
      tau[g] * sum(y[which_1,g] - k1_new[g] * pst[which_1])
    })
    nuc_1 <- nuc_1 / lamc_1
    c1_new <- rnorm(G, nuc_1, 1 / sqrt(lamc_1))
    
    ## updates for pst
    pst_updates <- sapply(seq_len(N), function(i) {
      k <- c <- NULL # awhhh yeah look at that R code
      if(gamma_new[i] == 0) {
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
    })
    pst_new <- rnorm(N, pst_updates[1,], 1 / sqrt(pst_updates[2,]))
    
    ## updates for tau
    ## create a mu vector first for convenience
    mu <- sapply(seq_len(G), function(g) {
      (1 - gamma_new) * (c0_new[g] + k0_new[g] * pst_new) +
        gamma_new * (c1_new[g] + k1_new[g] * pst_new)
    })
    
    alpha_new <- rep(alpha + N / 2, N)
    beta_new <- beta + colSums((y - mu)^2) / 2
    tau_new <- rgamma(G, alpha_new, beta_new)
    
    ## accept new parameters
    gamma <- gamma_new
    k0 <- k0_new
    k1 <- k1_new
    c0 <- c0_new
    c1 <- c1_new
    pst <- pst_new
    tau <- tau_new
    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      k0_trace[sample_pos,] <- k0
      k1_trace[sample_pos,] <- k1
      c0_trace[sample_pos,] <- c0
      c1_trace[sample_pos,] <- c1
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
    }
  }
  return(list(k0_trace = k0_trace, k1_trace = k1_trace,
              c0_trace = c0_trace, c1_trace = c1_trace,
              tau_trace = tau_trace, gamma_trace = gamma_trace,
              pst_trace = pst_trace))
}


