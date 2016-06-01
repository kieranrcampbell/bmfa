## Gibbs sampling for mixture of factor analyzers

library(matrixStats)

rbernoulli <- function(pi) rbinom(length(pi), 1, pi)
  
  #sapply(pi, function(p) sample(c(0,1), 1, 
   #                                                      prob = c(1-p,p)))

sigmoid <- function(k, delta, t) {
  2  / (1 + exp(-k * (t - delta)))
}

#' Turn a matrix's columns into informative names
mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}

mu_cg <- function(k, phi, delta, t) {
  phi * sigmoid(k, delta, t)
}

#' Sum dnorm logged with precision rather than sd
sdnl <- function(x, mean = 0, tau = 1, log = TRUE) sum(dnorm(x, mean, 1 / sqrt(tau), log))

#' 
eval_likelihood <- function(y, pst, k0, k1, phi0, phi1, delta, 
                            tau, gamma, tau_k, tau_phi, tau_delta,
                            alpha, beta) {
  # mean assuming all cells are on branch 0
  mu0 <- t(sapply(pst, function(t) mu_cg(k0, phi0, delta, t)))
  # mean assuming all cells are on branch 1
  mu1 <- t(sapply(pst, function(t) mu_cg(k1, phi1, delta, t)))
  # likelihood assuming all cells are on branch 0
  ll0 <- dnorm(y, mu0, 1 / sqrt(tau), log = TRUE)
  # likelihood assuming all cells are on branch 1
  ll1 <- dnorm(y, mu1, 1 / sqrt(tau), log = TRUE)
  
  ll <- sum(ll0[gamma == 0, ]) + sum(ll1[gamma == 1, ])
  
  ## Now find priors for everything
  lp <- sdnl(k0, 0, tau_k) + sdnl(k1, 0, tau_k)
  lp <- lp + sdnl(phi0, phi1, tau_phi) + sdnl(phi1, phi0, tau_phi)
  lp <- lp + sdnl(delta, 0.5, tau_delta)
  lp <- lp + sum(dgamma(tau, shape = alpha, rate = beta, log = TRUE))
  lp <- lp + sdnl(pst, 0.5, 1)
  
  return(lp + ll)
}


#' MH-within-Gibbs for nonlinear MFA
#' @param y Cell-by-gene
mh_gibbs <- function(y, iter = 2000, thin = 1, burn = iter / 2,
                     proposals = list(kappa_t = 0.05, kappa_k = 0.5, kappa_delta = 0.1)) {

  G <- ncol(y)
  N <- nrow(y)
  message(paste("Sampling for", N, "cells and", G, "genes"))

  ## assignments for each cell
  gamma <- sample(0:1, N, replace = TRUE)
  
  ## c & k parameters
  k0 <- k1 <- rep(0, G)
  phi0 <- phi1 <- rep(1, G)
  delta <- rep(0.5, G)
  
  ## precision parameters
  alpha <- 2
  beta <- 1
  tau <- rgamma(G, shape = alpha, rate = beta)
  
  ## percision hyperpriors
  tau_k <- tau_phi <- tau_delta <- 1
  
  ## pseudotime parameters
  r <- 1
  pst <- rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  k0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k0")
  k1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k1")
  phi0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "phi0")
  phi1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "phi1")
  delta_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "delta")
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  
  ## MH control - to tidy up
  accept_reject <- list(joint = rep(NA, iter), t = rep(NA, iter),
                        delta = rep(NA, iter), k = rep(NA, iter))
  
  for(it in 1:iter) {
    ## Sanity checks - remove later
    stopifnot(!(any(pst < 0) | any(pst > 1)))
    
    ## update for gamma
    pi <- sapply(seq_len(N), function(i) {
      y_i <- y[i,]
      comp0 <- sum(dnorm(y_i, mean = mu_cg(k0, phi0, delta, pst[i]), 1 / sqrt(tau), log = TRUE))
      comp1 <- sum(dnorm(y_i, mean = mu_cg(k1, phi1, delta, pst[i]), 1 / sqrt(tau), log = TRUE))
      pi_i <- comp0 - logSumExp(c(comp0, comp1))
      return(exp(pi_i))
    })
    gamma_new <- rbernoulli(1 - pi)
    
    ## set which and N parameters
    which_0 <- which(gamma_new == 0)
    which_1 <- which(gamma_new == 1)
    N_0 <- length(which_0)
    N_1 <- length(which_1)
    

    ## MH updates for t, k & delta ---------------------------------------------

    
    
    ## updates for phi
    
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


