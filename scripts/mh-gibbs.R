## Gibbs sampling for mixture of factor analyzers

library(matrixStats)
library(truncnorm)

rbernoulli <- function(pi) rbinom(length(pi), 1, pi)

sigmoid <- function(t, k, delta) {
  2  / (1 + exp(-k * (t - delta)))
}

#' Turn a matrix's columns into informative names
mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}

mu_cg <- function(k, phi, delta, t) {
  phi * sigmoid(t, k, delta)
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

test_likelihood <- function() {
  Y <- matrix(1:6, ncol = 2)
  t <- c(0.3, 0.6, 0.9)
  gamma <- c(0, 1, 0)
  tau <- c(0.5, 0.5)
  phi0 <- c(1, 1)
  phi1 <- c(2, 2)
  k0 <- c(1, 1)
  k1 <- c(-1, -1)
  delta <- c(0.5, 0.5)
  alpha <- 2; beta <- 1
  tau_k <- tau_phi <- tau_delta <- 1

  mu1 <- phi0 * sigmoid(t[1], k0, delta)
  mu2 <- phi1 * sigmoid(t[2], k1, delta)
  mu3 <- phi0 * sigmoid(t[3], k0, delta)

  ll <- sum(dnorm(Y[1, ], mu1, 1 / sqrt(tau), log = TRUE)) +
    sum(dnorm(Y[2, ], mu2, 1 / sqrt(tau), log = TRUE)) +
    sum(dnorm(Y[3, ], mu3, 1 / sqrt(tau), log = TRUE))
  
  lp <- sum(dnorm(c(k0, k1), 0, 1 / sqrt(tau_k), log = TRUE)) +
    2 * sum(dnorm(phi0, phi1, 1 / sqrt(tau_phi), log = TRUE)) +
    sum(dnorm(delta, 0.5, 1 / sqrt(tau_delta), log = TRUE)) +
    sum(dgamma(tau, shape = alpha, rate = beta, log = TRUE)) +
    sum(dnorm(t, 0.5, 1, log = TRUE))
  
  stopifnot( (ll + lp) ==   
  eval_likelihood(Y, t, k0, k1, phi0, phi1, delta, tau, gamma, 
                  tau_k, tau_phi, tau_delta, alpha, beta))
}

#' Propose a pseudotime based on current pseudotimes (t). If
#' which is a particular index, then the vector returned contains
#' only that index changed
propose_t <- function(t, kappa_t, which = seq_along(t)) {
  tp <- t
  tp[which] <- rtruncnorm(length(which), a = 0, b = 1, mean = t[which], sd = kappa_t)
  return(tp)
}

#' Proposal for either kappa or delta
propose_kd <- function(k, kappa_k, which = seq_along(k)) {
  kp <- k
  kp[which] <- rnorm(length(which), mean = k[which], sd = kappa_k)
}

#' The "Hastings" correction required by using truncated normal distributions
truncation_correction <- function(tp, t, kappa_t) {
  sum(log(dtruncnorm(t, 0, 1, tp, kappa_t))) - 
    sum(log(dtruncnorm(tp, 0, 1, t, kappa_t)))
}

#' Transition ratio; _p denotes proposed variable
Q <- function(y, pst, pst_p,
              k0, k0_p, k1, k1_p,
              phi0, phi0_p, phi1, phi1_p, 
              delta, delta_p, 
              tau, gamma, tau_k, tau_phi, tau_delta,
              alpha, beta, kappa_t) {
  q <- eval_likelihood(y, pst_p, k0_p, k1_p, phi0_p, phi1_p, delta_p, 
  tau, gamma, tau_k, tau_phi, tau_delta,
  alpha, beta) - 
    eval_likelihood(y, pst, k0, k1, phi0, phi1, delta, 
                    tau, gamma, tau_k, tau_phi, tau_delta,
                    alpha, beta) +
    truncation_correction(pst_p, pst, kappa_t)
  return(q)
}

sample_pst <- function(how = c("individual", "joint"),
                       y, pst, k0, k1, phi0, phi1, delta, 
                       tau, gamma, tau_k, tau_phi, tau_delta,
                       alpha, beta, kappa_t) {
  how <- match.arg(how)
  if(how == "individual") {
    t <- pst
    acceptance <- rep(0, seq_along(t))
    for(i in seq_along(t)) {
      tp <- propose_t(t, kappa_t, which = i)
      q <- Q(y, t, tp,
             k0, k0, k1, k1,
             phi0, phi0, phi1, phi1, 
             delta, delta, 
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
      if(q > log(runif(1))) {
        t <- tp
        acceptance[i] <- 1
      }
    }
    return(list(pst = t, acceptance = mean(acceptance)))
  } else {
    tp <- propose_t(pst, kappa_t)
    q <- Q(y, t, tp,
           k0, k0, k1, k1,
           phi0, phi0, phi1, phi1, 
           delta, delta, 
           tau, gamma, tau_k, tau_phi, tau_delta,
           alpha, beta, kappa_t)
    accept <- q > log(runif(1))
    if(accept) {
      return(list(pst = tp, acceptance = 1))
    } else {
      return(list(pst = t, acceptance = 1))
    }
  }
}

sample_k <- function(how = c("individual", "joint"),
                       y, pst, k0, k1, phi0, phi1, delta, 
                       tau, gamma, tau_k, tau_phi, tau_delta,
                       alpha, beta, kappa_k, kappa_t) {
  how <- match.arg(how)
  if(how == "individual") {
    k <- k0
    acceptance <- rep(0, seq_along(k))
    for(i in seq_along(k)) {
      kp <- propose_t(k, kappa_k, which = i)
      q <- Q(y, t, t,
             k0, kp, k1, k1,
             phi0, phi0, phi1, phi1, 
             delta, delta, 
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
      if(q > log(runif(1))) {
        t <- tp
        acceptance[i] <- 1
      }
    }
    return(list(pst = t, acceptance = mean(acceptance)))
  } else {
    tp <- propose_t(pst, kappa_t)
    q <- Q(y, t, tp,
           k0, k0, k1, k1,
           phi0, phi0, phi1, phi1, 
           delta, delta, 
           tau, gamma, tau_k, tau_phi, tau_delta,
           alpha, beta, kappa_t)
    accept <- q > log(runif(1))
    if(accept) {
      return(list(pst = tp, acceptance = 1))
    } else {
      return(list(pst = t, acceptance = 1))
    }
  }
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
    pst_new_list <- sample_pst("joint", Y, pst, k0, k1, phi0, phi1, delta, 
                          tau, gamma, tau_k, tau_phi, tau_delta,
                          alpha, beta, kappa_t)
    pst_new <- pst_new_list$pst
    accept_reject$t[it] <- pst_new_list$acceptance
    
    
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
    phi0 <- phi0_new
    phi1 <- phi1_new
    delta <- delta_new
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


