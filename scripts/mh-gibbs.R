## Metropolis-Hastings-within-Gibbs for mixture of nonlinear factor analyzers

library(matrixStats)
library(truncnorm)
library(ggmcmc)
library(testthat)

rbernoulli <- function(pi) rbinom(length(pi), 1, pi)

sigmoid <- function(t, k, delta) {
  return( 2  / (1 + exp(-k * (t - delta))) )
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
eval_likelihood <- function(y, pst, k0, k1, phi0, phi1, delta0, delta1, 
                            tau, gamma, tau_k, tau_phi, tau_delta,
                            alpha, beta) {
  G <- ncol(y)
  N <- nrow(y)
  
  # mean assuming all cells are on branch 0
  mu0 <- t(sapply(pst, function(t) mu_cg(k0, phi0, delta0, t)))
  # mean assuming all cells are on branch 1
  mu1 <- t(sapply(pst, function(t) mu_cg(k1, phi1, delta1, t)))
  # likelihood assuming all cells are on branch 0
  ll0 <- sapply(seq_len(G), function(g) dnorm(y[,g], mu0[,g], 1 / sqrt(tau[g]), log = TRUE))
  # likelihood assuming all cells are on branch 1
  ll1 <- sapply(seq_len(G), function(g) dnorm(y[,g], mu1[,g], 1 / sqrt(tau[g]), log = TRUE))
  
  ll <- sum(ll0[gamma == 0, ]) + sum(ll1[gamma == 1, ])
  
  ## Now find priors for everything
  lp <- sdnl(k0, 0, tau_k) + sdnl(k1, 0, tau_k)
  lp <- lp + sdnl(phi0, phi1, tau_phi) + sdnl(phi1, phi0, tau_phi)
  lp <- lp + sdnl(delta0, 0.5, tau_delta) + sdnl(delta1, 0.5, tau_delta)
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
  tau_k <- tau_phi <- 1
  tau_delta <- 1

  mu1 <- phi0 * sigmoid(t[1], k0, delta)
  mu2 <- phi1 * sigmoid(t[2], k1, delta)
  mu3 <- phi0 * sigmoid(t[3], k0, delta)

  ll <- sum(dnorm(Y[1, ], mu1, 1 / sqrt(tau), log = TRUE)) +
    sum(dnorm(Y[2, ], mu2, 1 / sqrt(tau), log = TRUE)) +
    sum(dnorm(Y[3, ], mu3, 1 / sqrt(tau), log = TRUE))
  
  lp <- sum(dnorm(c(k0, k1), 0, 1 / sqrt(tau_k), log = TRUE)) +
    2 * sum(dnorm(phi0, phi1, 1 / sqrt(tau_phi), log = TRUE)) +
    2 * sum(dnorm(delta, 0.5, 1 / sqrt(tau_delta), log = TRUE)) +
    sum(dgamma(tau, shape = alpha, rate = beta, log = TRUE)) +
    sum(dnorm(t, 0.5, 1, log = TRUE)) 

  expect_equal( (ll + lp),  
  eval_likelihood(Y, t, k0, k1, phi0, phi1, delta, delta, tau, gamma, 
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
  return(kp)
}

#' The "Hastings" correction required by using truncated normal distributions
truncation_correction <- function(tp, t, kappa_t) {
  sum(log(dtruncnorm(t, 0, 1, tp, kappa_t))) - 
    sum(log(dtruncnorm(tp, 0, 1, t, kappa_t)))
}

#' Transition ratio; _p denotes proposed variable
Q <- function(y, pst, pst_p,
              k0, k0_p, k1, k1_p,
              phi0, phi1,
              delta0, delta0_p,
              delta1, delta1_p,
              tau, gamma, tau_k, tau_phi, tau_delta,
              alpha, beta, kappa_t,
              include_truncation_correction = TRUE) {
  q <- eval_likelihood(y, pst_p, k0_p, k1_p, phi0, phi1, 
                       delta0_p, delta1_p,
                       tau, gamma, tau_k, tau_phi, tau_delta, alpha, beta) - 
    eval_likelihood(y, pst, k0, k1, phi0, phi1, 
                    delta0, delta1, 
                    tau, gamma, tau_k, tau_phi, tau_delta,
                    alpha, beta) 
  
  if(include_truncation_correction) {
    q <- q + truncation_correction(pst_p, pst, kappa_t)
  }
  
  return(q)
}

sample_pst <- function(how = c("individual", "joint"),
                       y, pst, k0, k1, phi0, phi1, delta0, delta1, 
                       tau, gamma, tau_k, tau_phi, tau_delta,
                       alpha, beta, kappa_t) {
  how <- match.arg(how)
  if(how == "individual") {
    t <- pst
    acceptance <- rep(0, length(t))
    for(i in sample(seq_along(t))) {
      tp <- propose_t(t, kappa_t, which = i)
      q <- Q(y, t, tp,
             k0, k0, k1, k1,
             phi0, phi1,
             delta0, delta0,
             delta1, delta1,
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
      if(q > log(runif(1))) {
        t <- tp
        acceptance[i] <- 1
      }
    }
    return(list(pst = t, acceptance = mean(acceptance)))
  } else {
    t <- pst
    tp <- propose_t(t, kappa_t)
    q <- Q(y, t, tp,
           k0, k0, k1, k1,
           phi0, phi1, 
           delta0, delta0,
           delta1, delta1,
           tau, gamma, tau_k, tau_phi, tau_delta,
           alpha, beta, kappa_t)
    accept <- q > log(runif(1))
    if(accept) {
      return(list(pst = tp, acceptance = 1))
    } else {
      return(list(pst = t, acceptance = 0))
    }
  }
}

sample_k <- function(par = c("k0", "k1"), how = c("individual", "joint"), 
                     y, pst, k0, k1, phi0, phi1, delta0, delta1,
                     tau, gamma, tau_k, tau_phi, tau_delta,
                     alpha, beta, kappa_k, kappa_t) {
  how <- match.arg(how)
  par <- match.arg(par)
  if(how == "individual") {
    k <- switch(par, k0 = k0, k1 = k1)
    acceptance <- rep(0, length(k))
    for(i in seq_along(k)) {
      kp <- propose_kd(k, kappa_k, which = i)
      q <- NULL
      if(par == "k0") {
      q <- Q(y, pst, pst, k0, kp, k1, k1,
             phi0, phi1, delta0, delta0,
             delta1, delta1,
             tau, gamma, tau_k, tau_phi, tau_delta, alpha, beta, kappa_t)
      } else {
      q <- Q(y, pst, pst, k0, k0, k1, kp,
             phi0, phi1, delta0, delta0,
             delta1, delta1,
             tau, gamma, tau_k, tau_phi, tau_delta, alpha, beta, kappa_t)   
      }
      #print(c(i,q, kp))
      if(q > log(runif(1))) {
        k <- kp
        acceptance[i] <- 1
      }
    }
    return(list(k = k, acceptance = mean(acceptance)))
  } else {
    q <- NULL
    if(par == "k0") {
      k <- k0
      kp <- propose_kd(k0, kappa_t)
      q <- Q(y, pst, pst, k0, kp, k1, k1,
             phi0, phi1, delta0, delta0,
             delta1, delta1, 
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
    } else {
      k <- k1
      kp <- propose_kd(k1, kappa_t)
      q <- Q(y, pst, pst, k0, k0, k1, kp,
             phi0, phi1, delta0, delta0, 
             delta1, delta1,
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
    }
    accept <- q > log(runif(1))
    if(accept) {
      return(list(k = kp, acceptance = 1))
    } else {
      return(list(k = k, acceptance = 0))
    }
  }
}

sample_delta <- function(par = c("delta0", "delta1"),
                         how = c("individual", "joint"),
                      y, pst, k0, k1, phi0, phi1, delta0, delta1, 
                       tau, gamma, tau_k, tau_phi, tau_delta,
                       alpha, beta, kappa_delta, kappa_t) {
  par <- match.arg(par)
  how <- match.arg(how)
  if(how == "individual") {
    d <- switch(par, delta0 = delta0, delta1 = delta1)
    acceptance <- rep(0, length(d))
    for(i in seq_along(d)) {
      dp <- propose_kd(d, kappa_delta, which = i)
      q <- NULL
      if(par == "delta0") {
        q <- Q(y, pst, pst,
               k0, k0, k1, k1,
               phi0, phi1, 
               delta0, dp,
               delta1, delta1,
               tau, gamma, tau_k, tau_phi, tau_delta,
               alpha, beta, kappa_t)
      } else {
        q <- Q(y, pst, pst,
               k0, k0, k1, k1,
               phi0, phi1, 
               delta0, delta0,
               delta1, dp,
               tau, gamma, tau_k, tau_phi, tau_delta,
               alpha, beta, kappa_t)
      }
      if(q > log(runif(1))) {
        d <- dp
        acceptance[i] <- 1
      }
    }
    return(list(delta = d, acceptance = mean(acceptance)))
  } else {
    q <- d <- NULL
    if(par == "delta0") {
      d <- delta0
      dp <- propose_kd(d, kappa_delta)
      q <- Q(y, pst, pst,
             k0, k0, k1, k1,
             phi0, phi1, 
             d, dp,
             delta1, delta1,
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
    } else {
      d <- delta1
      dp <- propose_kd(d, kappa_delta)
      q <- Q(y, pst, pst,
             k0, k0, k1, k1,
             phi0, phi1, 
             delta0, delta0,
             d, dp,
             tau, gamma, tau_k, tau_phi, tau_delta,
             alpha, beta, kappa_t)
    }

    accept <- q > log(runif(1))
    if(accept) {
      return(list(delta = dp, acceptance = 1))
    } else {
      return(list(delta = d, acceptance = 0))
    }
  }
}

mcmc_status_message <- function(it, iter) {
  ## Want to status every 10% of iterations
  it10 <- round(0.1 * iter)
  if(it %% it10 == 0) {
    message(paste0("Iteration ", it, " of ", iter, " (", round(it / iter * 100), "%)"))
  }
}



#' MH-within-Gibbs for nonlinear MFA
#' @param y Cell-by-gene
#' @param fixed A named list, which can contain 'pst', 'k0', 'k1' or 'delta', where
#' if the name appears then those values are held fixed
mh_gibbs <- function(y, iter = 2000, thin = 1, burn = iter / 2,
                     proposals = list(kappa_t = 0.05, kappa_k = 0.5, kappa_delta = 0.1),
                     fixed = list(), ptype = c("joint", "individual"),
                     pst_init = NULL) {
  ptype <- match.arg(ptype)

  G <- ncol(y)
  N <- nrow(y)

  message(paste("Sampling for", N, "cells and", G, "genes"))

  ## assignments for each cell
  gamma <- sample(0:1, N, replace = TRUE)
  
  if(is.null(pst_init)) {
    pst <- rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  } else {
    pst <- pst_init
  }
  
  if(is.null(pst_init)) {
    k0 <- k1 <- rep(0, G)
  } else {
    k0 <- k1 <- apply(y, 2, function(yy) coef(lm(yy ~ pst_init))[2])
  }
  phi0 <- phi1 <- colMeans(y)
  delta0 <- delta1 <- rep(0.5, G)
  
  ## proposal stdevs
  kappa_t <- proposals$kappa_t
  kappa_k <- proposals$kappa_k
  kappa_delta <- proposals$kappa_delta
  
  ## precision parameters
  alpha <- 2
  beta <- 1
  tau <- 1 / matrixStats::colVars(y)
  
  ## percision hyperpriors
  tau_k <- tau_phi <- 1
  tau_delta <- 10
  
  ## pseudotime parameters
  r <- 1
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  k0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k0")
  k1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k1")
  phi0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "phi0")
  phi1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "phi1")
  delta0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "delta0")
  delta1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "delta1")
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  lp_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = 1), "lp__")
  
  ## MH control - to tidy up
  accept_reject <- list(pst = rep(NA, iter), 
                        delta0 = rep(NA, iter), delta1 = rep(NA, iter), 
                        k0 = rep(NA, iter), k1 = rep(NA, iter))
  
  for(it in 1:iter) {
    mcmc_status_message(it, iter)
        
    ## Sanity checks - remove later
    stopifnot(!(any(pst < 0) | any(pst > 1)))
    
    ## update for gamma
    if("gamma" %in% names(fixed)) {
      gamma_new <- fixed$gamma
    } else {
      pi <- sapply(seq_len(N), function(i) {
        y_i <- y[i,]
        comp0 <- sum(dnorm(y_i, mean = mu_cg(k0, phi0, delta0, pst[i]), 1 / sqrt(tau), log = TRUE))
        comp1 <- sum(dnorm(y_i, mean = mu_cg(k1, phi1, delta1, pst[i]), 1 / sqrt(tau), log = TRUE))
        pi_i <- comp0 - logSumExp(c(comp0, comp1))
        return(exp(pi_i))
      })
      gamma_new <- rbernoulli(1 - pi)
    }
    
    ## set which and N parameters
    which_0 <- which(gamma_new == 0)
    which_1 <- which(gamma_new == 1)
    N_0 <- length(which_0)
    N_1 <- length(which_1)
    

    ## MH updates for t, k & delta ---------------------------------------------
    
    ## Pseudotimes
    if("pst" %in% names(fixed)) {
      pst_new <- fixed$pst
    } else {
      pst_new_list <- sample_pst(ptype, y, pst, k0, k1, phi0, phi1, 
                                 delta0, delta1, tau, gamma, tau_k, tau_phi, 
                                 tau_delta, alpha, beta, kappa_t)
      pst_new <- pst_new_list$pst
      accept_reject$pst[it] <- pst_new_list$acceptance
    }
    
    ## k0
    if("k0" %in% names(fixed)) {
       k0_new <- fixed$k0
    } else {
      k0_new_list <- sample_k("k0", ptype, y, pst, k0, k1, phi0, phi1, 
                              delta0, delta1, tau, gamma, tau_k, tau_phi, 
                              tau_delta, alpha, beta, kappa_k, kappa_t)
      k0_new <- k0_new_list$k
      accept_reject$k0[it] <- k0_new_list$acceptance
    }
    
    ## k1
    if("k1" %in% names(fixed)) {
      k1_new <- fixed$k1
    } else {
      k1_new_list <- sample_k("k1", ptype, y, pst, k0, k1, phi0, phi1, 
                              delta0, delta1, tau, gamma, tau_k, tau_phi, 
                              tau_delta, alpha, beta, kappa_k, kappa_t)
      k1_new <- k1_new_list$k
      accept_reject$k1[it] <- k1_new_list$acceptance
    }
    
    
    
    ## delta0
    if("delta0" %in% names(fixed)) {
      delta0_new <- fixed$delta0
    } else {
      delta0_new_list <- sample_delta("delta0", "joint", y, pst, k0, k1, phi0, phi1, 
                                      delta0, delta1, tau, gamma, tau_k, tau_phi, 
                                      tau_delta, alpha, beta, kappa_delta, kappa_t)
      delta0_new <- delta0_new_list$delta
      accept_reject$delta0[it] <- delta0_new_list$acceptance
    }
    
    ## delta1
    if("delta1" %in% names(fixed)) {
      delta1_new <- fixed$delta1
    } else {
      delta1_new_list <- sample_delta("delta1", "joint", y, pst, k0, k1, phi0, phi1, 
                                      delta0, delta1, 
                                      tau, gamma, tau_k, tau_phi, tau_delta,
                                      alpha, beta, kappa_delta, kappa_t)
      delta1_new <- delta1_new_list$delta
      accept_reject$delta1[it] <- delta1_new_list$acceptance
    }
    
    #' Gibbs updates ---------------------------------
    
    #' Updates for phi ---------
    #' In this section, mu0 is defined slightly differently from usual,
    #' dropping the factor of phi beforehand
    mu0 <- t(sapply(pst_new, function(t) sigmoid(t, k0_new, delta0_new)))
    mu1 <- t(sapply(pst_new, function(t) sigmoid(t, k1_new, delta1_new)))
    
    mu0 <- mu0[gamma_new == 0, , drop = FALSE]
    mu1 <- mu1[gamma_new == 1, , drop = FALSE]
    
    if(length(tau) != length(colSums(mu0^2))) {
      message("Length mismatch!")
      print(c(length(tau), length(colSums(mu0^2))))
    }
    
    lam0 <- tau_phi + tau * colSums(mu0^2)
    nu0 <- (tau_phi * phi1 + tau * colSums(y[gamma_new == 0, ,drop=FALSE] * mu0)) / lam0
    phi0_new <- rnorm(G, nu0, 1 / sqrt(lam0))
    
    lam1 <- tau_phi + tau * colSums(mu1^2)
    nu1 <- (tau_phi * phi0_new + tau * colSums(y[gamma_new == 1, ,drop=FALSE] * mu1)) / lam1
    phi1_new <- rnorm(G, nu1, 1 / sqrt(lam1))
    
    
    #' Updates for tau ------
    #' Here we revert back to the usual definition of mu
    mu0 <- t(sapply(pst_new, function(t) mu_cg(k0_new, phi0_new, delta0_new, t)))
    mu1 <- t(sapply(pst_new, function(t) mu_cg(k1_new, phi1_new, delta0_new, t)))
    
    mu <- matrix(NA, nrow = N, ncol = G)
    mu[gamma_new == 0, ] <- mu0[gamma_new == 0, ]
    mu[gamma_new == 1, ] <- mu1[gamma_new == 1, ]
    
    alpha_new <- rep(alpha + N / 2, N)
    beta_new <- beta + colSums((y - mu)^2) / 2
    tau_new <- rgamma(G, alpha_new, beta_new)
    
    ## accept new parameters (this mainly here for debugging)
    gamma <- gamma_new
    k0 <- k0_new
    k1 <- k1_new
    phi0 <- phi0_new
    phi1 <- phi1_new
    delta0 <- delta0_new
    delta1 <- delta1_new
    pst <- pst_new
    tau <- tau_new
    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      k0_trace[sample_pos,] <- k0
      k1_trace[sample_pos,] <- k1
      phi0_trace[sample_pos,] <- phi0
      phi1_trace[sample_pos,] <- phi1
      delta0_trace[sample_pos,] <- delta0
      delta1_trace[sample_pos,] <- delta1
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
      lp_trace[sample_pos,] <- eval_likelihood(y, pst, k0, k1, phi0, phi1, 
                                               delta0, delta1, tau, gamma, 
                                               tau_k, tau_phi, tau_delta, alpha, beta)
    }
  }
  
  accept_reject <- lapply(accept_reject, `[`, burn:iter)
  
  return(list(traces = list(k0_trace = k0_trace, k1_trace = k1_trace,
                            phi0_trace = phi0_trace, phi1_trace = phi1_trace, 
                            delta0_trace = delta0_trace, delta1_trace = delta1_trace,
                            tau_trace = tau_trace, gamma_trace = gamma_trace,
                            pst_trace = pst_trace, lp_trace = lp_trace),
              accept = accept_reject))
}

to_ggmcmc <- function(g) {
  x <- do.call(cbind, g$traces)
  mcmcl <- mcmc.list(list(mcmc(x)))
  return(ggs(mcmcl))
}

