## Gibbs sampling for mixture of factor analyzers

library(matrixStats)

rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))

#' Turn a matrix's columns into informative names
mcmcify <- function(m, name) {
  colnames(m) <- paste0(name, "[", seq_len(ncol(m)),"]")
  return(m)
}

posterior <- function(y, c0, c1, k0, k1, pst, tau, gamma, theta, eta, tau_k, tau_c, r, alpha, beta,
                      theta_tilde, eta_tilde, tau_theta, tau_eta) {
  G <- ncol(y)
  N <- nrow(y)
  
  zero_ll <- sapply(seq_len(N), function(i) {
    sum(dnorm(y[i,], c0 + k0 + pst[i], 1 / sqrt(tau), log = TRUE))
  })
  
  one_ll <- sapply(seq_len(N), function(i) {
    sum(dnorm(y[i,], c1 + k1 + pst[i], 1 / sqrt(tau), log = TRUE))
  })
        
  ll <- sum(zero_ll[gamma == 0]) + sum(one_ll[gamma == 1])

  prior <- sum(dnorm(k0, theta[1], 1 / sqrt(tau_k), log = TRUE)) +
    sum(dnorm(k1, theta[2], 1 / sqrt(tau_k), log = TRUE)) +
    sum(dnorm(c0, eta[1], 1 / sqrt(tau_c), log = TRUE)) +
    sum(dnorm(c1, eta[2], 1 / sqrt(tau_c), log = TRUE)) +
    sum(dnorm(theta, theta_tilde, 1 / sqrt(tau_theta), log = TRUE)) +
    sum(dnorm(eta, eta_tilde, 1 / sqrt(tau_eta), log = TRUE)) +
    sum(dgamma(tau, alpha, beta, log = TRUE)) +
    sum(dnorm(pst, 0, 1 / r, log = TRUE))
  
  
  return( ll + prior )
}

mfa_gibbs_constr <- function(y, iter = 2000, thin = 1, burn = iter / 2, 
                             tau_k = 1, tau_c = 1, pc_initialise = 1,
                             collapse = TRUE) {
  # iter <- 2000
  
  N <- ncol(y)
  G <- nrow(y)
  message(paste("Sampling for", N, "cells and", G, "genes"))

  ## branching hierarchy
  eta_tilde <- theta_tilde <- 0
  tau_eta <- tau_theta <- 1

  eta <- rnorm(2, eta_tilde, 1 / sqrt(tau_eta))  
  theta <- rnorm(2, theta_tilde, 1 / sqrt(tau_theta))
  
  ## c & k parameters
  c0 <- c1 <- k0 <- k1 <- rep(0, G)
  
  ## precision parameters
  alpha <- 2
  beta <- 1
  tau <- rgamma(G, shape = alpha, rate = beta)
  
  
  ## pseudotime parameters
  r <- 1
  pst <-  prcomp(t(y))$x[,pc_initialise] # rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  pst <- pst / sd(pst) # make it a little more consistent with prior assumptions
  
  ## assignments for each cell
  gamma <- sample(0:1, N, replace = TRUE) # as.numeric( pst < mean(pst)  ) #
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  eta_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = 2), "eta")
  theta_trace <- mcmcify(matrix(NA, nrow = nsamples, ncol = 2), "theta")
  k0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k0")
  k1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "k1")
  c0_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "c0")
  c1_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "c1")
  tau_trace <- mcmcify(matrix(NA, nrow = G_dim[1], ncol = G_dim[2]), "tau")
  gamma_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "gamma")
  pst_trace <- mcmcify(matrix(NA, nrow = N_dim[1], ncol = N_dim[2]), "pst")
  lp_trace <- matrix(NA, nrow = nsamples, ncol = 1)
  colnames(lp_trace) <- "lp__"
  
  y <- t(y) # cell by gene! # whyyyy
  
  for(it in 1:iter) {

    
    ## set which and N parameters
    which_0 <- which(gamma == 0)
    which_1 <- which(gamma == 1)
    N_0 <- length(which_0)
    N_1 <- length(which_1)
    
    ## update for k0
    lamk_0 <- sapply(seq_len(G), function(g) {
      tau_k + tau[g] * sum(pst[which_0]^2)
    })
    
    nuk0 <- sapply(seq_len(G), function(g) {
      tau_k * theta[1] + tau[g] * sum(pst[which_0] * (y[which_0, g] - c0[g]))
    })
    nuk0 <- nuk0 / lamk_0
    
    k0_new <- rnorm(G, nuk0, 1 / sqrt(lamk_0))
    
    ## update for k1
    lamk_1 <- sapply(seq_len(G), function(g) {
      tau_k + tau[g] * sum(pst[which_1]^2)
    })
    
    nuk1 <- sapply(seq_len(G), function(g) {
      tau_k * theta[2] + tau[g] * sum(pst[which_1] * (y[which_1, g] - c1[g]))
    })
    nuk1 <- nuk1 / lamk_1
    k1_new <- rnorm(G, nuk1, 1 / sqrt(lamk_1))
  
    ## update for c0
    lamc_0 <- sapply(seq_len(G), function(g) {
      tau_c + N_0 * tau[g]
    })
    
    nuc_0 <- sapply(seq_len(G), function(g) {
      tau_c * eta[1] + tau[g] * sum(y[which_0,g] - k0_new[g] * pst[which_0])
    })
    nuc_0 <- nuc_0 / lamc_0
    c0_new <- rnorm(G, nuc_0, 1 / sqrt(lamc_0))
    
    ## update for c1
    lamc_1 <- sapply(seq_len(G), function(g) {
      tau_c + N_1 * tau[g]
    })
    
    nuc_1 <- sapply(seq_len(G), function(g) {
      tau_c * eta[2] + tau[g] * sum(y[which_1,g] - k1_new[g] * pst[which_1])
    })
    nuc_1 <- nuc_1 / lamc_1
    c1_new <- rnorm(G, nuc_1, 1 / sqrt(lamc_1))
    
    ## updates for pst
    pst_updates <- sapply(seq_len(N), function(i) {
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
    })
    pst_new <- rnorm(N, pst_updates[1,], 1 / sqrt(pst_updates[2,]))
    
    ## updates for tau
    ## create a mu vector first for convenience
    mu <- sapply(seq_len(G), function(g) {
      (1 - gamma) * (c0_new[g] + k0_new[g] * pst_new) +
        gamma * (c1_new[g] + k1_new[g] * pst_new)
    })
    
    alpha_new <- rep(alpha + N / 2, N)
    beta_new <- beta + colSums((y - mu)^2) / 2
    tau_new <- rgamma(G, alpha_new, beta_new)
    
    ## updates for theta (k)
    lambda_theta <- tau_theta + G * tau_k
    nu_theta <- tau_theta * theta_tilde + tau_k * c(sum(k0_new), sum(k1_new))
    nu_theta <- nu_theta / lambda_theta
    
    theta_new <- rnorm(2, nu_theta, 1 / sqrt(lambda_theta))
    
    ## updates for eta (c)
    lambda_eta <- tau_eta + G * tau_c
    nu_eta <- tau_eta * eta_tilde + tau_c * c(sum(c0_new), sum(c1_new))
    nu_eta <- nu_eta / lambda_eta
    
    eta_new <- rnorm(2, nu_eta, 1 / sqrt(lambda_eta))
    
    ## update for gamma
    
    if(!collapse) {
      pi <- sapply(seq_len(N), function(i) {
        y_i <- y[i,]
        comp0 <- sum(dnorm(y_i, mean = c0_new + k0_new * pst_new[i], 1 / sqrt(tau_new), log = TRUE))
        comp1 <- sum(dnorm(y_i, mean = c1_new + k1_new * pst_new[i], 1 / sqrt(tau_new), log = TRUE))
        pi_i <- comp0 - logSumExp(c(comp0, comp1))
        return(exp(pi_i))
      })
    } else {
      ## collapsed update for gamma
      pi <- sapply(seq_len(N), function(i) {
        y_i <- y[i,]
  
        # comp0_mean <- (eta_new[1] + theta_new[1] * pst_new[i]) * rep(1, G)
        # comp0_var <- 1 / tau + 1 / tau_c + k0^2 * 1 / tau_k
        # 
        # comp1_mean <- (eta_new[2] + theta_new[2] * pst_new[i]) * rep(1, G)
        # comp1_var <- 1 / tau + 1 / tau_c + k1^2 * 1 / tau_k
        # 
        # comp0 <- sum(dnorm(y_i, mean = comp0_mean, sqrt(comp0_var), log = TRUE))
        # comp1 <- sum(dnorm(y_i, mean = comp1_mean, sqrt(comp1_var), log = TRUE))
        
        comp0_mean <- eta_new[1] + k0 * pst_new[i]
        comp0_var <- 1 / tau + 1 / tau_c
        
        comp1_mean <- eta_new[2] + k1 * pst_new[i]
        comp1_var <- 1 / tau + 1 / tau_c
        
        comp0 <- sum(dnorm(y_i, mean = comp0_mean, sqrt(comp0_var), log = TRUE))
        comp1 <- sum(dnorm(y_i, mean = comp1_mean, sqrt(comp1_var), log = TRUE))
        
        pi_i <- comp0 - logSumExp(c(comp0, comp1))
        return(exp(pi_i))
      })
    }
    gamma <- rbernoulli(1 - pi)
    
    
    ## accept new parameters
    # gamma <- gamma_new
    k0 <- k0_new
    k1 <- k1_new
    c0 <- c0_new
    c1 <- c1_new
    pst <- pst_new
    tau <- tau_new
    eta <- eta_new
    theta <- theta_new
    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      k0_trace[sample_pos,] <- k0
      k1_trace[sample_pos,] <- k1
      c0_trace[sample_pos,] <- c0
      c1_trace[sample_pos,] <- c1
      tau_trace[sample_pos,] <- tau
      gamma_trace[sample_pos,] <- gamma
      pst_trace[sample_pos,] <- pst
      theta_trace[sample_pos,] <- theta
      eta_trace[sample_pos,] <- eta
      
      post <- posterior(y, c0, c1, k0, k1, pst,
                tau, gamma, theta, eta, tau_k, tau_c, r,
                alpha, beta, theta_tilde, 
                eta_tilde, tau_theta, tau_eta)
      lp_trace[sample_pos,] <- post 
    }
  }
  return(list(k0_trace = k0_trace, k1_trace = k1_trace,
              c0_trace = c0_trace, c1_trace = c1_trace,
              tau_trace = tau_trace, gamma_trace = gamma_trace,
              pst_trace = pst_trace, theta_trace = theta_trace,
              eta_trace = eta_trace, lp_trace = lp_trace))
}


to_ggmcmc <- function(g) {
  x <- do.call(cbind, g)
  mcmcl <- mcmc.list(list(mcmc(x)))
  return(ggs(mcmcl))
}


