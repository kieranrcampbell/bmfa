#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;


/***
 * Convenience functions to compute log normal density and log sum of exponentials
 */

// [[Rcpp::export]]
double log_d_norm(double x, double mu, double sigma) {
  double ldn = -0.5 * log(2 * PI);
  ldn -= log(sigma) + 1 / (2 * sigma * sigma) * (x - mu) * (x - mu);
  return ldn;
}

// [[Rcpp::export]]
double log_sum_exp(NumericVector x) {
  double y;
  y = log(sum(exp(x - max(x)))) + max(x);
  return(y);
}

// [[Rcpp::export]]
IntegerVector r_bernoulli_vec(NumericVector pi) {
  int N = pi.size();
  IntegerVector gamma(N);
  NumericVector rands = runif(N);
  for(int i = 0; i < N; i++) {
    if(rands[i] < pi[i]) {
      gamma[i] = 1;
    } else {
      gamma[i] = 0;
    }
  }
  return gamma;
}

/*******
 * K sampling here
 ******/

NumericVector calculate_nuk(NumericMatrix y, NumericVector pst, NumericVector c,
                       NumericVector tau, NumericVector theta, NumericVector tau_k,
                       LogicalVector which_l) {
  int N = pst.size();
  int G = c.size();
  
  NumericVector nu_k(G);

  for(int g = 0; g < G; g++) {
    double tmp = 0;
    for(int i = 0; i < N; i++) {
      if(which_l[i])
        tmp += pst[i] * ( y(i,g) - c[g] );
    }
    nu_k[g] = tau_k[g] * theta[g] + tau[g] * tmp;
  }  
  
  return nu_k;
}


NumericVector calculate_lamk(NumericVector tau_k, NumericVector tau, 
                             NumericVector pst, LogicalVector which_l) {
  int G = tau.size();
  int N = pst.size();
  
  NumericVector lamk(G);

  double pst_sum = 0;  
  for(int i = 0; i < N; i++) {
    if(which_l[i])
      pst_sum += pst[i] * pst[i];
  }
  
  lamk = tau * pst_sum + tau_k;

  return lamk;
}

// [[Rcpp::export]]
NumericVector sample_k(NumericMatrix y, NumericVector pst, NumericVector c,
                       NumericVector tau, NumericVector theta, NumericVector tau_k,
                       LogicalVector which_l) {
  int G = c.size();
  NumericVector nuk = calculate_nuk(y, pst, c, tau, theta, tau_k, which_l);
  NumericVector lamk = calculate_lamk(tau_k, tau, pst, which_l);
  
  for(int g = 0; g < G; g++) 
    nuk[g] /= lamk[g];
  
  NumericVector k_new(G);
  
  for(int g = 0; g < G; g++)
    k_new[g] = as<double>(rnorm(1, nuk[g], 1 / sqrt(lamk[g])));
  
  return k_new;
}


/*******
 * C sampling here
 ******/

//[[Rcpp::export]]
NumericVector calculate_nuc(NumericMatrix y, NumericVector pst, NumericVector k,
                            NumericVector tau, double eta, double tau_c,
                            LogicalVector which_l) {
  int G = k.size();
  int N = pst.size();
  
  NumericVector nuc(G);
  
  for(int g = 0; g < G; g++) {
    double inner_sum = 0;
    for(int i = 0; i < N; i++) {
      if(which_l[i])
        inner_sum += y(i,g) - k[g] * pst[i];
    }
    nuc[g] = tau_c * eta + tau[g] * inner_sum;
  }
  
  
  return nuc;
}

// [[Rcpp::export]]
NumericVector calculate_lamc(NumericVector tau, double tau_c, int N) {
  int G = tau.size();
  NumericVector lamc(G);
  for(int g = 0 ; g < G; g++)
    lamc[g] =  tau_c  + N * tau[g];
  return lamc;
}

// [[Rcpp::export]]
NumericVector sample_c(NumericMatrix y, NumericVector pst, NumericVector k,
                       NumericVector tau, double eta, double tau_c,
                       LogicalVector which_l, int N) {
  int G = k.size();
  NumericVector nuc = calculate_nuc(y, pst, k, tau, eta, tau_c, which_l);
  NumericVector lamc = calculate_lamc(tau, tau_c, N);
  
  for(int g = 0; g < G; g++) 
    nuc[g] /= lamc[g];
  
  NumericVector c_new(G);
  
  for(int g = 0; g < G; g++)
    c_new[g] = as<double>(rnorm(1, nuc[g], 1 / sqrt(lamc[g])));
  
  return c_new;
}

/** 
 * Pseudotime updates
 * */

/* 
 *First column is lambda, second is mean. We explicitly calculate this to check 
 *consistency with previous results.
 */

NumericMatrix pst_update_par(NumericMatrix y, NumericVector k0, NumericVector k1, NumericVector c0, NumericVector c1,
                          double r, NumericVector gamma, NumericVector tau) {
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix pst_parameters(N, 2); // what we return
  
  // Placehold vectors for current branch
  NumericVector k(G);
  NumericVector c(G);
  
  double lam_ti, nu_ti;
  
  for(int i = 0; i < N; i++) {
    nu_ti = 0;
    
    if(gamma[i] == 0) {
      k = k0;
      c = c0;
    } else {
      k = k1;
      c = c1;
    }
    
    lam_ti = pow(r, 2) + sum(tau * pow(k, 2));
    for(int g = 0; g < G; g++) 
      nu_ti += tau[g] * k[g] * (y(i,g) - c[g]);
    
    nu_ti /= lam_ti;
    pst_parameters(i, 0) = nu_ti;
    pst_parameters(i, 1) = lam_ti;
  }
  
  return pst_parameters;
}

// [[Rcpp::export]]
NumericVector sample_pst(NumericMatrix y, NumericVector k0, NumericVector k1, NumericVector c0, NumericVector c1,
                         double r, NumericVector gamma, NumericVector tau) {
  int N = y.nrow();
  
  NumericMatrix pst_pars(N, 2);
  
  pst_pars = pst_update_par(y, k0, k1, c0, c1, r, gamma, tau);
  
  NumericVector pst_new(N);
  for(int i = 0; i < N; i++) {
    pst_new[i] = as<double>(rnorm(1, pst_pars(i, 0), 1 / sqrt(pst_pars(i, 1))));
  }
  
  return pst_new;
}


NumericMatrix tau_params(NumericMatrix y, NumericVector c0, NumericVector c1, NumericVector k0, NumericVector k1,
                         NumericVector gamma, NumericVector pst, double alpha, double beta) {
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_beta(G, 2); // new alpha and beta: first column alpha, second beta
  NumericMatrix mu(N, G);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      if(gamma[i] == 0) {
        mu(i,g) = c0[g] + k0[g] * pst[i];
      } else {
        mu(i,g) = c1[g] + k1[g] * pst[i];
      }
    }
  }
  
  for(int g = 0; g < G; g++) {
    alpha_beta(g, 0) = alpha + N / 2;
    double beta_new = beta;
    for(int i = 0; i < N; i++)
      beta_new += pow(y(i,g) - mu(i,g), 2) / 2;
    
    alpha_beta(g, 1) = beta_new;
  }
  
  return alpha_beta;
  
  
}

// [[Rcpp::export]]
NumericVector sample_tau(NumericMatrix y, NumericVector c0, NumericVector c1, NumericVector k0, NumericVector k1,
                         NumericVector gamma, NumericVector pst, double alpha, double beta) {
  // N = y.nrow();
  int G = y.ncol();
  
  NumericVector tau(G);
  NumericMatrix alpha_beta(G, 2);
  alpha_beta = tau_params(y, c0, c1, k0, k1, gamma, pst, alpha, beta);
  
  for(int g = 0; g < G; g++) {
    tau[g] = as<double>(rgamma(1, alpha_beta(g, 0), 1 / alpha_beta(g, 1))); // !!! RCPP gamma parametrised by shape - scale
  }
  
  return tau;
}



// [[Rcpp::export]]
NumericVector calculate_pi(NumericMatrix y, NumericVector c0, NumericVector c1, NumericVector k0, NumericVector k1,
                           NumericVector gamma, NumericVector pst, NumericVector tau, 
                           NumericVector eta, double tau_c, bool collapse) {
  int N = y.nrow();
  int G = y.ncol();
  
  NumericVector pi(N);
  if(collapse == 0) {
    for(int i = 0; i < N; i++) {
      double comp0 = 0, comp1 = 0;
      for(int g = 0; g < G; g++) {
        double y_ = y(i,g);
        double comp0_mean = c0[g] + k0[g] * pst[i];
        double comp1_mean = c1[g] + k1[g] * pst[i];
        double sd = 1 / sqrt(tau[g]);
        

        comp0 += log_d_norm(y_, comp0_mean, sd);
        comp1 += log_d_norm(y_, comp1_mean, sd);
      }
      NumericVector comb(2);
      comb[0] = comp0; comb[1] = comp1;
      pi(i) = exp(comp0 - log_sum_exp(comb));
    }
  } else {
    for(int i = 0; i < N; i++) {
      double comp0 = 0, comp1 = 0;
      for(int g = 0; g < G; g++) {
        double y_ = y(i,g);
        double comp0_mean = eta[0] + k0[g] * pst[i];
        double comp1_mean = eta[1] + k1[g] * pst[i];
        double sd = sqrt(1 / tau_c + 1 / sqrt(tau[g]));
        
        comp0 += log_d_norm(y_, comp0_mean, sd);
        comp1 += log_d_norm(y_, comp1_mean, sd);
      }
      NumericVector comb(2);
      comb[0] = comp0; comb[1] = comp1;
      pi(i) = exp(comp0 - log_sum_exp(comb));
    }
  }
  return pi;
}


/***
 * x sampling
 ***/

// [[Rcpp::export]]
NumericMatrix sample_x(NumericMatrix x, LogicalMatrix is_dropout,
                       NumericVector c0, NumericVector c1, NumericVector k0, NumericVector k1,
                           NumericVector gamma, NumericVector pst, NumericVector tau, double lambda) {
  int N = x.nrow();
  int G = x.ncol();
  
  NumericVector k(G);
  NumericVector c(G);
  
  NumericMatrix x_new(N, G);
  

  for(int i = 0; i < N; i++) {
    if(gamma[i] == 0) {
      k = k0;
      c = c0;
    } else {
      k = k1;
      c = c1;
    }
    for(int g = 0; g < G; g++) {
      if(is_dropout(i,g) == true) {
        double mu_ig = c[g] + k[g] * pst[i]; 
        x_new(i,g) = as<double>(rnorm(1, mu_ig + lambda / (N * tau[g]), 1 / sqrt(tau[g])));
      } else {
        x_new(i,g) = x(i,g);
      }
    }
  }
  
  
  return x_new;
}
