#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

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

