// [[Rcpp::depends(EMC2)]]
#include <Rcpp.h>
#include "EMC2/userfun.hpp"

// Example: two params (q0, alpha) and three inputs (rt, B, weight)
Rcpp::NumericVector custom_kernel(Rcpp::NumericMatrix kernel_pars, 
                                  Rcpp::NumericMatrix input) {
  // assume kernels_pars(q0, alpha) and inputs(rt, B, weight)
  int n = input.nrow();
  Rcpp::NumericVector q(n);
  Rcpp::NumericVector B_t(n);
  Rcpp::NumericVector pe(n);
  // B_t = B + w*Q_t for trial 1
  q[0] = kernel_pars(0,0);
  B_t[0] = input(0, 1) + input(0, 2)*q[0];
  for (int i = 1; i < n; ++i) {
    // 'prediction error' = B_t/rt_t - Q_t 
    pe[i-1] = (B_t[i-1] / input(i-1,0)) - q[i-1];
    // update: Q_t = Q_t-1 + alpha_t-1 * pe_t-1
    q[i] = q[i-1] + kernel_pars(i-1,1) * pe[i-1];
    // Apply base to B to get the correct threshold for 
    // next trial (needed for PE calculation)
    B_t[i] = input(i,1) + input(i,2) * q[i];
  }
  return q;
}

// Export pointer maker for registration
// [[Rcpp::export]]
SEXP EMC2_make_custom_kernel_ptr();
EMC2_MAKE_PTR(custom_kernel);
