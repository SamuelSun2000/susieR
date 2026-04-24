#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include "mr_ash.h"

using namespace cpp11;
using namespace arma;

// Exported: random permutation index vector of length p * numiter,
// returned as an R integer vector (0-based Armadillo-style indices).
[[cpp11::register]]
integers random_order(int p, int numiter) {
  uvec o = random_order_impl(p, numiter);
  return as_integers(o);
}

// Exported: mr.ash coordinate-ascent in suff-stat form.
//
// pi, beta, r are passed by value (as cpp11 doubles); we materialise local
// Armadillo copies, mutate locally, and return them in the result list so the
// R caller reads the updated values from there.
[[cpp11::register]]
writable::list caisa_cpp(const doubles_matrix<>& X, const doubles& y,
                         const doubles& w, const doubles& sa2,
                         const doubles& pi_init, const doubles& beta_init,
                         const doubles& r_init, double sigma2,
                         const integers& o_r,
                         int maxiter, int miniter,
                         double convtol, double epstol, std::string method_q,
                         bool updatepi, bool updatesigma,
                         bool verbose) {

  // Convert inputs (cpp11 -> Armadillo). Read-only inputs keep view semantics;
  // pi, beta, r are mutated locally, so force owned deep copies via conv_to.
  mat X_mat   = as_Mat(X);
  vec y_vec   = as_Col(y);
  vec w_vec   = as_Col(w);
  vec sa2_vec = as_Col(sa2);
  vec pi      = conv_to<vec>::from(as_Col(pi_init));
  vec beta    = conv_to<vec>::from(as_Col(beta_init));
  vec r       = conv_to<vec>::from(as_Col(r_init));
  uvec o      = as_uvec(o_r);

  // ---------------------------------------------------------------------
  // DEFINE SIZES
  // ---------------------------------------------------------------------
  int n                   = X_mat.n_rows;
  int p                   = X_mat.n_cols;
  int K                   = sa2_vec.n_elem;

  // ---------------------------------------------------------------------
  // PREDEFINE LOCAL VARIABLES
  // ---------------------------------------------------------------------
  vec varobj(maxiter);
  int iter               = 0;
  int i                  = 0;
  int j;

  double a1;
  double a2;
  vec piold;
  vec betaold;

  // ---------------------------------------------------------------------
  // PRECALCULATE
  // ---------------------------------------------------------------------
  mat S2inv              = 1 / outerAddition(1 / sa2_vec, w_vec);
  S2inv.row(0).fill(epstol);

  // ---------------------------------------------------------------------
  // START LOOP : CYCLE THROUGH COORDINATE ASCENT UPDATES
  // ---------------------------------------------------------------------
  for (iter = 0; iter < maxiter; iter++) {

    // reset parameters
    a1                   = 0;
    a2                   = 0;
    piold                = pi;
    pi.fill(0);
    betaold              = beta;

    // ---------------------------------------------------------------------
    // RUN COORDINATE ASCENT UPDATES : INDEX 1 - INDEX P
    // ---------------------------------------------------------------------
    for (j = 0; j < p; j++) {
      updatebetaj(X_mat.col(o(i)), w_vec(o(i)), beta(o(i)), r, piold, pi,
                  sigma2, sa2_vec, S2inv.col(o(i)), a1, a2, o(i), p, epstol);
      i++;
    }

    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 1
    // ---------------------------------------------------------------------
    varobj(iter)          = arma::dot(r, r) - arma::dot(square(beta), w_vec) + a1;

    // ---------------------------------------------------------------------
    // UPDATE SIGMA2 IF REQUESTED
    // ---------------------------------------------------------------------
    if (updatesigma) {
      if (method_q == std::string("sigma_indep_q")) {
        sigma2            = varobj(iter) + p * (1.0 - pi(0)) * sigma2;
        sigma2            = sigma2 / (n + p * (1.0 - pi(0)));
      } else if (method_q == std::string("sigma_dep_q")) {
        sigma2            = varobj(iter) / n;
      }
    }

    if (updatepi) {
      piold               = pi;
    }

    // ---------------------------------------------------------------------
    // CALCULATE VARIATIONAL OBJECTIVE 2
    // ---------------------------------------------------------------------
    varobj(iter)          = varobj(iter) / sigma2 / 2.0 +
                            log(2.0 * M_PI * sigma2) / 2.0 * n -
                            dot(pi, log(piold + epstol)) * p + a2;

    for (j = 1; j < K; j++) {
      varobj(iter)       += pi(j) * log(sa2_vec(j)) * p / 2;
    }

    if (!updatepi) {
      pi                  = piold;
    }

    // ---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // ---------------------------------------------------------------------
    if (iter >= miniter - 1) {
      double beta_norm = arma::norm(beta, 2);
      if (arma::norm(betaold - beta, 2) < convtol * std::max(1.0, beta_norm)) {
        iter++;
        break;
      }

      if (iter > 0) {
        if (varobj(iter) > varobj(iter - 1)) {
          break;
        }
      }
    }
  }

  if (verbose) {
    Rprintf("Mr.ASH terminated at iteration %d: max|beta|=%.4e, sigma2=%.4e, pi0=%.4f\n",
            iter, arma::max(arma::abs(beta)), sigma2, pi(0));
  }

  // ---------------------------------------------------------------------
  // RETURN VALUES
  // ---------------------------------------------------------------------
  using namespace cpp11::literals;
  writable::list out({
    "beta"_nm    = as_doubles(beta),
    "sigma2"_nm  = cpp11::as_sexp(sigma2),
    "pi"_nm      = as_doubles(pi),
    "iter"_nm    = cpp11::as_sexp(iter),
    "varobj"_nm  = as_doubles(vec(varobj.subvec(0, iter - 1)))
  });

  return out;
}
