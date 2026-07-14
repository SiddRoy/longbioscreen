// ==============================================================
// objective.cpp
//
// Stan Math implementation of the Stage-2 conditional longitudinal
// likelihood in model.Rmd and likelihood.R. Only the VALUE is written;
// gradients are obtained with stan::math::gradient on a thin summed
// functor.
//
// Constant-term convention:
//   We omit the additive normalizing constant -n_i / 2 * log(2 * pi),
//   exactly as likelihood.R's copied cluster_ll() does.
// ==============================================================

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
// [[Rcpp::plugins(cpp17)]]

#ifdef _REENTRANT
#undef _REENTRANT
#endif
#define _REENTRANT 1

#include <stan/math.hpp>
#include <RcppEigen.h>
#include <vector>

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ------------------------------------------------------------
// Fixed longitudinal data. del=1 and del=0 groups are separate
// blocks, matching the data-setup rule in model.Rmd.
// Cluster starts are 0-based inclusive; ends are 0-based exclusive.
// ------------------------------------------------------------
struct LmmCM03Data {
  VectorXd y1;
  VectorXd t1;
  VectorXd Ti1;
  VectorXd fTi1;
  std::vector<int> start1;
  std::vector<int> end1;
  int J1;

  VectorXd y0;
  VectorXd t0;
  VectorXd Ti_til0;
  std::vector<int> start0;
  std::vector<int> end0;
  int J0;

  LmmCM03Data(const Map<VectorXd>& y1_,
              const Map<VectorXd>& t1_,
              const Map<VectorXd>& Ti1_,
              const Map<VectorXd>& fTi1_,
              const Rcpp::IntegerVector& cluster_start1_,
              const Rcpp::IntegerVector& cluster_end1_,
              const Map<VectorXd>& y0_,
              const Map<VectorXd>& t0_,
              const Map<VectorXd>& Ti_til0_,
              const Rcpp::IntegerVector& cluster_start0_,
              const Rcpp::IntegerVector& cluster_end0_)
    : y1(y1_), t1(t1_), Ti1(Ti1_), fTi1(fTi1_), J1(cluster_start1_.size()),
      y0(y0_), t0(t0_), Ti_til0(Ti_til0_), J0(cluster_start0_.size()) {
    start1.reserve(J1);
    end1.reserve(J1);
    for (int j = 0; j < J1; ++j) {
      start1.push_back(cluster_start1_[j] - 1);
      end1.push_back(cluster_end1_[j]);
    }

    start0.reserve(J0);
    end0.reserve(J0);
    for (int j = 0; j < J0; ++j) {
      start0.push_back(cluster_start0_[j] - 1);
      end0.push_back(cluster_end0_[j]);
    }
  }
};

// [[Rcpp::export]]
SEXP setdata_03_lmm_CM(
    const Eigen::Map<Eigen::VectorXd> y1,
    const Eigen::Map<Eigen::VectorXd> t1,
    const Eigen::Map<Eigen::VectorXd> Ti1,
    const Eigen::Map<Eigen::VectorXd> fTi1,
    const Rcpp::IntegerVector cluster_start1,
    const Rcpp::IntegerVector cluster_end1,
    const Eigen::Map<Eigen::VectorXd> y0,
    const Eigen::Map<Eigen::VectorXd> t0,
    const Eigen::Map<Eigen::VectorXd> Ti_til0,
    const Rcpp::IntegerVector cluster_start0,
    const Rcpp::IntegerVector cluster_end0) {
  Rcpp::XPtr<LmmCM03Data> ptr(
    new LmmCM03Data(
      y1, t1, Ti1, fTi1, cluster_start1, cluster_end1,
      y0, t0, Ti_til0, cluster_start0, cluster_end0
    ),
    true
  );
  return ptr;
}

// ------------------------------------------------------------
// Quadrature data kept separate from longitudinal data, so R can
// rebuild it at perturbed phi without recompiling or rebuilding the
// fixed longitudinal-data XPtr.
// ------------------------------------------------------------
struct LmmCM03Quad {
  MatrixXd Ti_mat;
  MatrixXd lw_mat;
  MatrixXd fTi_mat;
  int K;
  int J0;

  LmmCM03Quad(const Map<MatrixXd>& Ti_mat_,
              const Map<MatrixXd>& lw_mat_,
              const Map<MatrixXd>& fTi_mat_,
              const Map<VectorXd>& Ti_til0_)
    : Ti_mat(Ti_mat_), lw_mat(lw_mat_), fTi_mat(fTi_mat_),
      K(Ti_mat_.rows()), J0(Ti_mat_.cols()) {
    if (lw_mat_.rows() != K || fTi_mat_.rows() != K ||
        lw_mat_.cols() != J0 || fTi_mat_.cols() != J0 ||
        Ti_til0_.size() != J0) {
      Rcpp::stop("quadrature matrix dimensions are not aligned");
    }

    for (int j = 0; j < J0; ++j) {
      for (int k = 0; k < K; ++k) {
        if (!(Ti_mat(k, j) > Ti_til0_(j))) {
          Rcpp::stop("quadrature node is not above its truncation boundary");
        }
      }
    }
  }
};

// [[Rcpp::export]]
SEXP setquad_03_lmm_CM(
    const Eigen::Map<Eigen::MatrixXd> Ti_mat,
    const Eigen::Map<Eigen::MatrixXd> lw_mat,
    const Eigen::Map<Eigen::MatrixXd> fTi_mat,
    const Eigen::Map<Eigen::VectorXd> Ti_til0) {
  Rcpp::XPtr<LmmCM03Quad> ptr(
    new LmmCM03Quad(Ti_mat, lw_mat, fTi_mat, Ti_til0),
    true
  );
  return ptr;
}

// ------------------------------------------------------------
// Copied random-intercept cluster likelihood core. rho is the
// variance ratio sigma_b^2 / sigma^2, not sigma_b itself.
// ------------------------------------------------------------
template <typename T>
T cluster_ll_from_summaries_03(const int n,
                               const T& sum_r2,
                               const T& sse_val,
                               const T& sigma_b,
                               const T& sigma) {
  T rho = stan::math::square(sigma_b / sigma);
  T d0 = -T(n) * stan::math::log(sigma);
  T ssr = -0.5 / stan::math::square(sigma) *
    (sum_r2 + T(n) * rho * sse_val) / (T(1.0) + T(n) * rho);

  return d0 - 0.5 * stan::math::log1p(T(n) * rho) + ssr;
}

template <typename T>
T cluster_loglik_residual_vec_03(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& r,
    const T& sigma_b,
    const T& sigma) {
  const int n = r.size();

  T sum_r = 0.0;
  for (int k = 0; k < n; ++k) {
    sum_r += r(k);
  }
  T mean_r = sum_r / T(n);

  T sum_r2 = 0.0;
  T sse_val = 0.0;
  for (int k = 0; k < n; ++k) {
    sum_r2 += stan::math::square(r(k));
    sse_val += stan::math::square(r(k) - mean_r);
  }

  return cluster_ll_from_summaries_03(n, sum_r2, sse_val, sigma_b, sigma);
}

// [[Rcpp::export]]
double cluster_ll_cpp(const Eigen::Map<Eigen::VectorXd> r,
                      const double sigma_b,
                      const double sigma) {
  VectorXd r_vec = r;
  return cluster_loglik_residual_vec_03<double>(r_vec, sigma_b, sigma);
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
observed_residual_vec_03(const LmmCM03Data& d,
                         const int j,
                         const T& xi,
                         const T& beta0,
                         const T& beta1,
                         const T& beta_T,
                         const T& beta_g) {
  const int s = d.start1[j];
  const int e = d.end1[j];
  const int n = e - s;
  Eigen::Matrix<T, Eigen::Dynamic, 1> r(n);

  const double Ti = d.Ti1(j);
  const double fTi = d.fTi1(j);
  for (int row = s; row < e; ++row) {
    T g = stan::math::exp(xi * (d.t1(row) - Ti));
    T mu = beta0 + beta1 * d.t1(row) + beta_T * fTi + beta_g * g;
    r(row - s) = d.y1(row) - mu;
  }

  return r;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
missing_residual_vec_03(const LmmCM03Data& d,
                        const LmmCM03Quad& q,
                        const int j,
                        const int node,
                        const T& xi,
                        const T& beta0,
                        const T& beta1,
                        const T& beta_T,
                        const T& beta_g) {
  const int s = d.start0[j];
  const int e = d.end0[j];
  const int n = e - s;
  Eigen::Matrix<T, Eigen::Dynamic, 1> r(n);

  const double Ti = q.Ti_mat(node, j);
  const double fTi = q.fTi_mat(node, j);
  for (int row = s; row < e; ++row) {
    T g = stan::math::exp(xi * (d.t0(row) - Ti));
    T mu = beta0 + beta1 * d.t0(row) + beta_T * fTi + beta_g * g;
    r(row - s) = d.y0(row) - mu;
  }

  return r;
}

template <typename T>
T observed_cluster_loglik_03(const LmmCM03Data& d,
                             const int j,
                             const T& xi,
                             const T& beta0,
                             const T& beta1,
                             const T& beta_T,
                             const T& beta_g,
                             const T& sigma_b,
                             const T& sigma) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> r =
    observed_residual_vec_03(d, j, xi, beta0, beta1, beta_T, beta_g);
  return cluster_loglik_residual_vec_03<T>(r, sigma_b, sigma);
}

template <typename T>
T missing_cluster_loglik_03(const LmmCM03Data& d,
                            const LmmCM03Quad& q,
                            const int j,
                            const T& xi,
                            const T& beta0,
                            const T& beta1,
                            const T& beta_T,
                            const T& beta_g,
                            const T& sigma_b,
                            const T& sigma) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> log_terms(q.K);

  for (int node = 0; node < q.K; ++node) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> r =
      missing_residual_vec_03(d, q, j, node, xi, beta0, beta1, beta_T, beta_g);
    log_terms(node) =
      cluster_loglik_residual_vec_03<T>(r, sigma_b, sigma) + q.lw_mat(node, j);
  }

  return stan::math::log_sum_exp(log_terms);
}

// ------------------------------------------------------------
// Component per-cluster vectors. These are the single place each
// component's math lives.
// ------------------------------------------------------------
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
observed_cluster_loglik_vec_03(const LmmCM03Data& d,
                               const T& xi,
                               const T& beta0,
                               const T& beta1,
                               const T& beta_T,
                               const T& beta_g,
                               const T& sigma_b,
                               const T& sigma) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> out(d.J1);
  for (int j = 0; j < d.J1; ++j) {
    out(j) = observed_cluster_loglik_03(
      d, j, xi, beta0, beta1, beta_T, beta_g, sigma_b, sigma
    );
  }
  return out;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
missing_cluster_loglik_vec_03(const LmmCM03Data& d,
                              const LmmCM03Quad& q,
                              const T& xi,
                              const T& beta0,
                              const T& beta1,
                              const T& beta_T,
                              const T& beta_g,
                              const T& sigma_b,
                              const T& sigma) {
  if (q.J0 != d.J0) {
    Rcpp::stop("quadrature columns do not match datX0 clusters");
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> out(d.J0);
  for (int j = 0; j < d.J0; ++j) {
    out(j) = missing_cluster_loglik_03(
      d, q, j, xi, beta0, beta1, beta_T, beta_g, sigma_b, sigma
    );
  }
  return out;
}

struct LogLikVec {
  const LmmCM03Data& d;
  const LmmCM03Quad& q;
  LogLikVec(const LmmCM03Data& d_, const LmmCM03Quad& q_) : d(d_), q(q_) {}

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    T xi = theta(0);
    T beta0 = theta(1);
    T beta1 = theta(2);
    T beta_T = theta(3);
    T beta_g = theta(4);
    T sigma_b = stan::math::exp(theta(5));
    T sigma = stan::math::exp(theta(6));

    Eigen::Matrix<T, Eigen::Dynamic, 1> out(d.J1 + d.J0);
    out.segment(0, d.J1) = observed_cluster_loglik_vec_03(
      d, xi, beta0, beta1, beta_T, beta_g, sigma_b, sigma
    );
    out.segment(d.J1, d.J0) = missing_cluster_loglik_vec_03(
      d, q, xi, beta0, beta1, beta_T, beta_g, sigma_b, sigma
    );
    return out;
  }
};

struct LogLik {
  const LmmCM03Data& d;
  const LmmCM03Quad& q;
  LogLik(const LmmCM03Data& d_, const LmmCM03Quad& q_) : d(d_), q(q_) {}

  template <typename T>
  T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    LogLikVec vec_fn(d, q);
    return stan::math::sum(vec_fn(theta));
  }
};

// ------------------------------------------------------------
// Exported double wrappers for component tests.
// ------------------------------------------------------------

// [[Rcpp::export]]
Eigen::VectorXd observed_cluster_ll_cpp(
    SEXP data_ptr_sexp,
    const Eigen::Map<Eigen::VectorXd> theta_in) {
  Rcpp::XPtr<LmmCM03Data> data_ptr(data_ptr_sexp);
  const LmmCM03Data& d = *data_ptr;
  VectorXd theta = theta_in;
  return observed_cluster_loglik_vec_03<double>(
    d, theta(0), theta(1), theta(2), theta(3), theta(4),
    std::exp(theta(5)), std::exp(theta(6))
  );
}

// [[Rcpp::export]]
Eigen::VectorXd missing_cluster_ll_cpp(
    SEXP data_ptr_sexp,
    SEXP quad_ptr_sexp,
    const Eigen::Map<Eigen::VectorXd> theta_in) {
  Rcpp::XPtr<LmmCM03Data> data_ptr(data_ptr_sexp);
  Rcpp::XPtr<LmmCM03Quad> quad_ptr(quad_ptr_sexp);
  const LmmCM03Data& d = *data_ptr;
  const LmmCM03Quad& q = *quad_ptr;
  VectorXd theta = theta_in;
  return missing_cluster_loglik_vec_03<double>(
    d, q, theta(0), theta(1), theta(2), theta(3), theta(4),
    std::exp(theta(5)), std::exp(theta(6))
  );
}

// ------------------------------------------------------------
// Main value + gradient endpoint.
// ------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List ll_R_cpp(const Eigen::Map<Eigen::VectorXd> theta_in,
                    SEXP data_ptr_sexp,
                    SEXP quad_ptr_sexp) {
  Rcpp::XPtr<LmmCM03Data> data_ptr(data_ptr_sexp);
  Rcpp::XPtr<LmmCM03Quad> quad_ptr(quad_ptr_sexp);
  LogLik ll(*data_ptr, *quad_ptr);
  double fx;
  VectorXd grad_fx;
  VectorXd theta = theta_in;
  stan::math::gradient(ll, theta, fx, grad_fx);

  return Rcpp::List::create(
    Rcpp::Named("value") = fx,
    Rcpp::Named("gradient") = grad_fx
  );
}

// ------------------------------------------------------------
// Per-cluster value endpoint, double only, no autodiff.
// ------------------------------------------------------------
// [[Rcpp::export]]
Eigen::VectorXd ll_R_cpp_vec(const Eigen::Map<Eigen::VectorXd> theta_in,
                             SEXP data_ptr_sexp,
                             SEXP quad_ptr_sexp) {
  Rcpp::XPtr<LmmCM03Data> data_ptr(data_ptr_sexp);
  Rcpp::XPtr<LmmCM03Quad> quad_ptr(quad_ptr_sexp);
  LogLikVec ll(*data_ptr, *quad_ptr);
  VectorXd theta = theta_in;
  return ll.template operator()<double>(theta);
}
