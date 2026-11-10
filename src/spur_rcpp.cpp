// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <vector>

#ifdef __APPLE__
extern "C" void vvexp(double* y, const double* x, const int* n);
#endif

using Rcpp::Dimension;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

inline int resolve_spur_threads(int requested = -1) {
  if (requested > 0) return requested;
  const char* cran = std::getenv("_R_CHECK_LIMIT_CORES_");
  if (cran != nullptr && *cran != '\0') return 2;
  const char* env = std::getenv("SPUR_NUM_THREADS");
  if (env != nullptr && *env != '\0') {
    char* end;
    long val = std::strtol(env, &end, 10);
    if (end != env && *end == '\0' && val > 0) return static_cast<int>(val);
  }
  env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  if (env != nullptr && *env != '\0') {
    char* end;
    long val = std::strtol(env, &end, 10);
    if (end != env && *end == '\0' && val > 0) return static_cast<int>(val);
  }
  return -1;
}

template <typename WorkerType>
inline void spur_parallel_for(std::size_t begin, std::size_t end, WorkerType& worker, std::size_t grainSize = 1, int nthreads = -1) {
  const int nth = resolve_spur_threads(nthreads);
  if (nth == 1) {
    worker(begin, end);
  } else {
    RcppParallel::parallelFor(begin, end, worker, grainSize, nth);
  }
}

inline double fast_neg_half_pow_15(double x) {
  const double x2 = x * x;
  const double x4 = x2 * x2;
  const double x7 = x4 * x2 * x;
  return 1.0 / (x7 * std::sqrt(x));
}

inline double fast_neg_half_pow(double x, int q) {
  if (q == 15) {
    return fast_neg_half_pow_15(x);
  }
  int m = q >> 1;
  double acc = 1.0;
  double base = x;
  while (m > 0) {
    if (m & 1) acc *= base;
    base *= base;
    m >>= 1;
  }
  if (q & 1) {
    acc *= std::sqrt(x);
  }
  return 1.0 / acc;
}

// Compute R type-7 quantiles in-place on `buf` for `probs` (length n_probs)
void compute_quantiles_type7(std::vector<double>& buf,
                             const double* probs,
                             int n_probs,
                             double* out) {
  const int n = static_cast<int>(buf.size());
  if (n == 0) {
    for (int i = 0; i < n_probs; ++i) out[i] = NA_REAL;
    return;
  }

  struct Req {
    int lo_idx;
    int hi_idx;
    double h;
  };
  std::vector<Req> reqs(n_probs);
  std::vector<int> needed;
  needed.reserve(2 * n_probs);

  for (int i = 0; i < n_probs; ++i) {
    double p = probs[i];
    double index = 1.0 + static_cast<double>(n - 1) * p;
    int lo = static_cast<int>(std::floor(index));
    int hi = static_cast<int>(std::ceil(index));
    if (lo < 1) lo = 1;
    if (hi < 1) hi = 1;
    if (lo > n) lo = n;
    if (hi > n) hi = n;
    reqs[i] = {lo - 1, hi - 1, index - static_cast<double>(lo)};
    needed.push_back(lo - 1);
    if (hi != lo) needed.push_back(hi - 1);
  }

  std::sort(needed.begin(), needed.end());
  needed.erase(std::unique(needed.begin(), needed.end()), needed.end());

  int start = 0;
  for (int idx : needed) {
    if (idx == start) {
      auto min_it = std::min_element(buf.begin() + start, buf.end());
      std::iter_swap(buf.begin() + start, min_it);
    } else {
      std::nth_element(buf.begin() + start, buf.begin() + idx, buf.end());
    }
    start = idx + 1;
  }

  for (int i = 0; i < n_probs; ++i) {
    double qs = buf[reqs[i].lo_idx];
    double x_hi = buf[reqs[i].hi_idx];
    if (reqs[i].h > 0.0 && x_hi != qs) {
      qs = (1.0 - reqs[i].h) * qs + reqs[i].h * x_hi;
    }
    out[i] = qs;
  }
}

void demean_matrix_inplace(double* mat, int n) {
  if (n <= 0) return;
  const double inv_n = 1.0 / static_cast<double>(n);
  std::vector<double> row_sum(n, 0.0);

  for (int j = 0; j < n; ++j) {
    double* col = mat + static_cast<size_t>(j) * n;
    double csum = 0.0;
    for (int i = 0; i < n; ++i) {
      csum += col[i];
    }
    const double cmean = csum * inv_n;
    for (int i = 0; i < n; ++i) {
      double v = col[i] - cmean;
      col[i] = v;
      row_sum[i] += v;
    }
  }

  for (int i = 0; i < n; ++i) {
    row_sum[i] *= inv_n;
  }
  for (int j = 0; j < n; ++j) {
    double* col = mat + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      col[i] -= row_sum[i];
    }
  }
}

void symmetrize_lower_to_upper_blocked(double* S, int n) {
  constexpr int BLK = 64;
  for (int j0 = 0; j0 < n; j0 += BLK) {
    const int j_end = std::min(n, j0 + BLK);
    for (int i0 = j0; i0 < n; i0 += BLK) {
      const int i_end = std::min(n, i0 + BLK);
      if (i0 == j0) {
        for (int j = j0; j < j_end; ++j) {
          const double* src_col = S + static_cast<size_t>(j) * n;
          for (int i = j + 1; i < i_end; ++i) {
            S[j + static_cast<size_t>(i) * n] = src_col[i];
          }
        }
      } else {
        for (int j = j0; j < j_end; ++j) {
          const double* src_col = S + static_cast<size_t>(j) * n;
          for (int i = i0; i < i_end; ++i) {
            S[j + static_cast<size_t>(i) * n] = src_col[i];
          }
        }
      }
    }
  }
}

void exp_lower_triangle_inplace(const double* distmat, double* S, int n, double neg_c) {
  for (int i = 0; i < n; ++i) {
    const double* d_col = distmat + static_cast<size_t>(i) * n;
    double* s_col = S + static_cast<size_t>(i) * n;
    const double d_ii = d_col[i];
    s_col[i] = (d_ii == 0.0) ? 1.0 : std::exp(neg_c * d_ii);
    const int len = n - 1 - i;
    if (len > 0) {
      const double* src = d_col + i + 1;
      double* dst = s_col + i + 1;
#ifdef __APPLE__
      for (int k = 0; k < len; ++k) {
        dst[k] = neg_c * src[k];
      }
      vvexp(dst, dst, &len);
#else
      for (int k = 0; k < len; ++k) {
        dst[k] = std::exp(neg_c * src[k]);
      }
#endif
    }
  }
}

// Parallel worker for great-circle distance lower triangle
struct DistLatLonWorker : public RcppParallel::Worker {
  const double* lat;
  const double* lon;
  const double* cos_lat;
  double* out_ptr;
  int n;

  DistLatLonWorker(const double* lat, const double* lon, const double* cos_lat, double* out_ptr, int n)
      : lat(lat), lon(lon), cos_lat(cos_lat), out_ptr(out_ptr), n(n) {}

  void operator()(std::size_t begin, std::size_t end) override {
    const double c = M_PI / 180.0;
    const double inv_pi = 1.0 / M_PI;
    for (std::size_t i_sz = begin; i_sz < end; ++i_sz) {
      const int i = static_cast<int>(i_sz);
      double* out_col = out_ptr + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      const double lat_i = lat[i];
      const double lon_i = lon[i];
      const double cos_i = cos_lat[i];
      for (int j = i + 1; j < n; ++j) {
        const double d1 = 0.5 * c * (lat[j] - lat_i);
        const double d2 = 0.5 * c * (lon[j] - lon_i);
        const double sd1 = std::sin(d1);
        const double sd2 = std::sin(d2);
        const double inside = sd1 * sd1 + cos_i * (cos_lat[j] * (sd2 * sd2));
        out_col[j] = std::asin(std::sqrt(inside)) * inv_pi;
      }
    }
  }
};

// Parallel worker for projected covariance grid chunks [m_start, m_end)
struct OmGridChunkWorker : public RcppParallel::Worker {
  const double* d_ptr_mat;
  const double* r_ptr;
  const double* c_ptr;
  const uint8_t* step_mode;
  const double* step_dc;
  const std::vector<std::pair<int, int>>& chunks;
  double* out_ptr;
  int n;
  int q;

  OmGridChunkWorker(const double* d_ptr_mat,
                    const double* r_ptr,
                    const double* c_ptr,
                    const uint8_t* step_mode,
                    const double* step_dc,
                    const std::vector<std::pair<int, int>>& chunks,
                    double* out_ptr,
                    int n,
                    int q)
      : d_ptr_mat(d_ptr_mat),
        r_ptr(r_ptr),
        c_ptr(c_ptr),
        step_mode(step_mode),
        step_dc(step_dc),
        chunks(chunks),
        out_ptr(out_ptr),
        n(n),
        q(q) {}

  void operator()(std::size_t c_begin, std::size_t c_end) override {
    constexpr int BLK = 256;
    const arma::mat R_mat(const_cast<double*>(r_ptr), n, q, false, true);

    std::vector<double> D_blk_buf(BLK * BLK);
    std::vector<double> S_blk_buf(BLK * BLK);
    std::vector<double> E_blk_buf(BLK * BLK);

    for (std::size_t ci = c_begin; ci < c_end; ++ci) {
      const int m0 = chunks[ci].first;
      const int m1 = chunks[ci].second;
      const int m_len = m1 - m0;

      arma::cube SR_chunk(n, q, m_len, arma::fill::zeros);

      for (int j0 = 0; j0 < n; j0 += BLK) {
        const int bj = std::min(BLK, n - j0);
        for (int i0 = j0; i0 < n; i0 += BLK) {
          const int bi = std::min(BLK, n - i0);
          const int len = bi * bj;

          for (int j = 0; j < bj; ++j) {
            const double* d_col = d_ptr_mat + static_cast<size_t>(j0 + j) * n + i0;
            std::memcpy(D_blk_buf.data() + static_cast<size_t>(j) * bi, d_col, bi * sizeof(double));
          }

          arma::mat S_blk(S_blk_buf.data(), bi, bj, false, true);
          const arma::subview<double> R_i = R_mat.rows(i0, i0 + bi - 1);
          const arma::subview<double> R_j = R_mat.rows(j0, j0 + bj - 1);

          for (int m = m0; m < m1; ++m) {
            const uint8_t mode = step_mode[m];
            if (mode == 2) {
              double* s_ptr = S_blk_buf.data();
              const double* e_ptr = E_blk_buf.data();
              for (int k = 0; k < len; ++k) {
                s_ptr[k] *= e_ptr[k];
              }
            } else if (mode == 1) {
              const double neg_dc = -step_dc[m];
              const double* d_ptr = D_blk_buf.data();
              double* e_ptr = E_blk_buf.data();
              double* s_ptr = S_blk_buf.data();
#ifdef __APPLE__
              for (int k = 0; k < len; ++k) {
                e_ptr[k] = neg_dc * d_ptr[k];
              }
              vvexp(e_ptr, e_ptr, &len);
#else
              for (int k = 0; k < len; ++k) {
                e_ptr[k] = std::exp(neg_dc * d_ptr[k]);
              }
#endif
              if (i0 == j0) {
                for (int d = 0; d < bi; ++d) {
                  if (d_ptr[d + d * bi] == 0.0) e_ptr[d + d * bi] = 1.0;
                }
              }
              for (int k = 0; k < len; ++k) {
                s_ptr[k] *= e_ptr[k];
              }
            } else {
              const double neg_c = -c_ptr[m];
              const double* d_ptr = D_blk_buf.data();
              double* s_ptr = S_blk_buf.data();
#ifdef __APPLE__
              for (int k = 0; k < len; ++k) {
                s_ptr[k] = neg_c * d_ptr[k];
              }
              vvexp(s_ptr, s_ptr, &len);
#else
              for (int k = 0; k < len; ++k) {
                s_ptr[k] = std::exp(neg_c * d_ptr[k]);
              }
#endif
              if (i0 == j0) {
                for (int d = 0; d < bi; ++d) {
                  if (d_ptr[d + d * bi] == 0.0) s_ptr[d + d * bi] = 1.0;
                }
              }
            }

            arma::mat& SR_m = SR_chunk.slice(m - m0);
            if (i0 == j0) {
              SR_m.rows(i0, i0 + bi - 1) += S_blk * R_i;
            } else {
              SR_m.rows(i0, i0 + bi - 1) += S_blk * R_j;
              SR_m.rows(j0, j0 + bj - 1) += S_blk.t() * R_i;
            }
          }
        }
      }

      for (int m = m0; m < m1; ++m) {
        arma::mat om_m(out_ptr + static_cast<size_t>(m) * q * q, q, q, false, true);
        om_m = R_mat.t() * SR_chunk.slice(m - m0);
      }
    }
  }
};

// Parallel worker for Monte Carlo quadratic form tiles in rcpp_pow_qf_eval
struct PowQfTileWorker : public RcppParallel::Worker {
  const double* ho_ptr;
  const double* ha_ptr;
  const double* e_ptr;
  const double* qe_ptr;
  double* lr_o_ptr;
  double* lr_a_ptr;
  int q;
  int nrep;
  int tile_size;

  PowQfTileWorker(const double* ho_ptr,
                  const double* ha_ptr,
                  const double* e_ptr,
                  const double* qe_ptr,
                  double* lr_o_ptr,
                  double* lr_a_ptr,
                  int q,
                  int nrep,
                  int tile_size)
      : ho_ptr(ho_ptr),
        ha_ptr(ha_ptr),
        e_ptr(e_ptr),
        qe_ptr(qe_ptr),
        lr_o_ptr(lr_o_ptr),
        lr_a_ptr(lr_a_ptr),
        q(q),
        nrep(nrep),
        tile_size(tile_size) {}

  void operator()(std::size_t t_begin, std::size_t t_end) override {
    const arma::mat ho_mat(const_cast<double*>(ho_ptr), q, q, false, true);
    const arma::mat ha_mat(const_cast<double*>(ha_ptr), q, q, false, true);
    arma::mat ya_tile(q, tile_size);
    arma::mat yo_tile(q, tile_size);

    for (std::size_t t = t_begin; t < t_end; ++t) {
      const int col0 = static_cast<int>(t) * tile_size;
      const int b = std::min(tile_size, nrep - col0);
      const arma::mat e_tile(const_cast<double*>(e_ptr + static_cast<size_t>(col0) * q), q, b, false, true);

      arma::mat ya_sub(ya_tile.memptr(), q, b, false, true);
      arma::mat yo_sub(yo_tile.memptr(), q, b, false, true);
      ya_sub = ho_mat * e_tile;
      yo_sub = ha_mat * e_tile;

      for (int j = 0; j < b; ++j) {
        const double* ya_col = ya_sub.colptr(j);
        const double* yo_col = yo_sub.colptr(j);
        double qa_o = 0.0;
        double qo_a = 0.0;
        for (int i = 0; i < q; ++i) {
          qa_o += ya_col[i] * ya_col[i];
          qo_a += yo_col[i] * yo_col[i];
        }
        const double qe_k = qe_ptr[col0 + j];
        lr_o_ptr[col0 + j] = qe_k / qa_o;
        lr_a_ptr[col0 + j] = qo_a / qe_k;
      }
    }
  }
};

// Parallel worker for Monte Carlo LR grid (across ir = 0 .. n_rho - 1)
struct McLrGridWorker : public RcppParallel::Worker {
  const double* ch_om_ho_ptr;
  const double* ch_omi_ho_ptr;
  const double* ch_omi_ha_ptr;
  const double* emat_ptr;
  const double* lr_stat_ptr;
  const double* probs_ptr;
  double* cvalue_ptr;
  double* pvalue_ptr;
  int q;
  int nrep;
  int n_y;
  int n_probs;
  int n_rho;

  McLrGridWorker(const double* ch_om_ho_ptr,
                 const double* ch_omi_ho_ptr,
                 const double* ch_omi_ha_ptr,
                 const double* emat_ptr,
                 const double* lr_stat_ptr,
                 const double* probs_ptr,
                 double* cvalue_ptr,
                 double* pvalue_ptr,
                 int q,
                 int nrep,
                 int n_y,
                 int n_probs,
                 int n_rho)
      : ch_om_ho_ptr(ch_om_ho_ptr),
        ch_omi_ho_ptr(ch_omi_ho_ptr),
        ch_omi_ha_ptr(ch_omi_ha_ptr),
        emat_ptr(emat_ptr),
        lr_stat_ptr(lr_stat_ptr),
        probs_ptr(probs_ptr),
        cvalue_ptr(cvalue_ptr),
        pvalue_ptr(pvalue_ptr),
        q(q),
        nrep(nrep),
        n_y(n_y),
        n_probs(n_probs),
        n_rho(n_rho) {}

  void operator()(std::size_t ir_begin, std::size_t ir_end) override {
    constexpr int TILE = 512;
    const arma::mat ch_omi_ho_mat(const_cast<double*>(ch_omi_ho_ptr), q, q, false, true);
    const arma::mat ch_omi_ha_mat(const_cast<double*>(ch_omi_ha_ptr), q, q, false, true);

    arma::mat A_ho(q, q);
    arma::mat A_ha(q, q);
    arma::mat y_ho_ho_tile(q, TILE);
    arma::mat y_ho_ha_tile(q, TILE);
    std::vector<double> lr_ho(nrep);
    std::vector<double> cv_row(n_probs);
    std::vector<int> gt_counts(n_y);

    for (std::size_t ir_sz = ir_begin; ir_sz < ir_end; ++ir_sz) {
      const int ir = static_cast<int>(ir_sz);
      const arma::mat ch_om_ho(const_cast<double*>(ch_om_ho_ptr + static_cast<size_t>(ir) * q * q), q, q, false, true);
      std::fill(gt_counts.begin(), gt_counts.end(), 0);

      A_ho = ch_omi_ho_mat * ch_om_ho.t();
      A_ha = ch_omi_ha_mat * ch_om_ho.t();

      for (int col0 = 0; col0 < nrep; col0 += TILE) {
        const int b = std::min(TILE, nrep - col0);
        const arma::mat emat_tile(const_cast<double*>(emat_ptr + static_cast<size_t>(col0) * q), q, b, false, true);

        arma::mat y_ho_sub(y_ho_ho_tile.memptr(), q, b, false, true);
        arma::mat y_ha_sub(y_ho_ha_tile.memptr(), q, b, false, true);
        y_ho_sub = A_ho * emat_tile;
        y_ha_sub = A_ha * emat_tile;

        for (int j = 0; j < b; ++j) {
          const double* col_ho = y_ho_sub.colptr(j);
          const double* col_ha = y_ha_sub.colptr(j);
          double q_ho = 0.0;
          double q_ha = 0.0;
          for (int i = 0; i < q; ++i) {
            q_ho += col_ho[i] * col_ho[i];
            q_ha += col_ha[i] * col_ha[i];
          }
          const double lr_val = q_ho / q_ha;
          lr_ho[col0 + j] = lr_val;
          for (int jy = 0; jy < n_y; ++jy) {
            if (lr_val > lr_stat_ptr[jy]) gt_counts[jy]++;
          }
        }
      }

      compute_quantiles_type7(lr_ho, probs_ptr, n_probs, cv_row.data());
      for (int jp = 0; jp < n_probs; ++jp) {
        cvalue_ptr[ir + static_cast<size_t>(jp) * n_rho] = cv_row[jp];
      }
      for (int jy = 0; jy < n_y; ++jy) {
        pvalue_ptr[ir + static_cast<size_t>(jy) * n_rho] =
            static_cast<double>(gt_counts[jy]) / static_cast<double>(nrep);
      }
    }
  }
};

// Parallel worker for c_ci Monte Carlo across null grid points i = 0 .. n_c_ho - 1
struct CCiMcWorker : public RcppParallel::Worker {
  const double* ch_om_ho_ptr;
  const double* ch_omi_ho_ptr;
  const double* const_den_ho_ptr;
  const double* const_den_ha_ptr;
  const double* emat_ptr;
  const double* emat_pow_ptr;
  const double* W_ha_ptr;
  const double* X_ptr;
  double* pv_ptr;
  double den_ha_avg_X;
  int q;
  int p;
  int nrep;
  int n_c_ha;

  CCiMcWorker(const double* ch_om_ho_ptr,
              const double* ch_omi_ho_ptr,
              const double* const_den_ho_ptr,
              const double* const_den_ha_ptr,
              const double* emat_ptr,
              const double* emat_pow_ptr,
              const double* W_ha_ptr,
              const double* X_ptr,
              double* pv_ptr,
              double den_ha_avg_X,
              int q,
              int p,
              int nrep,
              int n_c_ha)
      : ch_om_ho_ptr(ch_om_ho_ptr),
        ch_omi_ho_ptr(ch_omi_ho_ptr),
        const_den_ho_ptr(const_den_ho_ptr),
        const_den_ha_ptr(const_den_ha_ptr),
        emat_ptr(emat_ptr),
        emat_pow_ptr(emat_pow_ptr),
        W_ha_ptr(W_ha_ptr),
        X_ptr(X_ptr),
        pv_ptr(pv_ptr),
        den_ha_avg_X(den_ha_avg_X),
        q(q),
        p(p),
        nrep(nrep),
        n_c_ha(n_c_ha) {}

  void operator()(std::size_t i_begin, std::size_t i_end) override {
    constexpr int TILE = 512;
    const arma::mat W_ha_mat(const_cast<double*>(W_ha_ptr), n_c_ha, p, false, true);
    const double inv_n_ha = 1.0 / static_cast<double>(n_c_ha);

    arma::mat e_tile_buf(q, TILE);
    arma::mat U_tile_buf(p, TILE);
    arma::mat Q_tile_buf(n_c_ha, TILE);

    for (std::size_t i_sz = i_begin; i_sz < i_end; ++i_sz) {
      const int i = static_cast<int>(i_sz);
      const arma::mat ch_null(const_cast<double*>(ch_om_ho_ptr + static_cast<size_t>(i) * q * q), q, q, false, true);
      const double* ch_omi = ch_omi_ho_ptr + static_cast<size_t>(i) * q * q;
      const double const_den = const_den_ho_ptr[i];

      double q_X_ho = 0.0;
      for (int r = 0; r < q; ++r) {
        double acc = 0.0;
        for (int c = 0; c < q; ++c) {
          acc += ch_omi[r + static_cast<size_t>(c) * q] * X_ptr[c];
        }
        q_X_ho += acc * acc;
      }
      const double den_ho_X = const_den * fast_neg_half_pow(q_X_ho, q);
      const double lr_X = den_ha_avg_X / den_ho_X;

      int count_gt = 0;

      for (int col0 = 0; col0 < nrep; col0 += TILE) {
        const int b = std::min(TILE, nrep - col0);
        const arma::mat emat_tile(const_cast<double*>(emat_ptr + static_cast<size_t>(col0) * q), q, b, false, true);

        arma::mat e_sub(e_tile_buf.memptr(), q, b, false, true);
        arma::mat U_sub(U_tile_buf.memptr(), p, b, false, true);
        arma::mat Q_sub(Q_tile_buf.memptr(), n_c_ha, b, false, true);

        e_sub = ch_null.t() * emat_tile;

        for (int col = 0; col < b; ++col) {
          const double* e_col = e_sub.colptr(col);
          double* u_col = U_sub.colptr(col);
          int k = 0;
          for (int cb = 0; cb < q; ++cb) {
            const double eb = e_col[cb];
            for (int ca = 0; ca <= cb; ++ca) {
              u_col[k++] = e_col[ca] * eb;
            }
          }
        }

        Q_sub = W_ha_mat * U_sub;

        for (int col = 0; col < b; ++col) {
          const double den_ho = const_den * emat_pow_ptr[col0 + col];
          const double* q_ha_col = Q_sub.colptr(col);
          double sum_ha = 0.0;
          if (q == 15) {
            for (int j = 0; j < n_c_ha; ++j) {
              sum_ha += const_den_ha_ptr[j] * fast_neg_half_pow_15(q_ha_col[j]);
            }
          } else {
            for (int j = 0; j < n_c_ha; ++j) {
              sum_ha += const_den_ha_ptr[j] * fast_neg_half_pow(q_ha_col[j], q);
            }
          }
          const double den_ha_avg = sum_ha * inv_n_ha;
          const double lr = den_ha_avg / den_ho;
          if (lr > lr_X) count_gt++;
        }
      }

      pv_ptr[i] = static_cast<double>(count_gt) / static_cast<double>(nrep);
    }
  }
};

}  // namespace

// [[Rcpp::export]]
NumericMatrix rcpp_distmat_lat_lon(NumericMatrix s, int nthreads = -1) {
  const int n = s.nrow();
  const double* s_ptr = s.begin();
  const double* lat = s_ptr;
  const double* lon = s_ptr + n;
  const double c = M_PI / 180.0;

  std::vector<double> cos_lat(n);
  for (int i = 0; i < n; ++i) {
    cos_lat[i] = std::cos(c * lat[i]);
  }

  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  DistLatLonWorker worker(lat, lon, cos_lat.data(), out_ptr, n);
  if (n >= 512) {
    spur_parallel_for(0, n, worker, 64, nthreads);
  } else {
    worker(0, n);
  }
  symmetrize_lower_to_upper_blocked(out_ptr, n);
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_distmat_euclidean(NumericMatrix s) {
  const int n = s.nrow();
  const int d = s.ncol();
  const double* s_ptr = s.begin();

  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  if (d == 2) {
    const double* s0 = s_ptr;
    const double* s1 = s_ptr + n;
    for (int i = 0; i < n; ++i) {
      double* out_col = out_ptr + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      const double x0 = s0[i];
      const double y0 = s1[i];
      for (int j = i + 1; j < n; ++j) {
        const double dx = s0[j] - x0;
        const double dy = s1[j] - y0;
        out_col[j] = std::sqrt(dx * dx + dy * dy);
      }
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double* out_col = out_ptr + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      for (int j = i + 1; j < n; ++j) {
        double sum_sq = 0.0;
        for (int k = 0; k < d; ++k) {
          const double diff = s_ptr[j + static_cast<size_t>(k) * n] - s_ptr[i + static_cast<size_t>(k) * n];
          sum_sq += diff * diff;
        }
        out_col[j] = std::sqrt(sum_sq);
      }
    }
  }
  symmetrize_lower_to_upper_blocked(out_ptr, n);
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_demean_sigma(NumericMatrix sigma) {
  const int n = sigma.nrow();
  NumericMatrix out(Dimension(n, n));
  std::memcpy(out.begin(), sigma.begin(), static_cast<size_t>(n) * n * sizeof(double));
  demean_matrix_inplace(out.begin(), n);
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_get_sigma_dm(NumericMatrix distmat, double c) {
  const int n = distmat.nrow();
  NumericMatrix out(Dimension(n, n));
  exp_lower_triangle_inplace(distmat.begin(), out.begin(), n, -c);
  symmetrize_lower_to_upper_blocked(out.begin(), n);
  demean_matrix_inplace(out.begin(), n);
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_get_sigma_lbm(NumericMatrix distmat) {
  const int n = distmat.nrow();
  const double* d_ptr = distmat.begin();
  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  for (int j = 0; j < n; ++j) {
    const double d_0j = d_ptr[static_cast<size_t>(j) * n];
    const double* d_col = d_ptr + static_cast<size_t>(j) * n;
    double* out_col = out_ptr + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      out_col[i] = 0.5 * (d_ptr[i] + d_0j - d_col[i]);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_get_sigma_lbm_dm(NumericMatrix distmat) {
  NumericMatrix out = rcpp_get_sigma_lbm(distmat);
  demean_matrix_inplace(out.begin(), distmat.nrow());
  return out;
}

// [[Rcpp::export]]
NumericVector rcpp_get_cbar_vec(NumericVector rhobar_vec, NumericMatrix distmat) {
  const int n_rho = rhobar_vec.size();
  const int n = distmat.nrow();
  const double* d_ptr = distmat.begin();

  const size_t m = static_cast<size_t>(n) * (n - 1) / 2;

  struct EvalPoint {
    double c;
    double v;
  };
  struct HistBin {
    double count;
    double mean;
    double lo;
    double hi;
    double w_lo;
    double w_hi;
  };
  static int cached_n = 0;
  static std::vector<double> cached_vd_raw;
  static std::vector<HistBin> cached_bins;
  static std::vector<EvalPoint> history;
  static std::unordered_map<uint64_t, double> cbar_cache;

  bool match = (cached_n == n && cached_vd_raw.size() == m);
  if (match && m > 0) {
    size_t pos = 0;
    for (int j = 0; j < n - 1 && match; ++j) {
      const double* col = d_ptr + static_cast<size_t>(j) * n + (j + 1);
      const size_t len = static_cast<size_t>(n - 1 - j);
      if (std::memcmp(cached_vd_raw.data() + pos, col, len * sizeof(double)) != 0) {
        match = false;
        break;
      }
      pos += len;
    }
  }

  if (!match) {
    cached_n = n;
    cached_vd_raw.resize(m);
    double min_d = std::numeric_limits<double>::infinity();
    double max_d = -std::numeric_limits<double>::infinity();
    size_t pos = 0;
    for (int j = 0; j < n - 1; ++j) {
      const double* col = d_ptr + static_cast<size_t>(j) * n + (j + 1);
      const size_t len = static_cast<size_t>(n - 1 - j);
      std::memcpy(cached_vd_raw.data() + pos, col, len * sizeof(double));
      for (size_t k = 0; k < len; ++k) {
        const double d = col[k];
        if (d < min_d) min_d = d;
        if (d > max_d) max_d = d;
      }
      pos += len;
    }

    cached_bins.clear();
    if (m > 0) {
      constexpr int NBINS = 1024;
      size_t b_count[NBINS] = {0};
      double b_sum[NBINS] = {0.0};
      double b_min[NBINS];
      double b_max[NBINS];
      for (int b = 0; b < NBINS; ++b) {
        b_min[b] = std::numeric_limits<double>::infinity();
        b_max[b] = -std::numeric_limits<double>::infinity();
      }
      const double span = max_d - min_d;
      const double inv_span = (span > 0.0) ? (static_cast<double>(NBINS) / span) : 0.0;
      const double* raw = cached_vd_raw.data();
      for (size_t k = 0; k < m; ++k) {
        const double d = raw[k];
        int b = static_cast<int>((d - min_d) * inv_span);
        if (b < 0) b = 0;
        else if (b >= NBINS) b = NBINS - 1;
        b_count[b]++;
        b_sum[b] += d;
        if (d < b_min[b]) b_min[b] = d;
        if (d > b_max[b]) b_max[b] = d;
      }
      cached_bins.reserve(NBINS);
      for (int b = 0; b < NBINS; ++b) {
        if (b_count[b] > 0) {
          const double cnt = static_cast<double>(b_count[b]);
          const double mu = b_sum[b] / cnt;
          const double lo = b_min[b];
          const double hi = b_max[b];
          double w_lo = 1.0;
          double w_hi = 0.0;
          if (hi > lo) {
            w_hi = (mu - lo) / (hi - lo);
            w_lo = (hi - mu) / (hi - lo);
          }
          cached_bins.push_back({cnt, mu, lo, hi, w_lo, w_hi});
        }
      }
    }

    history.clear();
    history.reserve(512);
    cbar_cache.clear();
    cbar_cache.reserve(64);
  }

  const double* vd = cached_vd_raw.data();
  const double inv_m = 1.0 / static_cast<double>(m);

  constexpr int CHUNK = 2048;
  double arg_buf[CHUNK];
  double exp_buf[CHUNK];

  auto eval_convex_bounds = [&](double c, double& lb_out, double& ub_out) {
    const int nb = static_cast<int>(cached_bins.size());
    const double neg_c = -c;
    double lb_sum = 0.0;
    double ub_sum = 0.0;
    for (int idx = 0; idx < nb; ++idx) {
      const auto& bin = cached_bins[idx];
      const double e_mu = std::exp(neg_c * bin.mean);
      lb_sum += bin.count * e_mu;
      if (bin.hi > bin.lo) {
        const double e_lo = std::exp(neg_c * bin.lo);
        const double e_hi = std::exp(neg_c * bin.hi);
        ub_sum += bin.count * (bin.w_lo * e_lo + bin.w_hi * e_hi);
      } else {
        ub_sum += bin.count * e_mu;
      }
    }
    lb_out = lb_sum * inv_m;
    ub_out = ub_sum * inv_m;
  };

  auto eval_mean_exp = [&](double c) -> double {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end() && it->c == c) {
      return it->v;
    }
    double sum = 0.0;
    const double neg_c = -c;
#ifdef __APPLE__
    size_t k = 0;
    while (k < m) {
      const int b = static_cast<int>(std::min<size_t>(CHUNK, m - k));
      for (int i = 0; i < b; ++i) {
        arg_buf[i] = neg_c * vd[k + i];
      }
      vvexp(exp_buf, arg_buf, &b);
      for (int i = 0; i < b; ++i) {
        sum += exp_buf[i];
      }
      k += b;
    }
#else
    for (size_t k = 0; k < m; ++k) {
      sum += std::exp(neg_c * vd[k]);
    }
#endif
    double val = sum * inv_m;
    history.insert(it, {c, val});
    return val;
  };

  auto is_below_rhobar = [&](double c, double rhobar) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) return it->v < rhobar;
      if (it->v >= rhobar) return false;
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v < rhobar) return true;
    }
    if (!cached_bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (ub + 1e-12 < rhobar) return true;
      if (lb - 1e-12 > rhobar) return false;
    }
    return eval_mean_exp(c) < rhobar;
  };

  auto is_above_rhobar = [&](double c, double rhobar) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) return it->v > rhobar;
      if (it->v > rhobar) return true;
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v <= rhobar) return false;
    }
    if (!cached_bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (lb - 1e-12 > rhobar) return true;
      if (ub + 1e-12 < rhobar) return false;
    }
    return eval_mean_exp(c) > rhobar;
  };

  NumericVector out(n_rho);
  for (int r = 0; r < n_rho; ++r) {
    const double rhobar = rhobar_vec[r];
    uint64_t rho_bits;
    std::memcpy(&rho_bits, &rhobar, sizeof(double));
    auto c_it = cbar_cache.find(rho_bits);
    if (c_it != cbar_cache.end()) {
      out[r] = c_it->second;
      continue;
    }
    double c0 = 10.0;
    double c1 = 10.0;

    bool i1 = false;
    int jj = 0;
    while (!i1) {
      i1 = is_above_rhobar(c0, rhobar);
      if (!i1) {
        c1 = c0;
        c0 = c0 * 0.5;
        jj++;
      }
      if (jj > 500) Rcpp::stop("rhobar too large");
    }

    i1 = false;
    jj = 0;
    while (!i1) {
      i1 = is_below_rhobar(c1, rhobar);
      if (!i1) {
        c0 = c1;
        c1 = 2.0 * c1;
        jj++;
      }
      if (c1 > 10000.0) i1 = true;
      if (jj > 500) Rcpp::stop("rhobar too small");
    }

    while ((c1 - c0) > 0.001) {
      double cm = std::sqrt(c0 * c1);
      if (is_below_rhobar(cm, rhobar)) {
        c1 = cm;
      } else {
        c0 = cm;
      }
    }

    const double cbar_val = std::sqrt(c0 * c1);
    cbar_cache.emplace(rho_bits, cbar_val);
    out[r] = cbar_val;
  }

  return out;
}

// [[Rcpp::export]]
NumericVector rcpp_om_grid(NumericMatrix distmat, NumericMatrix r_proj, NumericVector c_vec, int nthreads = -1) {
  const int n = distmat.nrow();
  const int q = r_proj.ncol();
  const int n_c = c_vec.size();
  const double* d_ptr_mat = distmat.begin();
  const double* r_ptr = r_proj.begin();
  const double* c_ptr = c_vec.begin();

  NumericVector out(static_cast<size_t>(q) * q * n_c);
  if (n_c == 1) {
    out.attr("dim") = Dimension(q, q);
  } else {
    out.attr("dim") = Dimension(q, q, n_c);
  }
  double* out_ptr = out.begin();

  constexpr int REFRESH_STRIDE = 8;
  std::vector<uint8_t> step_mode(n_c, 0);
  std::vector<double> step_dc(n_c, 0.0);
  {
    double active_dc = -1.0;
    int step_in_run = 0;
    for (int m = 0; m < n_c; ++m) {
      if (m == 0) {
        step_mode[0] = 0;
        continue;
      }
      const double dc = c_ptr[m] - c_ptr[m - 1];
      if (dc > 0.0 && active_dc > 0.0 && std::fabs(dc - active_dc) <= 1e-12) {
        if ((step_in_run % REFRESH_STRIDE) == 0) {
          step_mode[m] = 0;
          // At refresh boundary, next step must re-ensure E_blk if run continues in a new parallel chunk
          if (m + 1 < n_c && std::fabs((c_ptr[m + 1] - c_ptr[m]) - active_dc) <= 1e-12) {
            step_mode[m + 1] = 1;
            step_dc[m + 1] = active_dc;
          }
        } else if (step_mode[m] != 1) {
          step_mode[m] = 2;
        }
        step_in_run++;
      } else if (dc > 0.0 && m + 1 < n_c && std::fabs((c_ptr[m + 1] - c_ptr[m]) - dc) <= 1e-12) {
        active_dc = dc;
        step_mode[m] = 1;
        step_dc[m] = active_dc;
        step_in_run = 2;
      } else {
        active_dc = -1.0;
        step_mode[m] = 0;
        step_in_run = 0;
      }
    }
  }

  // Partition [0, n_c) into independent chunks starting at mode 0 boundaries
  std::vector<std::pair<int, int>> chunks;
  {
    int start = 0;
    for (int m = 1; m < n_c; ++m) {
      if (step_mode[m] == 0) {
        chunks.emplace_back(start, m);
        start = m;
      }
    }
    if (start < n_c) {
      chunks.emplace_back(start, n_c);
    }
  }

  OmGridChunkWorker worker(d_ptr_mat, r_ptr, c_ptr, step_mode.data(), step_dc.data(),
                           chunks, out_ptr, n, q);
  if (chunks.size() > 1) {
    spur_parallel_for(0, chunks.size(), worker, 1, nthreads);
  } else {
    worker(0, chunks.size());
  }

  return out;
}

// [[Rcpp::export]]
double rcpp_pow_qf_eval(NumericMatrix ho, NumericMatrix ha, NumericMatrix e, NumericVector qe, int nthreads = -1) {
  const int q = e.nrow();
  const int nrep = e.ncol();

  std::vector<double> lr_o(nrep);
  std::vector<double> lr_a(nrep);

  constexpr int TILE = 512;
  const int n_tiles = (nrep + TILE - 1) / TILE;

  PowQfTileWorker worker(ho.begin(), ha.begin(), e.begin(), qe.begin(),
                         lr_o.data(), lr_a.data(), q, nrep, TILE);
  if (n_tiles >= 8) {
    spur_parallel_for(0, n_tiles, worker, 8, nthreads);
  } else {
    worker(0, n_tiles);
  }

  const double prob95 = 0.95;
  double cv = 0.0;
  compute_quantiles_type7(lr_o, &prob95, 1, &cv);

  int count_gt = 0;
  for (int k = 0; k < nrep; ++k) {
    if (lr_a[k] > cv) count_gt++;
  }

  return static_cast<double>(count_gt) / static_cast<double>(nrep);
}

// [[Rcpp::export]]
List rcpp_mc_lr_grid(NumericVector ch_om_ho_all,
                     NumericMatrix ch_omi_ho,
                     NumericMatrix ch_omi_ha,
                     NumericMatrix emat,
                     NumericVector lr_stat,
                     NumericVector probs,
                     int nthreads = -1) {
  const int q = emat.nrow();
  const int nrep = emat.ncol();
  const int n_y = lr_stat.size();
  const int n_probs = probs.size();
  const int n_rho = ch_om_ho_all.size() / (q * q);

  NumericMatrix cvalue_mat(Dimension(n_rho, n_probs));
  NumericMatrix pvalue_mat(Dimension(n_rho, n_y));

  McLrGridWorker worker(ch_om_ho_all.begin(), ch_omi_ho.begin(), ch_omi_ha.begin(),
                        emat.begin(), lr_stat.begin(), probs.begin(),
                        cvalue_mat.begin(), pvalue_mat.begin(),
                        q, nrep, n_y, n_probs, n_rho);
  if (n_rho > 1) {
    spur_parallel_for(0, n_rho, worker, 1, nthreads);
  } else {
    worker(0, n_rho);
  }

  return List::create(cvalue_mat, pvalue_mat);
}

// [[Rcpp::export]]
NumericMatrix rcpp_c_ci_mc(NumericVector ch_om_ho_all,
                           NumericVector ch_omi_ho_all,
                           NumericVector const_den_ho,
                           NumericVector ch_omi_ha_all,
                           NumericVector omi_ha_all,
                           NumericVector const_den_ha,
                           NumericMatrix emat,
                           NumericVector X,
                           int nthreads = -1) {
  const int q = emat.nrow();
  const int nrep = emat.ncol();
  const int n_c_ho = const_den_ho.size();
  const int n_c_ha = const_den_ha.size();

  const double* ch_om_ho_ptr = ch_om_ho_all.begin();
  const double* ch_omi_ho_ptr = ch_omi_ho_all.begin();
  const double* const_den_ho_ptr = const_den_ho.begin();
  const double* ch_omi_ha_ptr = ch_omi_ha_all.begin();
  const double* omi_ha_ptr = omi_ha_all.begin();
  const double* const_den_ha_ptr = const_den_ha.begin();
  const double* emat_ptr = emat.begin();
  const double* X_ptr = X.begin();

  double sum_den_ha_X = 0.0;
  for (int j = 0; j < n_c_ha; ++j) {
    const double* ch_omi_ha = ch_omi_ha_ptr + static_cast<size_t>(j) * q * q;
    double q_X = 0.0;
    for (int r = 0; r < q; ++r) {
      double acc = 0.0;
      for (int c = 0; c < q; ++c) {
        acc += ch_omi_ha[r + static_cast<size_t>(c) * q] * X_ptr[c];
      }
      q_X += acc * acc;
    }
    sum_den_ha_X += const_den_ha_ptr[j] * fast_neg_half_pow(q_X, q);
  }
  const double den_ha_avg_X = sum_den_ha_X / static_cast<double>(n_c_ha);

  const int p = q * (q + 1) / 2;
  std::vector<int> idx_a(p);
  std::vector<int> idx_b(p);
  {
    int k = 0;
    for (int b = 0; b < q; ++b) {
      for (int a = 0; a <= b; ++a) {
        idx_a[k] = a;
        idx_b[k] = b;
        k++;
      }
    }
  }

  std::vector<double> W_ha(static_cast<size_t>(n_c_ha) * p);
  for (int j = 0; j < n_c_ha; ++j) {
    const double* omi_ha = omi_ha_ptr + static_cast<size_t>(j) * q * q;
    for (int k = 0; k < p; ++k) {
      const int a = idx_a[k];
      const int b = idx_b[k];
      const double w = (a == b) ? omi_ha[a + static_cast<size_t>(b) * q]
                                : 2.0 * omi_ha[a + static_cast<size_t>(b) * q];
      W_ha[j + static_cast<size_t>(k) * n_c_ha] = w;
    }
  }

  NumericMatrix pv_mat(Dimension(n_c_ho, 1));
  double* pv_ptr = pv_mat.begin();

  std::vector<double> emat_pow(nrep);
  for (int col = 0; col < nrep; ++col) {
    const double* e_col = emat_ptr + static_cast<size_t>(col) * q;
    double q_ho = 0.0;
    for (int r = 0; r < q; ++r) {
      q_ho += e_col[r] * e_col[r];
    }
    emat_pow[col] = fast_neg_half_pow(q_ho, q);
  }

  CCiMcWorker worker(ch_om_ho_ptr, ch_omi_ho_ptr, const_den_ho_ptr, const_den_ha_ptr,
                     emat_ptr, emat_pow.data(), W_ha.data(), X_ptr, pv_ptr,
                     den_ha_avg_X, q, p, nrep, n_c_ha);
  spur_parallel_for(0, n_c_ho, worker, 1, nthreads);

  return pv_mat;
}

// [[Rcpp::export]]
NumericMatrix rcpp_lbm_gls_from_sigma(NumericMatrix sigma) {
  const int n = sigma.nrow();
  const arma::mat sigma_mat(sigma.begin(), n, n, false, true);

  arma::vec vals;
  arma::mat vecs;
  if (!arma::eig_sym(vals, vecs, sigma_mat, "dc")) {
    Rcpp::stop("arma::eig_sym (Divide-and-Conquer dsyevd) failed.");
  }

  int first_keep = n;
  for (int j = 0; j < n; ++j) {
    if (vals[j] > 1e-10) {
      first_keep = j;
      break;
    }
  }
  const int k_keep = n - first_keep;

  NumericMatrix out(Dimension(n, n));
  if (k_keep <= 0) {
    return out;
  }

  arma::mat V_keep = vecs.cols(first_keep, n - 1);
  for (int j = 0; j < k_keep; ++j) {
    const double scale = 1.0 / std::sqrt(std::sqrt(vals[first_keep + j]));
    V_keep.col(j) *= scale;
  }

  arma::mat out_mat(out.begin(), n, n, false, true);
  out_mat = V_keep * V_keep.t();
  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_get_sigma_residual_proj(NumericMatrix distmat,
                                           double c,
                                           NumericMatrix x_in,
                                           NumericMatrix p_x) {
  const int n = distmat.nrow();
  const int k = x_in.ncol();

  NumericMatrix s_mat(Dimension(n, n));
  double* S_ptr = s_mat.begin();

  exp_lower_triangle_inplace(distmat.begin(), S_ptr, n, -c);
  symmetrize_lower_to_upper_blocked(S_ptr, n);

  arma::mat S(S_ptr, n, n, false, true);
  const arma::mat X_in(x_in.begin(), n, k, false, true);
  const arma::mat P_X(p_x.begin(), n, k, false, true);

  S -= P_X * (X_in.t() * S);
  S -= (S * X_in) * P_X.t();

  return s_mat;
}

// [[Rcpp::export]]
NumericMatrix rcpp_nn_matrix(NumericMatrix distmat) {
  const int n = distmat.nrow();
  const double* d_ptr = distmat.begin();

  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  for (int i = 0; i < n; ++i) {
    const double* col_i = d_ptr + static_cast<size_t>(i) * n;
    double min_val = 1e10;
    for (int j = 0; j < n; ++j) {
      const double d = (col_i[j] == 0.0) ? 1e10 : col_i[j];
      if (d < min_val) min_val = d;
    }
    int count = 0;
    for (int j = 0; j < n; ++j) {
      const double d = (col_i[j] == 0.0) ? 1e10 : col_i[j];
      if (d == min_val) count++;
    }
    const double neg_w = -1.0 / static_cast<double>(count);
    for (int j = 0; j < n; ++j) {
      const double d = (col_i[j] == 0.0) ? 1e10 : col_i[j];
      if (d == min_val) {
        out_ptr[i + static_cast<size_t>(j) * n] = neg_w;
      }
    }
    out_ptr[i + static_cast<size_t>(i) * n] += 1.0;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_iso_matrix(NumericMatrix distmat, double radius) {
  const int n = distmat.nrow();
  const double* d_ptr = distmat.begin();

  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  for (int i = 0; i < n; ++i) {
    const double* col_i = d_ptr + static_cast<size_t>(i) * n;
    int count = 0;
    for (int j = 0; j < n; ++j) {
      const double d = (col_i[j] == 0.0) ? 1e10 : col_i[j];
      if (d <= radius) count++;
    }
    if (count > 0) {
      const double neg_w = -1.0 / static_cast<double>(count);
      for (int j = 0; j < n; ++j) {
        const double d = (col_i[j] == 0.0) ? 1e10 : col_i[j];
        if (d <= radius) {
          out_ptr[i + static_cast<size_t>(j) * n] = neg_w;
        }
      }
      out_ptr[i + static_cast<size_t>(i) * n] += 1.0;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix rcpp_cluster_matrix(NumericVector cluster) {
  const int n = cluster.size();
  const double* c_ptr = cluster.begin();

  std::vector<int> counts(n, 0);
  for (int i = 0; i < n; ++i) {
    const double ci = c_ptr[i];
    int cnt = 0;
    for (int j = 0; j < n; ++j) {
      if (c_ptr[j] == ci) cnt++;
    }
    counts[i] = cnt;
  }

  NumericMatrix out(Dimension(n, n));
  double* out_ptr = out.begin();

  for (int j = 0; j < n; ++j) {
    const double cj = c_ptr[j];
    double* col = out_ptr + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      double val = (i == j) ? 1.0 : 0.0;
      if (c_ptr[i] == cj) {
        val -= 1.0 / static_cast<double>(counts[i]);
      }
      col[i] = val;
    }
  }

  return out;
}
