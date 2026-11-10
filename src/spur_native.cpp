#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <vector>

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#ifndef FCONE
#define FCONE
#endif

#ifdef __APPLE__
extern "C" void vvexp(double* y, const double* x, const int* n);
#endif

namespace {

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

  // Subtract column means and accumulate row sums of demeaned columns
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

  // Subtract row means
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

}  // namespace

extern "C" {

SEXP cpp_distmat_lat_lon(SEXP s_sexp) {
  const int n = Rf_nrows(s_sexp);
  const double* s = REAL(s_sexp);
  const double* lat = s;
  const double* lon = s + n;
  const double c = M_PI / 180.0;
  const double inv_pi = 1.0 / M_PI;

  std::vector<double> cos_lat(n);
  for (int i = 0; i < n; ++i) {
    cos_lat[i] = std::cos(c * lat[i]);
  }

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  for (int i = 0; i < n; ++i) {
    out[i + static_cast<size_t>(i) * n] = 0.0;
    const double lat_i = lat[i];
    const double lon_i = lon[i];
    const double cos_i = cos_lat[i];
    double* out_col = out + static_cast<size_t>(i) * n;
    for (int j = i + 1; j < n; ++j) {
      const double d1 = 0.5 * c * (lat[j] - lat_i);
      const double d2 = 0.5 * c * (lon[j] - lon_i);
      const double sd1 = std::sin(d1);
      const double sd2 = std::sin(d2);
      const double inside = sd1 * sd1 + cos_i * (cos_lat[j] * (sd2 * sd2));
      out_col[j] = std::asin(std::sqrt(inside)) * inv_pi;
    }
  }
  symmetrize_lower_to_upper_blocked(out, n);

  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_distmat_euclidean(SEXP s_sexp) {
  const int n = Rf_nrows(s_sexp);
  const int d = Rf_ncols(s_sexp);
  const double* s = REAL(s_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  if (d == 2) {
    const double* s0 = s;
    const double* s1 = s + n;
    for (int i = 0; i < n; ++i) {
      double* out_col = out + static_cast<size_t>(i) * n;
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
      double* out_col = out + static_cast<size_t>(i) * n;
      out_col[i] = 0.0;
      for (int j = i + 1; j < n; ++j) {
        double sum_sq = 0.0;
        for (int k = 0; k < d; ++k) {
          const double diff = s[j + static_cast<size_t>(k) * n] - s[i + static_cast<size_t>(k) * n];
          sum_sq += diff * diff;
        }
        out_col[j] = std::sqrt(sum_sq);
      }
    }
  }
  symmetrize_lower_to_upper_blocked(out, n);

  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_demean_sigma(SEXP sigma_sexp) {
  const int n = Rf_nrows(sigma_sexp);
  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  std::memcpy(REAL(out_sexp), REAL(sigma_sexp), static_cast<size_t>(n) * n * sizeof(double));
  demean_matrix_inplace(REAL(out_sexp), n);
  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_get_sigma_dm(SEXP distmat_sexp, SEXP c_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double c = Rf_asReal(c_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  exp_lower_triangle_inplace(distmat, out, n, -c);
  symmetrize_lower_to_upper_blocked(out, n);

  demean_matrix_inplace(out, n);
  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_get_sigma_lbm(SEXP distmat_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  for (int j = 0; j < n; ++j) {
    const double d_0j = distmat[static_cast<size_t>(j) * n];
    const double* d_col = distmat + static_cast<size_t>(j) * n;
    double* out_col = out + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      out_col[i] = 0.5 * (distmat[i] + d_0j - d_col[i]);
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_get_sigma_lbm_dm(SEXP distmat_sexp) {
  SEXP out_sexp = PROTECT(cpp_get_sigma_lbm(distmat_sexp));
  demean_matrix_inplace(REAL(out_sexp), Rf_nrows(distmat_sexp));
  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_get_cbar_vec(SEXP rhobar_sexp, SEXP distmat_sexp) {
  const int n_rho = Rf_length(rhobar_sexp);
  const double* rhobar_vec = REAL(rhobar_sexp);
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);

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

  // Quick signature check + exact lower-triangle verification
  bool match = (cached_n == n && cached_vd_raw.size() == m);
  if (match && m > 0) {
    size_t pos = 0;
    for (int j = 0; j < n - 1 && match; ++j) {
      const double* col = distmat + static_cast<size_t>(j) * n + (j + 1);
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
      const double* col = distmat + static_cast<size_t>(j) * n + (j + 1);
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

  // Exact Jensen lower bound and Edmundson-Madansky upper bound from 1024 histogram bins
  auto eval_convex_bounds = [&](double c, double& lb_out, double& ub_out) {
    const int nb = static_cast<int>(cached_bins.size());
    const double neg_c = -c;
    double lb_sum = 0.0;
    double ub_sum = 0.0;
#ifdef __APPLE__
    // 3 values per bin (mean, lo, hi) <= 3072 total
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
#else
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
#endif
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

  // Returns true iff f(c) < rhobar, using history bounds + 1024-bin convexity bounds before full scan
  auto is_below_rhobar = [&](double c, double rhobar) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) {
        return it->v < rhobar;
      }
      if (it->v >= rhobar) {
        return false;
      }
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v < rhobar) {
        return true;
      }
    }

    if (!cached_bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (ub + 1e-12 < rhobar) {
        return true;
      }
      if (lb - 1e-12 > rhobar) {
        return false;
      }
    }

    return eval_mean_exp(c) < rhobar;
  };

  auto is_above_rhobar = [&](double c, double rhobar) -> bool {
    auto it = std::lower_bound(
        history.begin(), history.end(), c,
        [](const EvalPoint& pt, double val) { return pt.c < val; });
    if (it != history.end()) {
      if (it->c == c) {
        return it->v > rhobar;
      }
      if (it->v > rhobar) {
        return true;
      }
    }
    if (it != history.begin()) {
      auto prev = it - 1;
      if (prev->v <= rhobar) {
        return false;
      }
    }

    if (!cached_bins.empty()) {
      double lb = 0.0, ub = 1.0;
      eval_convex_bounds(c, lb, ub);
      if (lb - 1e-12 > rhobar) {
        return true;
      }
      if (ub + 1e-12 < rhobar) {
        return false;
      }
    }

    return eval_mean_exp(c) > rhobar;
  };

  SEXP out_sexp = PROTECT(Rf_allocVector(REALSXP, n_rho));
  double* out = REAL(out_sexp);

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
      if (jj > 500) {
        UNPROTECT(1);
        Rf_error("rhobar too large");
      }
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
      if (c1 > 10000.0) {
        i1 = true;
      }
      if (jj > 500) {
        UNPROTECT(1);
        Rf_error("rhobar too small");
      }
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

  UNPROTECT(1);
  return out_sexp;
}

// Compute R_proj^T * exp(-c * distmat) * R_proj for each c in c_vec
SEXP cpp_om_grid(SEXP distmat_sexp, SEXP r_proj_sexp, SEXP c_vec_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const int q = Rf_ncols(r_proj_sexp);
  const int n_c = Rf_length(c_vec_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double* r_proj = REAL(r_proj_sexp);
  const double* c_vec = REAL(c_vec_sexp);

  SEXP out_sexp;
  if (n_c == 1) {
    out_sexp = PROTECT(Rf_allocMatrix(REALSXP, q, q));
  } else {
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(dim)[0] = q;
    INTEGER(dim)[1] = q;
    INTEGER(dim)[2] = n_c;
    out_sexp = PROTECT(Rf_allocArray(REALSXP, dim));
  }
  double* out = REAL(out_sexp);

  constexpr int BLK = 256;
  std::vector<double> SR_all(static_cast<size_t>(n) * q * n_c, 0.0);
  std::vector<double> D_blk(BLK * BLK);
  std::vector<double> S_blk(BLK * BLK);
  std::vector<double> E_blk(BLK * BLK);

  // Precompute grid step pattern across m = 0 .. n_c - 1
  constexpr int REFRESH_STRIDE = 8;
  std::vector<uint8_t> step_mode(n_c, 0);  // 0 = exact exp, 1 = init E_blk + mul, 2 = mul E_blk
  std::vector<double> step_dc(n_c, 0.0);
  {
    double active_dc = -1.0;
    int step_in_run = 0;
    for (int m = 0; m < n_c; ++m) {
      if (m == 0) {
        step_mode[0] = 0;
        continue;
      }
      const double dc = c_vec[m] - c_vec[m - 1];
      if (dc > 0.0 && active_dc > 0.0 && std::fabs(dc - active_dc) <= 1e-12) {
        if ((step_in_run % REFRESH_STRIDE) == 0) {
          step_mode[m] = 0;
        } else {
          step_mode[m] = 2;
        }
        step_in_run++;
      } else if (dc > 0.0 && m + 1 < n_c && std::fabs((c_vec[m + 1] - c_vec[m]) - dc) <= 1e-12) {
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

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta_one = 1.0;
  const double beta_zero = 0.0;

  for (int j0 = 0; j0 < n; j0 += BLK) {
    const int bj = std::min(BLK, n - j0);
    for (int i0 = j0; i0 < n; i0 += BLK) {
      const int bi = std::min(BLK, n - i0);
      const int len = bi * bj;

      // Copy block (i0..i0+bi-1, j0..j0+bj-1) of distmat into contiguous L1/L2 buffer D_blk
      for (int j = 0; j < bj; ++j) {
        const double* d_col = distmat + static_cast<size_t>(j0 + j) * n + i0;
        std::memcpy(D_blk.data() + static_cast<size_t>(j) * bi, d_col, bi * sizeof(double));
      }

      const double* R_i = r_proj + i0;
      const double* R_j = r_proj + j0;

      for (int m = 0; m < n_c; ++m) {
        const uint8_t mode = step_mode[m];
        if (mode == 2) {
          double* s_ptr = S_blk.data();
          const double* e_ptr = E_blk.data();
          for (int k = 0; k < len; ++k) {
            s_ptr[k] *= e_ptr[k];
          }
        } else if (mode == 1) {
          const double neg_dc = -step_dc[m];
          const double* d_ptr = D_blk.data();
          double* e_ptr = E_blk.data();
          double* s_ptr = S_blk.data();
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
          const double neg_c = -c_vec[m];
          const double* d_ptr = D_blk.data();
          double* s_ptr = S_blk.data();
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

        double* SR_m = SR_all.data() + static_cast<size_t>(m) * n * q;
        if (i0 == j0) {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &q, &bj, &alpha,
                          S_blk.data(), &bi, R_i, &n, &beta_one, SR_m + i0, &n FCONE FCONE);
        } else {
          F77_CALL(dgemm)(&trans_n, &trans_n, &bi, &q, &bj, &alpha,
                          S_blk.data(), &bi, R_j, &n, &beta_one, SR_m + i0, &n FCONE FCONE);
          F77_CALL(dgemm)(&trans_t, &trans_n, &bj, &q, &bi, &alpha,
                          S_blk.data(), &bi, R_i, &n, &beta_one, SR_m + j0, &n FCONE FCONE);
        }
      }
    }
  }

  for (int m = 0; m < n_c; ++m) {
    const double* SR_m = SR_all.data() + static_cast<size_t>(m) * n * q;
    double* om_m = out + static_cast<size_t>(m) * q * q;
    F77_CALL(dgemm)(&trans_t, &trans_n, &q, &q, &n, &alpha,
                    r_proj, &n, SR_m, &n, &beta_zero, om_m, &q FCONE FCONE);
  }

  UNPROTECT(n_c == 1 ? 1 : 2);
  return out_sexp;
}

// Evaluate power of quadratic form given ho (q x q), ha (q x q), e (q x nrep), qe (nrep)
SEXP cpp_pow_qf_eval(SEXP ho_sexp, SEXP ha_sexp, SEXP e_sexp, SEXP qe_sexp) {
  const int q = Rf_nrows(e_sexp);
  const int nrep = Rf_ncols(e_sexp);
  const double* ho = REAL(ho_sexp);
  const double* ha = REAL(ha_sexp);
  const double* e = REAL(e_sexp);
  const double* qe = REAL(qe_sexp);

  std::vector<double> lr_o(nrep);
  std::vector<double> lr_a(nrep);

  constexpr int TILE = 512;
  std::vector<double> ya_tile(static_cast<size_t>(q) * TILE);
  std::vector<double> yo_tile(static_cast<size_t>(q) * TILE);

  const char trans_n = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;

  for (int col0 = 0; col0 < nrep; col0 += TILE) {
    const int b = std::min(TILE, nrep - col0);
    const double* e_tile = e + static_cast<size_t>(col0) * q;

    F77_CALL(dgemm)(&trans_n, &trans_n, &q, &b, &q, &alpha,
                    ho, &q, e_tile, &q, &beta, ya_tile.data(), &q FCONE FCONE);
    F77_CALL(dgemm)(&trans_n, &trans_n, &q, &b, &q, &alpha,
                    ha, &q, e_tile, &q, &beta, yo_tile.data(), &q FCONE FCONE);

    for (int j = 0; j < b; ++j) {
      const double* ya_col = ya_tile.data() + static_cast<size_t>(j) * q;
      const double* yo_col = yo_tile.data() + static_cast<size_t>(j) * q;
      double qa_o = 0.0;
      double qo_a = 0.0;
      for (int i = 0; i < q; ++i) {
        qa_o += ya_col[i] * ya_col[i];
        qo_a += yo_col[i] * yo_col[i];
      }
      const double qe_k = qe[col0 + j];
      lr_o[col0 + j] = qe_k / qa_o;
      lr_a[col0 + j] = qo_a / qe_k;
    }
  }

  const double prob95 = 0.95;
  double cv = 0.0;
  compute_quantiles_type7(lr_o, &prob95, 1, &cv);

  int count_gt = 0;
  for (int k = 0; k < nrep; ++k) {
    if (lr_a[k] > cv) count_gt++;
  }

  return Rf_ScalarReal(static_cast<double>(count_gt) / static_cast<double>(nrep));
}

// Evaluate Monte Carlo LR grid for spur_i0_test / spur_i0_resid_test / spur_i1_test / spur_i1_resid_test
SEXP cpp_mc_lr_grid(SEXP ch_om_ho_sexp,
                    SEXP ch_omi_ho_sexp,
                    SEXP ch_omi_ha_sexp,
                    SEXP emat_sexp,
                    SEXP lr_stat_sexp,
                    SEXP probs_sexp) {
  const int q = Rf_nrows(emat_sexp);
  const int nrep = Rf_ncols(emat_sexp);
  const int n_y = Rf_length(lr_stat_sexp);
  const int n_probs = Rf_length(probs_sexp);
  const int total_ch = Rf_length(ch_om_ho_sexp);
  const int n_rho = total_ch / (q * q);

  const double* ch_om_ho_all = REAL(ch_om_ho_sexp);
  const double* ch_omi_ho = REAL(ch_omi_ho_sexp);
  const double* ch_omi_ha = REAL(ch_omi_ha_sexp);
  const double* emat = REAL(emat_sexp);
  const double* lr_stat = REAL(lr_stat_sexp);
  const double* probs = REAL(probs_sexp);

  SEXP cvalue_mat_sexp = PROTECT(Rf_allocMatrix(REALSXP, n_rho, n_probs));
  SEXP pvalue_mat_sexp = PROTECT(Rf_allocMatrix(REALSXP, n_rho, n_y));
  double* cvalue_mat = REAL(cvalue_mat_sexp);
  double* pvalue_mat = REAL(pvalue_mat_sexp);

  constexpr int TILE = 512;
  std::vector<double> A_ho(static_cast<size_t>(q) * q);
  std::vector<double> A_ha(static_cast<size_t>(q) * q);
  std::vector<double> y_ho_ho_tile(static_cast<size_t>(q) * TILE);
  std::vector<double> y_ho_ha_tile(static_cast<size_t>(q) * TILE);
  std::vector<double> lr_ho(nrep);
  std::vector<double> cv_row(n_probs);
  std::vector<int> gt_counts(n_y);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;

  for (int ir = 0; ir < n_rho; ++ir) {
    const double* ch_om_ho = ch_om_ho_all + static_cast<size_t>(ir) * q * q;
    std::fill(gt_counts.begin(), gt_counts.end(), 0);

    // Pre-multiply A_ho = ch_omi_ho %*% t(ch_om_ho) and A_ha = ch_omi_ha %*% t(ch_om_ho)
    F77_CALL(dgemm)(&trans_n, &trans_t, &q, &q, &q, &alpha,
                    ch_omi_ho, &q, ch_om_ho, &q, &beta, A_ho.data(), &q FCONE FCONE);
    F77_CALL(dgemm)(&trans_n, &trans_t, &q, &q, &q, &alpha,
                    ch_omi_ha, &q, ch_om_ho, &q, &beta, A_ha.data(), &q FCONE FCONE);

    for (int col0 = 0; col0 < nrep; col0 += TILE) {
      const int b = std::min(TILE, nrep - col0);
      const double* emat_tile = emat + static_cast<size_t>(col0) * q;

      F77_CALL(dgemm)(&trans_n, &trans_n, &q, &b, &q, &alpha,
                      A_ho.data(), &q, emat_tile, &q, &beta, y_ho_ho_tile.data(), &q FCONE FCONE);
      F77_CALL(dgemm)(&trans_n, &trans_n, &q, &b, &q, &alpha,
                      A_ha.data(), &q, emat_tile, &q, &beta, y_ho_ha_tile.data(), &q FCONE FCONE);

      for (int j = 0; j < b; ++j) {
        const double* col_ho = y_ho_ho_tile.data() + static_cast<size_t>(j) * q;
        const double* col_ha = y_ho_ha_tile.data() + static_cast<size_t>(j) * q;
        double q_ho = 0.0;
        double q_ha = 0.0;
        for (int i = 0; i < q; ++i) {
          q_ho += col_ho[i] * col_ho[i];
          q_ha += col_ha[i] * col_ha[i];
        }
        const double lr_val = q_ho / q_ha;
        lr_ho[col0 + j] = lr_val;
        for (int jy = 0; jy < n_y; ++jy) {
          if (lr_val > lr_stat[jy]) gt_counts[jy]++;
        }
      }
    }

    compute_quantiles_type7(lr_ho, probs, n_probs, cv_row.data());
    for (int jp = 0; jp < n_probs; ++jp) {
      cvalue_mat[ir + static_cast<size_t>(jp) * n_rho] = cv_row[jp];
    }
    for (int jy = 0; jy < n_y; ++jy) {
      pvalue_mat[ir + static_cast<size_t>(jy) * n_rho] =
          static_cast<double>(gt_counts[jy]) / static_cast<double>(nrep);
    }
  }

  SEXP out_sexp = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_sexp, 0, cvalue_mat_sexp);
  SET_VECTOR_ELT(out_sexp, 1, pvalue_mat_sexp);
  UNPROTECT(3);
  return out_sexp;
}

// Fast Monte Carlo evaluation for c_ci (spur_persistence)
SEXP cpp_c_ci_mc(SEXP ch_om_ho_sexp,
                 SEXP ch_omi_ho_sexp,
                 SEXP const_den_ho_sexp,
                 SEXP ch_omi_ha_sexp,
                 SEXP omi_ha_sexp,
                 SEXP const_den_ha_sexp,
                 SEXP emat_sexp,
                 SEXP X_sexp) {
  const int q = Rf_nrows(emat_sexp);
  const int nrep = Rf_ncols(emat_sexp);
  const int n_c_ho = Rf_length(const_den_ho_sexp);
  const int n_c_ha = Rf_length(const_den_ha_sexp);

  const double* ch_om_ho_all = REAL(ch_om_ho_sexp);
  const double* ch_omi_ho_all = REAL(ch_omi_ho_sexp);
  const double* const_den_ho = REAL(const_den_ho_sexp);
  const double* ch_omi_ha_all = REAL(ch_omi_ha_sexp);
  const double* omi_ha_all = REAL(omi_ha_sexp);
  const double* const_den_ha = REAL(const_den_ha_sexp);
  const double* emat = REAL(emat_sexp);
  const double* X = REAL(X_sexp);

  // 1) Precompute den_ha_avg_X across j = 0 .. n_c_ha - 1
  double sum_den_ha_X = 0.0;
  for (int j = 0; j < n_c_ha; ++j) {
    const double* ch_omi_ha = ch_omi_ha_all + static_cast<size_t>(j) * q * q;
    double q_X = 0.0;
    for (int r = 0; r < q; ++r) {
      double acc = 0.0;
      for (int c = 0; c < q; ++c) {
        acc += ch_omi_ha[r + static_cast<size_t>(c) * q] * X[c];
      }
      q_X += acc * acc;
    }
    sum_den_ha_X += const_den_ha[j] * fast_neg_half_pow(q_X, q);
  }
  const double den_ha_avg_X = sum_den_ha_X / static_cast<double>(n_c_ha);

  // 2) Build W_ha (n_c_ha x p) where p = q*(q+1)/2 for vech quadratic form Q = W_ha * U
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
    const double* omi_ha = omi_ha_all + static_cast<size_t>(j) * q * q;
    for (int k = 0; k < p; ++k) {
      const int a = idx_a[k];
      const int b = idx_b[k];
      const double w = (a == b) ? omi_ha[a + static_cast<size_t>(b) * q]
                                : 2.0 * omi_ha[a + static_cast<size_t>(b) * q];
      W_ha[j + static_cast<size_t>(k) * n_c_ha] = w;
    }
  }

  SEXP pv_mat_sexp = PROTECT(Rf_allocMatrix(REALSXP, n_c_ho, 1));
  double* pv_mat = REAL(pv_mat_sexp);

  // Precompute (colSums(emat^2))^(-q/2), since ||ch_omi_i * t(ch_om_i) * emat_k||^2 == ||emat_k||^2
  std::vector<double> emat_pow(nrep);
  for (int col = 0; col < nrep; ++col) {
    const double* e_col = emat + static_cast<size_t>(col) * q;
    double q_ho = 0.0;
    for (int r = 0; r < q; ++r) {
      q_ho += e_col[r] * e_col[r];
    }
    emat_pow[col] = fast_neg_half_pow(q_ho, q);
  }

  constexpr int TILE = 512;
  std::vector<double> e_tile(static_cast<size_t>(q) * TILE);
  std::vector<double> U_tile(static_cast<size_t>(p) * TILE);
  std::vector<double> Q_tile(static_cast<size_t>(n_c_ha) * TILE);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;
  const double inv_n_ha = 1.0 / static_cast<double>(n_c_ha);

  for (int i = 0; i < n_c_ho; ++i) {
    const double* ch_null = ch_om_ho_all + static_cast<size_t>(i) * q * q;
    const double* ch_omi = ch_omi_ho_all + static_cast<size_t>(i) * q * q;
    const double const_den = const_den_ho[i];

    // den_ho_X
    double q_X_ho = 0.0;
    for (int r = 0; r < q; ++r) {
      double acc = 0.0;
      for (int c = 0; c < q; ++c) {
        acc += ch_omi[r + static_cast<size_t>(c) * q] * X[c];
      }
      q_X_ho += acc * acc;
    }
    const double den_ho_X = const_den * fast_neg_half_pow(q_X_ho, q);
    const double lr_X = den_ha_avg_X / den_ho_X;

    int count_gt = 0;

    for (int col0 = 0; col0 < nrep; col0 += TILE) {
      const int b = std::min(TILE, nrep - col0);
      const double* emat_tile = emat + static_cast<size_t>(col0) * q;

      // e_tile = t(ch_null) %*% emat_tile (q x b)
      F77_CALL(dgemm)(&trans_t, &trans_n, &q, &b, &q, &alpha,
                      ch_null, &q, emat_tile, &q, &beta, e_tile.data(), &q FCONE FCONE);

      // Build U_tile (p x b)
      for (int col = 0; col < b; ++col) {
        const double* e_col = e_tile.data() + static_cast<size_t>(col) * q;
        double* u_col = U_tile.data() + static_cast<size_t>(col) * p;
        int k = 0;
        for (int cb = 0; cb < q; ++cb) {
          const double eb = e_col[cb];
          for (int ca = 0; ca <= cb; ++ca) {
            u_col[k++] = e_col[ca] * eb;
          }
        }
      }

      // Q_tile = W_ha * U_tile (n_c_ha x p) * (p x b) = (n_c_ha x b)
      F77_CALL(dgemm)(&trans_n, &trans_n, &n_c_ha, &b, &p, &alpha,
                      W_ha.data(), &n_c_ha, U_tile.data(), &p, &beta, Q_tile.data(), &n_c_ha FCONE FCONE);

      for (int col = 0; col < b; ++col) {
        const double den_ho = const_den * emat_pow[col0 + col];
        const double* q_ha_col = Q_tile.data() + static_cast<size_t>(col) * n_c_ha;
        double sum_ha = 0.0;
        if (q == 15) {
          for (int j = 0; j < n_c_ha; ++j) {
            sum_ha += const_den_ha[j] * fast_neg_half_pow_15(q_ha_col[j]);
          }
        } else {
          for (int j = 0; j < n_c_ha; ++j) {
            sum_ha += const_den_ha[j] * fast_neg_half_pow(q_ha_col[j], q);
          }
        }
        const double den_ha_avg = sum_ha * inv_n_ha;
        const double lr = den_ha_avg / den_ho;
        if (lr > lr_X) count_gt++;
      }
    }

    pv_mat[i] = static_cast<double>(count_gt) / static_cast<double>(nrep);
  }

  UNPROTECT(1);
  return pv_mat_sexp;
}

// Fast LBM-GLS matrix from demeaned LBM covariance matrix via dsyevd + dsyrk
SEXP cpp_lbm_gls_from_sigma(SEXP sigma_sexp) {
  const int n = Rf_nrows(sigma_sexp);
  std::vector<double> vecs(static_cast<size_t>(n) * n);
  std::vector<double> vals(n);
  std::memcpy(vecs.data(), REAL(sigma_sexp), static_cast<size_t>(n) * n * sizeof(double));

  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  int liwork = -1;
  double work_query = 0.0;
  int iwork_query = 0;

  F77_CALL(dsyevd)(&jobz, &uplo, &n, vecs.data(), &n, vals.data(),
                   &work_query, &lwork, &iwork_query, &liwork, &info FCONE FCONE);
  if (info != 0) {
    Rf_error("LAPACK dsyevd workspace query failed with info = %d", info);
  }
  lwork = static_cast<int>(work_query);
  liwork = iwork_query;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);

  F77_CALL(dsyevd)(&jobz, &uplo, &n, vecs.data(), &n, vals.data(),
                   work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE);
  if (info != 0) {
    Rf_error("LAPACK dsyevd failed with info = %d", info);
  }

  // dsyevd returns eigenvalues in ascending order; keep values > 1e-10
  int first_keep = n;
  for (int j = 0; j < n; ++j) {
    if (vals[j] > 1e-10) {
      first_keep = j;
      break;
    }
  }
  const int k_keep = n - first_keep;

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);
  if (k_keep <= 0) {
    std::memset(out, 0, static_cast<size_t>(n) * n * sizeof(double));
    UNPROTECT(1);
    return out_sexp;
  }

  // Scale kept eigenvector columns by val^(-1/4)
  for (int j = first_keep; j < n; ++j) {
    const double scale = 1.0 / std::sqrt(std::sqrt(vals[j]));
    double* col = vecs.data() + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      col[i] *= scale;
    }
  }

  // out = V_keep * V_keep^T via dsyrk
  const char trans = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;
  const double* v_keep = vecs.data() + static_cast<size_t>(first_keep) * n;
  F77_CALL(dsyrk)(&uplo, &trans, &n, &k_keep, &alpha,
                  v_keep, &n, &beta, out, &n FCONE FCONE);

  // Mirror upper triangle to lower triangle
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      out[j + static_cast<size_t>(i) * n] = out[i + static_cast<size_t>(j) * n];
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

// Compute M * exp(-c * distmat) * M^T in O(n^2 k) using P_X = X * (X^T X)^(-1)
SEXP cpp_get_sigma_residual_proj(SEXP distmat_sexp,
                                 SEXP c_sexp,
                                 SEXP x_in_sexp,
                                 SEXP p_x_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const int k = Rf_ncols(x_in_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double c = Rf_asReal(c_sexp);
  const double* x_in = REAL(x_in_sexp);
  const double* p_x = REAL(p_x_sexp);

  SEXP s_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* S = REAL(s_sexp);

  exp_lower_triangle_inplace(distmat, S, n, -c);
  symmetrize_lower_to_upper_blocked(S, n);

  std::vector<double> XS(static_cast<size_t>(k) * n);
  std::vector<double> SX(static_cast<size_t>(n) * k);

  const char trans_n = 'N';
  const char trans_t = 'T';
  const double one = 1.0;
  const double neg_one = -1.0;
  const double zero = 0.0;

  // XS = X_in^T * S (k x n)
  F77_CALL(dgemm)(&trans_t, &trans_n, &k, &n, &n, &one,
                  x_in, &n, S, &n, &zero, XS.data(), &k FCONE FCONE);
  // S = S - P_X * XS (n x n)
  F77_CALL(dgemm)(&trans_n, &trans_n, &n, &n, &k, &neg_one,
                  p_x, &n, XS.data(), &k, &one, S, &n FCONE FCONE);
  // SX = S * X_in (n x k)
  F77_CALL(dgemm)(&trans_n, &trans_n, &n, &k, &n, &one,
                  S, &n, x_in, &n, &zero, SX.data(), &n FCONE FCONE);
  // S = S - SX * P_X^T (n x n)
  F77_CALL(dgemm)(&trans_n, &trans_t, &n, &n, &k, &neg_one,
                  SX.data(), &n, p_x, &n, &one, S, &n FCONE FCONE);

  UNPROTECT(1);
  return s_sexp;
}

SEXP cpp_nn_matrix(SEXP distmat_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);
  std::memset(out, 0, static_cast<size_t>(n) * n * sizeof(double));

  // distmat is symmetric, so column i equals row i
  for (int i = 0; i < n; ++i) {
    const double* col_i = distmat + static_cast<size_t>(i) * n;
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
        out[i + static_cast<size_t>(j) * n] = neg_w;
      }
    }
    out[i + static_cast<size_t>(i) * n] += 1.0;
  }

  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_iso_matrix(SEXP distmat_sexp, SEXP radius_sexp) {
  const int n = Rf_nrows(distmat_sexp);
  const double* distmat = REAL(distmat_sexp);
  const double radius = Rf_asReal(radius_sexp);

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);
  std::memset(out, 0, static_cast<size_t>(n) * n * sizeof(double));

  for (int i = 0; i < n; ++i) {
    const double* col_i = distmat + static_cast<size_t>(i) * n;
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
          out[i + static_cast<size_t>(j) * n] = neg_w;
        }
      }
      out[i + static_cast<size_t>(i) * n] += 1.0;
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

SEXP cpp_cluster_matrix(SEXP cluster_sexp) {
  const int n = Rf_length(cluster_sexp);
  const double* cluster = REAL(cluster_sexp);

  std::vector<int> counts(n, 0);
  for (int i = 0; i < n; ++i) {
    const double ci = cluster[i];
    int cnt = 0;
    for (int j = 0; j < n; ++j) {
      if (cluster[j] == ci) cnt++;
    }
    counts[i] = cnt;
  }

  SEXP out_sexp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
  double* out = REAL(out_sexp);

  for (int j = 0; j < n; ++j) {
    const double cj = cluster[j];
    double* col = out + static_cast<size_t>(j) * n;
    for (int i = 0; i < n; ++i) {
      double val = (i == j) ? 1.0 : 0.0;
      if (cluster[i] == cj) {
        val -= 1.0 / static_cast<double>(counts[i]);
      }
      col[i] = val;
    }
  }

  UNPROTECT(1);
  return out_sexp;
}

static const R_CallMethodDef CallEntries[] = {
    {"cpp_distmat_lat_lon", (DL_FUNC)&cpp_distmat_lat_lon, 1},
    {"cpp_distmat_euclidean", (DL_FUNC)&cpp_distmat_euclidean, 1},
    {"cpp_demean_sigma", (DL_FUNC)&cpp_demean_sigma, 1},
    {"cpp_get_sigma_dm", (DL_FUNC)&cpp_get_sigma_dm, 2},
    {"cpp_get_sigma_lbm", (DL_FUNC)&cpp_get_sigma_lbm, 1},
    {"cpp_get_sigma_lbm_dm", (DL_FUNC)&cpp_get_sigma_lbm_dm, 1},
    {"cpp_get_cbar_vec", (DL_FUNC)&cpp_get_cbar_vec, 2},
    {"cpp_om_grid", (DL_FUNC)&cpp_om_grid, 3},
    {"cpp_pow_qf_eval", (DL_FUNC)&cpp_pow_qf_eval, 4},
    {"cpp_mc_lr_grid", (DL_FUNC)&cpp_mc_lr_grid, 6},
    {"cpp_c_ci_mc", (DL_FUNC)&cpp_c_ci_mc, 8},
    {"cpp_lbm_gls_from_sigma", (DL_FUNC)&cpp_lbm_gls_from_sigma, 1},
    {"cpp_get_sigma_residual_proj", (DL_FUNC)&cpp_get_sigma_residual_proj, 4},
    {"cpp_nn_matrix", (DL_FUNC)&cpp_nn_matrix, 1},
    {"cpp_iso_matrix", (DL_FUNC)&cpp_iso_matrix, 2},
    {"cpp_cluster_matrix", (DL_FUNC)&cpp_cluster_matrix, 1},
    {NULL, NULL, 0}};

void R_init_spuR(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

}  // extern "C"
