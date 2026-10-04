#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cstdint>
using namespace Rcpp;

// The autoregressive structure and its derivatives, one exported kernel per
// order. The free values are eta_0 (the scale, through its link) and
// eta_1..eta_q (the partial autocorrelations r_k, each through its own
// link). With gamma_h = s(eta_0) rho_h(r) and phi = phi(r) the Levinson-Durbin
// coefficients,
//
//   phi_k^(k) = r_k,   phi_j^(k) = phi_j^(k-1) - r_k phi_{k-j}^(k-1),
//   rho_k = sum_{j<=k} phi_j^(k) rho_{k-j},   rho_h = sum_j phi_j^(q) rho_{h-j}
//
// (rho_0 = 1, the second form past the order). Step k is affine in r_k and
// the earlier coefficients do not involve it, so phi is MULTILINEAR in r:
// a derivative of phi that names some r_k twice is zero. The derivatives in
// r are carried for every multiset L of pacf indices up to the order of the
// kernel, the product phi_j rho_{h-j} differentiated by Leibniz over the
// sub-multisets B of L with no repeated index (the only ones where phi's
// derivative is not zero), with weight prod_k L_k over the indices of B.
// The links enter once, at the end: each r_k depends on eta_k alone, so
//
//   d^c_eta rho = sum_{1 <= L_k <= c_k} prod_k B_{c_k, L_k}(r_k', r_k'', ...) d^L_r rho,
//
// with B the partial Bell polynomials, and d^m_eta0 gamma = s^(m) rho. An
// order-K kernel returns the components of order K alone, for the index
// tuples it is given; the derivatives in r of lower orders are its
// ingredients. Measured against the mechanical twin ar_taylor_jet_cpp in the
// tests.

namespace {

// partial Bell polynomial B_{n,j}(x1, ..., x4), 1 <= j <= n <= 4
inline double bell(int n, int j, const double* x) {
  switch (n * 10 + j) {
  case 11: return x[1];
  case 21: return x[2];
  case 22: return x[1] * x[1];
  case 31: return x[3];
  case 32: return 3.0 * x[1] * x[2];
  case 33: return x[1] * x[1] * x[1];
  case 41: return x[4];
  case 42: return 4.0 * x[1] * x[3] + 3.0 * x[2] * x[2];
  case 43: return 6.0 * x[1] * x[1] * x[2];
  case 44: return x[1] * x[1] * x[1] * x[1];
  }
  return 0.0;
}

template <int K>
struct ArCore {
  int p, q;
  // the multisets of pacf indices of size 0..K, as count vectors
  std::vector<std::vector<int>> ms;
  std::unordered_map<std::uint64_t, int> id;
  // per multiset: its sub-multisets with no repeated index
  struct Sub { int b, rest; double w; };
  std::vector<std::vector<Sub>> sub;
  std::vector<std::vector<double>> rho;   // rho[h][multiset]
  std::vector<std::vector<double>> phi;   // phi[j - 1][multiset], order q

  std::uint64_t key(const std::vector<int>& c) const {
    std::uint64_t k = 0;
    for (int v = q - 1; v >= 0; --v) k = k * (K + 1) + (std::uint64_t) c[v];
    return k;
  }
  int find(const std::vector<int>& c) const {
    auto it = id.find(key(c));
    return it == id.end() ? -1 : it->second;
  }
  void enumerate(std::vector<int>& c, int v, int left) {
    if (v == q) {
      id[key(c)] = (int) ms.size();
      ms.push_back(c);
      return;
    }
    for (int a = 0; a <= left; ++a) {
      c[v] = a;
      enumerate(c, v + 1, left - a);
    }
    c[v] = 0;
  }

  ArCore(int p_, int q_, const NumericMatrix& seeds) : p(p_), q(q_) {
    std::vector<int> c(q, 0);
    enumerate(c, 0, K);
    const int M = (int) ms.size();
    sub.resize(M);
    for (int i = 0; i < M; ++i) {
      std::vector<int> supp;
      for (int v = 0; v < q; ++v) if (ms[i][v] > 0) supp.push_back(v);
      const int ns = (int) supp.size();
      for (int mask = 0; mask < (1 << ns); ++mask) {
        std::vector<int> b(q, 0), r = ms[i];
        double w = 1.0;
        for (int t = 0; t < ns; ++t) if (mask & (1 << t)) {
          b[supp[t]] = 1;
          r[supp[t]] -= 1;
          w *= ms[i][supp[t]];
        }
        sub[i].push_back(Sub{find(b), find(r), w});
      }
    }
    const int z = find(std::vector<int>(q, 0));
    std::vector<int> e(q, 0);
    std::vector<int> unit(q);
    for (int v = 0; v < q; ++v) { e[v] = 1; unit[v] = find(e); e[v] = 0; }

    const int H = std::max(p, q + 1);
    rho.assign(H, std::vector<double>(M, 0.0));
    rho[0][z] = 1.0;
    std::vector<std::vector<double>> old, nw;
    for (int k = 1; k <= q; ++k) {
      const double rk = seeds(k, 0);
      const int v = k - 1;
      nw.assign(k, std::vector<double>(M, 0.0));
      for (int j = 1; j < k; ++j) {
        const std::vector<double>& pj = old[j - 1];
        const std::vector<double>& pkj = old[k - j - 1];
        for (int i = 0; i < M; ++i) {
          const int lk = ms[i][v];
          if (lk >= 2) continue;
          if (lk == 1) {
            std::vector<int> r = ms[i];
            r[v] = 0;
            nw[j - 1][i] = -pkj[find(r)];
          } else {
            nw[j - 1][i] = pj[i] - rk * pkj[i];
          }
        }
      }
      nw[k - 1][z] = rk;
      if (K >= 1) nw[k - 1][unit[v]] = 1.0;
      // rho_k: the j = k term is r_k rho_0, the others by Leibniz
      std::vector<double>& rh = rho[k];
      rh[z] = rk;
      if (K >= 1) rh[unit[v]] = 1.0;
      for (int j = 1; j < k; ++j) {
        const std::vector<double>& pj = nw[j - 1];
        const std::vector<double>& rprev = rho[k - j];
        for (int i = 0; i < M; ++i) {
          double s = 0.0;
          for (const Sub& sb : sub[i]) s += sb.w * pj[sb.b] * rprev[sb.rest];
          rh[i] += s;
        }
      }
      old.swap(nw);
    }
    phi = old;
    for (int h = q + 1; h < p; ++h) {
      std::vector<double>& rh = rho[h];
      for (int j = 1; j <= q; ++j) {
        const std::vector<double>& pj = phi[j - 1];
        const std::vector<double>& rprev = rho[h - j];
        for (int i = 0; i < M; ++i) {
          double s = 0.0;
          for (const Sub& sb : sub[i]) s += sb.w * pj[sb.b] * rprev[sb.rest];
          rh[i] += s;
        }
      }
    }
  }

  // the order-K components for the given tuples (1-based free indices, 1 the
  // scale), gamma (p rows) and phi (q rows)
  List components(const NumericMatrix& seeds, const IntegerMatrix& tuples) const {
    const int nt = tuples.nrow();
    NumericMatrix g(p, nt), f(q, nt);
    double sd[5];
    for (int a = 0; a <= K; ++a) sd[a] = seeds(0, a);
    std::vector<double> rdv((size_t) q * 5, 0.0);
    for (int v = 0; v < q; ++v)
      for (int a = 1; a <= K; ++a) rdv[(size_t) v * 5 + a] = seeds(v + 1, a);
    for (int t = 0; t < nt; ++t) {
      int m = 0;
      std::vector<int> c(q, 0);
      for (int a = 0; a < K; ++a) {
        const int ix = tuples(t, a) - 1;
        if (ix == 0) ++m; else c[ix - 1] += 1;
      }
      std::vector<int> supp;
      for (int v = 0; v < q; ++v) if (c[v] > 0) supp.push_back(v);
      // every L with 1 <= L_v <= c_v on the support, 0 elsewhere
      std::vector<int> L(q, 0);
      for (int v : supp) L[v] = 1;
      std::vector<std::pair<int, double>> terms;
      while (true) {
        double w = 1.0;
        for (int v : supp) w *= bell(c[v], L[v], &rdv[(size_t) v * 5]);
        terms.push_back({find(L), w});
        size_t s = 0;
        while (s < supp.size() && L[supp[s]] == c[supp[s]]) { L[supp[s]] = 1; ++s; }
        if (s == supp.size()) break;
        L[supp[s]] += 1;
      }
      for (int h = 0; h < p; ++h) {
        double acc = 0.0;
        for (const auto& tw : terms) acc += tw.second * rho[h][tw.first];
        g(h, t) = sd[m] * acc;
      }
      for (int j = 0; j < q; ++j) {
        double acc = 0.0;
        if (m == 0) for (const auto& tw : terms) acc += tw.second * phi[j][tw.first];
        f(j, t) = acc;
      }
    }
    return List::create(Named("gamma") = g, Named("phi") = f);
  }
};

void check_seeds(int q, const NumericMatrix& seeds, int K) {
  if (seeds.nrow() != q + 1 || seeds.ncol() < K + 1)
    stop("'seeds' must have q + 1 rows and at least K + 1 columns.");
}

template <int K>
List ar_order(int p, int q, NumericMatrix seeds, IntegerMatrix tuples) {
  check_seeds(q, seeds, K);
  if (tuples.ncol() != K) stop("'tuples' must have one column per index.");
  ArCore<K> core(p, q, seeds);
  return core.components(seeds, tuples);
}

}  // namespace

// the autocovariances gamma_0..gamma_{p-1} and the coefficients phi_1..phi_q
// [[Rcpp::export]]
List ar_value_cpp(int p, int q, NumericMatrix seeds) {
  check_seeds(q, seeds, 0);
  ArCore<0> core(p, q, seeds);
  NumericVector g(p), f(q);
  for (int h = 0; h < p; ++h) g[h] = seeds(0, 0) * core.rho[h][0];
  for (int j = 0; j < q; ++j) f[j] = core.phi[j][0];
  return List::create(Named("gamma") = g, Named("phi") = f);
}

// first derivatives, one column per tuple of `tuples` (n x 1)
// [[Rcpp::export]]
List ar_d1_cpp(int p, int q, NumericMatrix seeds, IntegerMatrix tuples) {
  return ar_order<1>(p, q, seeds, tuples);
}

// second derivatives, one column per row of `tuples` (m x 2)
// [[Rcpp::export]]
List ar_d2_cpp(int p, int q, NumericMatrix seeds, IntegerMatrix tuples) {
  return ar_order<2>(p, q, seeds, tuples);
}

// third derivatives, one column per row of `tuples` (m x 3)
// [[Rcpp::export]]
List ar_d3_cpp(int p, int q, NumericMatrix seeds, IntegerMatrix tuples) {
  return ar_order<3>(p, q, seeds, tuples);
}

// fourth derivatives, one column per row of `tuples` (m x 4)
// [[Rcpp::export]]
List ar_d4_cpp(int p, int q, NumericMatrix seeds, IntegerMatrix tuples) {
  return ar_order<4>(p, q, seeds, tuples);
}

// The Toeplitz matrices M_ij = g[|i - j|, t], one per column of g, with the
// dimnames given: the assembly every derivative order of the structure ends in.
// [[Rcpp::export]]
List ar_toeplitz_cpp(NumericMatrix g, CharacterVector nm) {
  const int p = g.nrow(), nt = g.ncol();
  List out(nt);
  List dn = List::create(nm, nm);
  for (int t = 0; t < nt; ++t) {
    NumericMatrix m(p, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < p; ++i) m(i, j) = g(i > j ? i - j : j - i, t);
    m.attr("dimnames") = dn;
    out[t] = m;
  }
  return out;
}
