#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// The Levinson-Durbin recursion with its derivatives to fourth order,
// propagated as plain arrays. Each tracked quantity carries its value and
// the full symmetric derivative tensors with respect to the n = q + 1 free
// values; the recursion is sums and products only, so the propagation rules
// are the product rule written out per order below, term by term. Nothing
// is differenced and no generic differentiation machinery is involved.

struct T4 {
  double v;
  std::vector<double> d1, d2, d3, d4;
  explicit T4(int n)
    : v(0.0), d1(n, 0.0), d2((size_t)n * n, 0.0),
      d3((size_t)n * n * n, 0.0), d4((size_t)n * n * n * n, 0.0) {}
};

// c = a * b: Leibniz, every subset of the differentiation indices going to
// one factor and its complement to the other.
static T4 t4_mul(const T4& a, const T4& b, int n) {
  T4 c(n);
  c.v = a.v * b.v;
  for (int i = 0; i < n; ++i)
    c.d1[i] = a.d1[i] * b.v + a.v * b.d1[i];
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      size_t ij = (size_t)i * n + j;
      c.d2[ij] = a.d2[ij] * b.v + a.d1[i] * b.d1[j] + a.d1[j] * b.d1[i] +
        a.v * b.d2[ij];
    }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        size_t ijk = ((size_t)i * n + j) * n + k;
        size_t ij = (size_t)i * n + j, ik = (size_t)i * n + k,
               jk = (size_t)j * n + k;
        c.d3[ijk] = a.d3[ijk] * b.v +
          a.d2[ij] * b.d1[k] + a.d2[ik] * b.d1[j] + a.d2[jk] * b.d1[i] +
          a.d1[i] * b.d2[jk] + a.d1[j] * b.d2[ik] + a.d1[k] * b.d2[ij] +
          a.v * b.d3[ijk];
      }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l) {
          size_t ijkl = (((size_t)i * n + j) * n + k) * n + l;
          size_t ijk = ((size_t)i * n + j) * n + k,
                 ijl = ((size_t)i * n + j) * n + l,
                 ikl = ((size_t)i * n + k) * n + l,
                 jkl = ((size_t)j * n + k) * n + l;
          size_t ij = (size_t)i * n + j, ik = (size_t)i * n + k,
                 il = (size_t)i * n + l, jk = (size_t)j * n + k,
                 jl = (size_t)j * n + l, kl = (size_t)k * n + l;
          c.d4[ijkl] = a.d4[ijkl] * b.v +
            a.d3[ijk] * b.d1[l] + a.d3[ijl] * b.d1[k] +
            a.d3[ikl] * b.d1[j] + a.d3[jkl] * b.d1[i] +
            a.d2[ij] * b.d2[kl] + a.d2[ik] * b.d2[jl] +
            a.d2[il] * b.d2[jk] + a.d2[jk] * b.d2[il] +
            a.d2[jl] * b.d2[ik] + a.d2[kl] * b.d2[ij] +
            a.d1[i] * b.d3[jkl] + a.d1[j] * b.d3[ikl] +
            a.d1[k] * b.d3[ijl] + a.d1[l] * b.d3[ijk] +
            a.v * b.d4[ijkl];
        }
  return c;
}

static T4 t4_add(const T4& a, const T4& b, int n) {
  T4 c(n);
  c.v = a.v + b.v;
  for (size_t m = 0; m < a.d1.size(); ++m) c.d1[m] = a.d1[m] + b.d1[m];
  for (size_t m = 0; m < a.d2.size(); ++m) c.d2[m] = a.d2[m] + b.d2[m];
  for (size_t m = 0; m < a.d3.size(); ++m) c.d3[m] = a.d3[m] + b.d3[m];
  for (size_t m = 0; m < a.d4.size(); ++m) c.d4[m] = a.d4[m] + b.d4[m];
  return c;
}

static T4 t4_sub(const T4& a, const T4& b, int n) {
  T4 c(n);
  c.v = a.v - b.v;
  for (size_t m = 0; m < a.d1.size(); ++m) c.d1[m] = a.d1[m] - b.d1[m];
  for (size_t m = 0; m < a.d2.size(); ++m) c.d2[m] = a.d2[m] - b.d2[m];
  for (size_t m = 0; m < a.d3.size(); ++m) c.d3[m] = a.d3[m] - b.d3[m];
  for (size_t m = 0; m < a.d4.size(); ++m) c.d4[m] = a.d4[m] - b.d4[m];
  return c;
}

// A free value's image under its link: derivatives sit on one variable's
// diagonal because each depends on exactly one free value.
static T4 t4_seed(int var, const NumericVector& g, int n) {
  T4 s(n);
  s.v = g[0];
  s.d1[var] = g[1];
  s.d2[(size_t)var * n + var] = g[2];
  s.d3[((size_t)var * n + var) * n + var] = g[3];
  s.d4[(((size_t)var * n + var) * n + var) * n + var] = g[4];
  return s;
}

static void t4_pack(const T4& x, NumericMatrix& out, int row) {
  int col = 0;
  out(row, col++) = x.v;
  for (size_t m = 0; m < x.d1.size(); ++m) out(row, col++) = x.d1[m];
  for (size_t m = 0; m < x.d2.size(); ++m) out(row, col++) = x.d2[m];
  for (size_t m = 0; m < x.d3.size(); ++m) out(row, col++) = x.d3[m];
  for (size_t m = 0; m < x.d4.size(); ++m) out(row, col++) = x.d4[m];
}

// seeds: (q + 1) x 5, row 0 the scale and rows 1..q the partial
// autocorrelations, each row the link inverse and its four derivatives at
// the free value. Returns gamma (p lags) and phi (q coefficients), each row
// a packed derivative record of 1 + n + n^2 + n^3 + n^4 numbers.
// [[Rcpp::export]]
List ar_taylor_cpp(int p, int q, NumericMatrix seeds) {
  int n = q + 1;
  size_t comps = 1 + (size_t)n + (size_t)n * n + (size_t)n * n * n +
    (size_t)n * n * n * n;

  T4 scale = t4_seed(0, seeds(0, _), n);
  std::vector<T4> r;
  r.reserve(q);
  for (int k = 0; k < q; ++k) r.push_back(t4_seed(k + 1, seeds(k + 1, _), n));

  // Levinson-Durbin on the derivative records:
  //   phi_k^(k) = r_k,  phi_j^(k) = phi_j^(k-1) - r_k phi_{k-j}^(k-1),
  //   rho_k = r_k + sum_{j<k} phi_j^(k) rho_{k-j},
  // and the Yule-Walker continuation beyond the order.
  std::vector<T4> phi, rho;
  rho.reserve(p);
  T4 one(n);
  one.v = 1.0;
  rho.push_back(one);
  for (int k = 1; k <= q; ++k) {
    std::vector<T4> nw;
    nw.reserve(k);
    for (int j = 1; j < k; ++j)
      nw.push_back(t4_sub(phi[j - 1], t4_mul(r[k - 1], phi[k - j - 1], n), n));
    nw.push_back(r[k - 1]);
    phi.swap(nw);
    T4 acc = r[k - 1];
    for (int j = 1; j < k; ++j)
      acc = t4_add(acc, t4_mul(phi[j - 1], rho[k - j], n), n);
    rho.push_back(acc);
  }
  for (int h = q + 1; h <= p - 1; ++h) {
    T4 acc(n);
    for (int j = 1; j <= q; ++j)
      acc = t4_add(acc, t4_mul(phi[j - 1], rho[h - j], n), n);
    rho.push_back(acc);
  }

  NumericMatrix gamma(p, (int)comps), phim(q, (int)comps);
  for (int h = 0; h < p; ++h) t4_pack(t4_mul(scale, rho[h], n), gamma, h);
  for (int j = 0; j < q; ++j) t4_pack(phi[j], phim, j);
  return List::create(Named("gamma") = gamma, Named("phi") = phim);
}
