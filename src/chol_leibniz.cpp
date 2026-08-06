#include <Rcpp.h>
#include <string>
using namespace Rcpp;

// The Leibniz assembly of the log-Cholesky derivatives, written against the
// structure of the factor rather than as matrix products. Every derivative
// of L is a single-entry matrix E_k -- exp(eta_k) at a diagonal position,
// 1 below the diagonal, zero at every other entry, and zero outright when a
// mixed or repeated off-diagonal direction is asked -- so each Leibniz term
// is one of three cheap updates:
//
//   E_a E_b^T   one entry, v_a v_b at (r_a, r_b), and only when c_a = c_b;
//   E_k L^T     one row,   v_k L[, c_k] into row r_k;
//   L E_k^T     one column, the same vector into column r_k.
//
// The R twin (.chol_leibniz_r) builds the same components through dense
// p x p products; the two are compared at machine precision in the tests.

// [[Rcpp::export]]
List chol_leibniz_cpp(int p, IntegerMatrix tuples, IntegerVector row,
                      IntegerVector col, LogicalVector ondiag,
                      NumericMatrix L, CharacterVector free_names,
                      CharacterVector dim_names) {
    int ncomp = tuples.nrow(), order = tuples.ncol();
    List out(ncomp);
    CharacterVector nms(ncomp);
    List dn = List::create(dim_names, dim_names);

    for (int i = 0; i < ncomp; ++i) {
        NumericMatrix m(p, p);
        int nsub = 1 << order;
        for (int b = 0; b < nsub; ++b) {
            // the part of the tuple sent to the left factor
            int ks = -1, cntS = 0;
            bool ok = true;
            for (int j = 0; j < order; ++j) {
                if (!(b & (1 << j))) continue;
                int k = tuples(i, j) - 1;
                if (cntS == 0) ks = k;
                else if (k != ks) { ok = false; break; }
                ++cntS;
            }
            if (!ok || (cntS > 1 && !ondiag[ks])) continue;
            // and the complement, sent to the right factor
            int kc = -1, cntC = 0;
            for (int j = 0; j < order; ++j) {
                if (b & (1 << j)) continue;
                int k = tuples(i, j) - 1;
                if (cntC == 0) kc = k;
                else if (k != kc) { ok = false; break; }
                ++cntC;
            }
            if (!ok || (cntC > 1 && !ondiag[kc])) continue;

            if (cntS == 0 && cntC > 0) {
                int r = row[kc] - 1, c = col[kc] - 1;
                double v = ondiag[kc] ? L(r, c) : 1.0;
                for (int ii = 0; ii < p; ++ii) m(ii, r) += v * L(ii, c);
            } else if (cntC == 0 && cntS > 0) {
                int r = row[ks] - 1, c = col[ks] - 1;
                double v = ondiag[ks] ? L(r, c) : 1.0;
                for (int jj = 0; jj < p; ++jj) m(r, jj) += v * L(jj, c);
            } else if (cntS > 0 && cntC > 0) {
                int ra = row[ks] - 1, ca = col[ks] - 1;
                int rb = row[kc] - 1, cb = col[kc] - 1;
                if (ca == cb) {
                    double va = ondiag[ks] ? L(ra, ca) : 1.0;
                    double vb = ondiag[kc] ? L(rb, cb) : 1.0;
                    m(ra, rb) += va * vb;
                }
            }
        }
        m.attr("dimnames") = dn;
        out[i] = m;

        std::string nm;
        for (int j = 0; j < order; ++j) {
            if (j) nm += ":";
            nm += as<std::string>(free_names[tuples(i, j) - 1]);
        }
        nms[i] = nm;
    }
    out.attr("names") = nms;
    return out;
}
