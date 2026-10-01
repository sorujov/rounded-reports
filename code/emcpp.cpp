#include <Rcpp.h>
using namespace Rcpp;
// E1, |C1|, E2, |C2| at level lam for cell vector p (via partial sums g_k = S_k - lam*k).
static void ek(const std::vector<double>& cp, int K, double lam, double* out) {
  int n0 = K + 1;
  std::vector<double> g(n0), lbest(n0), rbest(n0);
  std::vector<int> llen(n0), rlen(n0);
  g[0] = 0; for (int k = 1; k <= K; ++k) g[k] = cp[k-1] - lam * k;
  double cmin = g[0]; int amin = 0; double lv = 0; int lvl = 0;
  for (int i = 0; i < n0; ++i) {
    if (g[i] <= cmin) { cmin = g[i]; amin = i; }
    double gain = g[i] - cmin; int len = i - amin;
    if (gain >= lv) { lv = gain; lvl = len; }
    lbest[i] = lv; llen[i] = lvl;
  }
  double cmax = g[n0-1]; int amax = n0-1; double rv = 0; int rvl = 0;
  for (int i = n0-1; i >= 0; --i) {
    if (g[i] >= cmax) { cmax = g[i]; amax = i; }
    double gain = cmax - g[i]; int len = amax - i;
    if (gain >= rv) { rv = gain; rvl = len; }
    rbest[i] = rv; rlen[i] = rvl;
  }
  double e1 = lbest[n0-1]; int m1 = llen[n0-1];
  double e2 = -1e300; int m2 = 0;
  for (int s = 0; s < n0-1; ++s) { double t = lbest[s] + rbest[s+1]; if (t > e2) { e2 = t; m2 = llen[s] + rlen[s+1]; } }
  if (e1 >= e2) { e2 = e1; m2 = m1; }
  out[0] = e1; out[1] = m1; out[2] = e2; out[3] = m2;
}
static long NCALL = 0;
static void brk(const std::vector<double>& cp, int K, int which, double a, double b, std::vector<double>& acc, int depth) {
  if (++NCALL > 20L * (K + 10)) return;
  double fa[4], fb[4], fl[4];
  ek(cp, K, a, fa); ek(cp, K, b, fb);
  double va = fa[2*which], sa = -fa[2*which+1], vb = fb[2*which], sb = -fb[2*which+1];
  if (std::fabs(sa - sb) < 1e-12 || depth > 200) return;
  double l = (vb - va + sa * a - sb * b) / (sa - sb);
  if (!(l > a + 1e-16 && l < b - 1e-16)) return;
  ek(cp, K, l, fl);
  double tang = va + sa * (l - a);
  acc.push_back(l);
  if (std::fabs(fl[2*which] - tang) <= 1e-13 * std::max(1.0, std::fabs(tang))) return;
  brk(cp, K, which, a, l, acc, depth + 1); brk(cp, K, which, l, b, acc, depth + 1);
}
// [[Rcpp::export]]
NumericVector em_exact_cpp(NumericVector p) {
  int K = p.size(); std::vector<double> cp(K); double s = 0, mx = 0;
  for (int i = 0; i < K; ++i) { s += p[i]; cp[i] = s; if (p[i] > mx) mx = p[i]; }
  double lmax = mx * (1 + 1e-9); NCALL = 0;
  std::vector<double> cand; cand.push_back(0); cand.push_back(lmax);
  brk(cp, K, 0, 0, lmax, cand, 0); brk(cp, K, 1, 0, lmax, cand, 0);
  if (NCALL > 20L * (K + 10)) { for (int j = 1; j < 4000; ++j) cand.push_back(lmax * j / 4000.0); }
  double best = 0, lb = NA_REAL, f[4];
  for (double l : cand) { ek(cp, K, l, f); double d = f[2] - f[0]; if (d > best) { best = d; lb = l; } }
  return NumericVector::create(best, lb);
}
