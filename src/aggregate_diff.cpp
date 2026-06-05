#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List collect_diff_weight_cpp_dense(
    NumericMatrix b12,
    NumericMatrix d12,
    NumericMatrix A,
    NumericMatrix B
) {
  // Dimensions
  const int nA = b12.nrow() - 1;   // i,j = 0..nA
  const int s  = d12.nrow() - 1;   // k,l = 0..s

  // Worst-case allocation (safe upper bound)
  const R_xlen_t maxN =
    (R_xlen_t)(nA + 1) * (nA + 1) * (s + 1) * (s + 1);

  NumericVector diff(maxN);
  NumericVector weight(maxN);

  R_xlen_t idx = 0;

  // ---- Main loops ----
  for (int i = 0; i <= nA; ++i) {
    for (int j = 0; j <= nA; ++j) {

      const double bij = b12(i, j);
      if (bij == 0.0) continue;        // ✅ safe skip

      for (int k = 0; k <= s; ++k) {
        const double Aik = A(i, k);

        for (int l = 0; l <= s; ++l) {
          const double dkl = d12(k, l);
          if (dkl == 0.0) continue;    // ✅ safe skip

          diff[idx]   = Aik - B(j, l);
          weight[idx] = bij * dkl;
          ++idx;
        }
      }
    }
  }

  // ---- Trim to used length ----
  return List::create(
    Named("diff")   = diff[Range(0, idx - 1)],
                          Named("weight") = weight[Range(0, idx - 1)]
  );
}
