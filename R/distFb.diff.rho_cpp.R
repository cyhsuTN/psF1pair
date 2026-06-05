
#' @title Distribution of the Difference between Two F Beta Scores
#' @description A function to show a distribution of the difference between two F beta scores.
#' @import Rcpp
#' @importFrom stats fft median quantile rgamma

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param rho.N Correlation of two classifiers' predictions
#' among negative cases.
#' @param rho.P Correlation of two classifiers' predictions
#' among positive cases.
#' @param beta beta.

#' @examples
#' \dontrun{
#' fbscore <- F1.cond.two.b.rho(N=50, s=20,
#'                               pp1=9/11, ps1=0.75,
#'                               pp2=0.90, ps2=0.85,
#'                               rho.N=0, rho.P=0,
#'                               beta=1)
#' psF1::plotF1(fbscore, type=c("Density", "CDF")[2], xlim=c(-1,1))
#' }

#' @export
F1.cond.two.b.rho <- function(N, s,
                              pp1=9/11, ps1=0.75,
                              pp2=9/11, ps2=0.75,
                              rho.N=0, rho.P=0,
                              beta=1) {

  one.plus.beta2 <- 1 + beta^2
  s.beta2 <- s * beta^2

  pa1 <- paCalculation(N, s, pp1, ps1, beta)$pa
  pa2 <- paCalculation(N, s, pp2, ps2, beta)$pa

  nA <- N - s

  b12 <- joint_binomial_pmf_fft_sparse(
    nA, pa1, pa2, rho.N,
    prune_tol = 1e-10,
    as_sparse = TRUE
  )
  d12 <- joint_binomial_pmf_fft_sparse(
    s, ps1, ps2, rho.P,
    prune_tol = 1e-10,
    as_sparse = TRUE
  )

  i_vec <- 0:nA
  k_vec <- 0:s
  A <- outer(i_vec, k_vec, function(i, k)
    one.plus.beta2 * k / (k + i + s.beta2)
  )
  B <- A

  #sb <- Matrix::summary(b12)
  #sd <- Matrix::summary(d12)

  #if (nrow(sb) == 0L || nrow(sd) == 0L)
  #  stop("Empty support")

  # ---- C++ loop (fast) ----
  tmp <- collect_diff_weight_cpp_dense(
    b12 = as.matrix(b12),
    d12 = as.matrix(d12),
    A   = A,
    B   = B
  )

  # ---- R aggregation (exact) ----
  key <- signif(tmp$diff, 15)
  o   <- order(key)

  key_sorted <- key[o]
  w_sorted   <- tmp$weight[o]

  grp <- c(TRUE, diff(key_sorted) != 0)
  gid <- cumsum(grp)

  pf  <- rowsum(w_sorted, gid, reorder = FALSE)
  f1s <- key_sorted[grp]

  data.frame(f1s = f1s, pf1 = as.numeric(pf))

}





