
F1.cond.b <- function(N, s, pp=9/11, ps=0.75, beta=1) {

  #pa <- paCalculation_v1(N, s, pp, ps)$pa
  pa <- paCalculation(N, s, pp, ps, beta)$pa

  d <- 0:s
  b <- 0:(N-s)
  d_probs <- dbinom(d, s, prob=ps)
  b_probs <- dbinom(b, N-s, prob=pa)

  idx.d <- (d_probs > 1E-10)
  idx.b <- (b_probs > 1E-10)

  d <- d[idx.d]
  b <- b[idx.b]

  one.plus.beta2 <- (1+beta^2)
  s.beta2 <- s*beta^2

  val1 <- outer(d, b, function(d, b) ifelse(d == 0, 0, one.plus.beta2 * d / (d + b + s.beta2)))
  val2 <- outer(d, b, function(d, b) d_probs[d + 1] * b_probs[b + 1])

  pf1 <- tapply(as.vector(val2), as.vector(val1), sum)
  f1s <- as.numeric(names(pf1))
  return(data.frame(f1s=f1s, pf1=as.numeric(pf1) ))

}

paCalculation <- function(N, s, pp=9/11, ps=0.75, beta=1) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)
  s.beta2 <- s*beta2

  Fb <- one.plus.beta2*pp*ps/(pp*beta2+ps)
  if(is.na(Fb)) Fb <- 0

  if(ps>1 | pp>1) {
    stop(paste0("pp > 1 or ps > 1"))
  }

  if(ps<0 | pp<0) {
    stop(paste0("pp < 0 or ps < 0"))
  }

  d <- 0:s
  d_probs <- dbinom(d, s, prob=ps)
  idx.d <- d_probs > 1E-10
  d <- d[idx.d]

  test1 <- function(pa) {
    b <- 0:(N-s)
    b_probs <- dbinom(b, N-s, prob=pa)

    idx.b <- b_probs > 1E-10
    b <- b[idx.b]

    cond.pre <- sum(outer(d, b, function(d, b) {
      ifelse(d == 0, 0, one.plus.beta2 * d / (d + b + s.beta2) * d_probs[d + 1] * b_probs[b + 1])
    }))

    cond.pre
  }

  if(pp == 1) {
    pa <- 0
    details <- NULL
  } else if(test1(0) < Fb) {
    print(paste0("Warning: no value of pa satisfies Fb with ", "pp=", round(pp,2), " and ps=", round(ps,2),
                 ". Force pa = 0."))
    pa <- 0; details <- NULL
  } else if(test1(1) > Fb) {
    print(paste0("Warning: no value of pa satisfies Fb with ", "pp=", round(pp,2), " and ps=", round(ps,2),
                 ". Force pa = 1."))
    pa <- 1; details <- NULL
  } else {
    f.pa <- function(b) test1(b) - Fb
    solve.pa <- uniroot(f.pa, interval = c(0, 1))
    pa <- solve.pa$root
    details <- solve.pa
  }

  list(pa = pa, details = details)
}



dFb_variance_approx <- function(N, s,
                                pp1=9/11, ps1=0.75,
                                pp2=9/11, ps2=0.75,
                                rho.N=0, rho.P=0,
                                beta=1) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)
  s.beta2 <- s*beta2

  pa1 <- paCalculation(N, s, pp1, ps1, beta)$pa
  pa2 <- paCalculation(N, s, pp2, ps2, beta)$pa

  nA <- N - s

  w_fun <- function(x1, x2, y1, y2) {
    den1 <- x1 + y1 + beta2
    den2 <- x2 + y2 + beta2
    one.plus.beta2*y1/den1 - one.plus.beta2*y2/den2
  }

  partial_w_fun <- function(x1, x2, y1, y2) {
    den1 <- (x1 + y1 + beta2)^2
    den2 <- (x2 + y2 + beta2)^2
    partial_w <- c(-one.plus.beta2*y1/den1,
                   one.plus.beta2*y2/den2,
                   one.plus.beta2*(x1+beta2)/den1,
                   -one.plus.beta2*(x2+beta2)/den2)
    partial_w
  }

  ss1 <- sqrt( c(pa1*(1-pa1), pa2*(1-pa2)) )
  ss2 <- sqrt( c(ps1*(1-ps1), ps2*(1-ps2)) )
  Sigma1 <- (nA/s) * (ss1 %o% ss1) * ((1-rho.N)*diag(2) + rho.N)
  Sigma2 <- (ss2 %o% ss2) * ((1-rho.P)*diag(2) + rho.P)
  Sigma <- Matrix::bdiag(Sigma1, Sigma2)

  Ew <- w_fun(x1=nA/s*pa1, x2=nA/s*pa2, y1=ps1, y2=ps2)
  partial_w <- partial_w_fun(x1=nA/s*pa1, x2=nA/s*pa2, y1=ps1, y2=ps2)
  var_w <- as.numeric(t(partial_w) %*% Sigma %*% partial_w / s)

  list(Ew=Ew, Varw=var_w)

  #out <- data.frame(f1s = f1s, pf1 = as.numeric(pf))
  #(dmean <- sum(out$f1s * out$pf1))
  #(dvar <- sum(out$f1s^2 * out$pf1) - dmean^2)

}


##########################
# --- helpers: 2D FFT / IFFT ---------------------------------------------------
# Forward 2D FFT
fft2 <- function(X) {
  # FFT along columns, then along rows
  Xc <- apply(X, 2, fft)                 # along rows (each column)
  t(apply(Xc, 1, fft))                   # along columns (each row); keep dims
}

# Inverse 2D FFT (normalized)
ifft2 <- function(F) {
  Xc <- apply(F, 2, function(col) fft(col, inverse = TRUE))
  Xr <- t(apply(Xc, 1, function(row) fft(row, inverse = TRUE)))
  Xr / (nrow(F) * ncol(F))
}

# --- new: joint pmf via 2D FFT with sparse output -----------------------------
# Returns an (n+1) x (n+1) matrix; by default sparse dgCMatrix.
# You can control pruning via prune_tol; negatives below tol_neg are clamped to 0.
joint_binomial_pmf_fft_sparse <- function(
    n, p1, p2, rho,
    prune_tol = 1e-10,     # set tiny coefficients to zero (sparsify)
    tol_neg   = 1e-10,     # clamp small negative numerical noise
    pad_fast  = TRUE,      # pad to "fast" FFT length (stats::nextn)
    as_sparse = TRUE       # return a dgCMatrix if TRUE, else dense numeric matrix
) {
  # ---- joint Bernoulli cell probs ----
  pi11 <- p1 * p2 + rho * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  pi10 <- p1 - pi11
  pi01 <- p2 - pi11
  pi00 <- 1 - p1 - p2 + pi11
  probs <- c(pi11, pi10, pi01, pi00)

  if (any(probs < -tol_neg)) {
    stop("Invalid rho: produces negative joint probabilities (beyond tol_neg).")
  }
  # clamp small negatives (numerical)
  probs[probs < 0] <- 0
  pi11 <- probs[1]; pi10 <- probs[2]; pi01 <- probs[3]; pi00 <- probs[4]

  # target size
  base_m <- n + 1L
  if (pad_fast) {
    # Next fast length for FFT (has small prime factors); prevents wrap-around
    M <- stats::nextn(base_m)  # requires stats (loaded by default)
  } else {
    M <- base_m
  }

  # Build kernel K on MxM grid: entries at (0,0), (1,0), (0,1), (1,1)
  K <- matrix(0, M, M)
  K[1, 1] <- pi00
  K[2, 1] <- pi10
  K[1, 2] <- pi01
  K[2, 2] <- pi11

  # n-fold convolution via FFT
  KF   <- fft2(K)
  KF_n <- KF^n
  pmf_full <- Re(ifft2(KF_n))  # numerical real part

  # Trim to (n+1)x(n+1)
  pmf <- pmf_full[seq_len(base_m), seq_len(base_m)]

  # Numerical cleanup: clamp tiny negatives, normalize
  pmf[pmf < 0 & pmf > -tol_neg] <- 0
  s <- sum(pmf)
  if (abs(s - 1) > 1e-10) {
    # Small drift; renormalize to be safe
    pmf <- pmf / s
  }

  # Prune tiny coefficients and convert to sparse if requested
  if (as_sparse || prune_tol > 0) {
    pmf[pmf < prune_tol] <- 0
  }

  if (as_sparse) {
    # return sparsified dgCMatrix
    Matrix::Matrix(pmf, sparse = TRUE)
  } else {
    pmf
  }
}


joint_binomial_pmf_all_fast <- function(n, p1, p2, rho, tol = 1e-12) {

  # ----- Compute joint Bernoulli probabilities -----
  pi11 <- p1*p2 + rho * sqrt(p1*(1-p1)*p2*(1-p2))
  pi10 <- p1 - pi11
  pi01 <- p2 - pi11
  pi00 <- 1 - p1 - p2 + pi11

  probs <- c(pi11, pi10, pi01, pi00)
  if (any(probs < -tol)) stop("Invalid rho: produces negative joint probabilities")

  # Clamp tiny negatives to 0
  probs[probs < 0] <- 0
  pi11 <- probs[1]; pi10 <- probs[2]; pi01 <- probs[3]; pi00 <- probs[4]

  # ----- Dynamic programming -----
  pmf <- matrix(0, n + 1, n + 1)
  pmf[1, 1] <- 1

  for (t in 1:n) {
    new <- matrix(0, n + 1, n + 1)

    # old support is 0:(t-1), i.e. indices 1:t
    old <- pmf[1:t, 1:t]

    # (0,0)
    new[1:t, 1:t] <- new[1:t, 1:t] + pi00 * old
    # (1,0) shifts x
    new[2:(t+1), 1:t] <- new[2:(t+1), 1:t] + pi10 * old
    # (0,1) shifts y
    new[1:t, 2:(t+1)] <- new[1:t, 2:(t+1)] + pi01 * old
    # (1,1) shifts both
    new[2:(t+1), 2:(t+1)] <- new[2:(t+1), 2:(t+1)] + pi11 * old

    pmf <- new
  }

  rownames(pmf) <- 0:n
  colnames(pmf) <- 0:n
  pmf
}


pmulti <- function(p1, p2, rho) {
  pi11 <- p1*p2 + rho * sqrt(p1*(1-p1)*p2*(1-p2))
  pi10 <- p1 - pi11
  pi01 <- p2 - pi11
  pi00 <- 1 - p1 - p2 + pi11
  c(pi11, pi10, pi01, pi00)
}



# Bayesian phi (φ) correlation for a 2x2 table using Dirichlet posterior smoothing
# - Accepts raw binary vectors x,y OR cell counts a,b,c,d
# - Uses Jeffreys prior by default (alpha = 0.5), good for small samples
# - Returns posterior draws, summaries, and optional plots
bayes_phi <- function(x = NULL, y = NULL,
                      a = NULL, b = NULL, c = NULL, d = NULL,
                      alpha = 0.5,         # Dirichlet prior concentration per cell (0.5 = Jeffreys)
                      ndraws = 20000,      # posterior draws
                      seed = 123) {

  stopifnot(is.numeric(alpha), alpha > 0)
  set.seed(seed)

  # ---- 1) Get counts a,b,c,d ----
  if (!is.null(x) && !is.null(y)) {
    if (length(x) != length(y)) stop("x and y must have the same length.")
    if (!all(x %in% c(0,1), na.rm = TRUE)) stop("x must be binary in {0,1}.")
    if (!all(y %in% c(0,1), na.rm = TRUE)) stop("y must be binary in {0,1}.")
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]

    a <- sum(x == 1 & y == 1)  # p11
    b <- sum(x == 1 & y == 0)  # p10
    c <- sum(x == 0 & y == 1)  # p01
    d <- sum(x == 0 & y == 0)  # p00
  } else {
    # Using counts directly
    if (any(vapply(list(a,b,c,d), is.null, logical(1)))) {
      stop("Provide either (x, y) OR the four counts (a, b, c, d).")
    }
    if (any(c(a,b,c,d) < 0) || any(!is.finite(c(a,b,c,d)))) {
      stop("Counts a,b,c,d must be nonnegative finite numbers.")
    }
  }

  counts <- c(a = as.numeric(a), b = as.numeric(b), c = as.numeric(c), d = as.numeric(d))
  n <- sum(counts)
  if (n == 0) stop("Total count n is zero.")

  # ---- 2) Posterior for cell probabilities is Dirichlet(alpha + counts) ----
  post_alpha <- counts + alpha

  # Dirichlet sampler
  rdirichlet <- function(n, alpha_vec) {
    k <- length(alpha_vec)
    out <- matrix(NA_real_, nrow = n, ncol = k)
    for (j in 1:k) out[, j] <- rgamma(n, shape = alpha_vec[j], rate = 1)
    out <- out / rowSums(out)
    colnames(out) <- names(alpha_vec)
    out
  }

  P <- rdirichlet(ndraws, post_alpha)  # columns: a b c d (p11, p10, p01, p00)
  p11 <- P[, "a"]; p10 <- P[, "b"]; p01 <- P[, "c"]; p00 <- P[, "d"]

  # ---- 3) Compute phi for each draw ----
  # Marginal probabilities
  px1 <- p11 + p10
  px0 <- p01 + p00
  py1 <- p11 + p01
  py0 <- p10 + p00

  num <- p11 * p00 - p10 * p01
  den <- sqrt(px1 * px0 * py1 * py0)

  phi <- num / den

  # Numerical guard for extremely rare cases where den ~ 0
  phi[!is.finite(phi)] <- NA_real_

  # ---- 4) Summaries ----
  phi_mean   <- mean(phi, na.rm = TRUE)
  phi_median <- median(phi, na.rm = TRUE)
  phi_ci     <- quantile(phi, probs = c(0.025, 0.975), na.rm = TRUE)

  # Posterior for margins (handy to inspect)
  px1_mean <- mean(px1); py1_mean <- mean(py1)

  res <- list(
    input_counts = counts,
    n = n,
    prior_alpha_per_cell = alpha,
    posterior_alpha = post_alpha,
    #draws = data.frame(phi = phi, p11 = p11, p10 = p10, p01 = p01, p00 = p00,
    #                   px1 = px1, py1 = py1),
    summary = list(
      phi_mean = unname(phi_mean),
      phi_median = unname(phi_median),
      phi_ci95 = unname(phi_ci),
      px1_mean = px1_mean,
      py1_mean = py1_mean
    )
  )

  class(res) <- "bayes_phi_result"
  return(res)
}
