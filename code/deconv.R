## ---------------------------------------------------------------------------
## The mechanism-efficient route: the latent histogram at a width finer than
## the certified one.
##
## Section 2 identifies the reporting mechanism without touching the income
## distribution.  Once it is in hand the report distribution is a known linear
## map of the latent histogram at any width w that divides every active level,
##
##      q = A(w) p,      A(w)[v, k] = sum_l pi_l(x_k) 1{ the level-l interval
##                                    containing cell k reports v },
##
## and A(w) is square and, for the ladders that occur in practice, extremely
## well conditioned: for the ACS ladder at w = 100 dollars over a thousand
## cells its condition number is about 50, under floor and nearest rounding
## alike.  So p is point identified at w = 100 rather than at the least common
## multiple of 50,000, and the price is not ill-posedness but a variance
## inflation of order kappa(A).  The certified width of Section 3 is the
## endpoint of this at which A is the identity on coarse cells and nothing at
## all need be known about pi; every finer width buys resolution with an
## assumption.  Reporting kappa(w) alongside w is what makes the trade
## explicit.
## ---------------------------------------------------------------------------

## A(w).  `pi_fun(x)` returns a length(x) by length(levels) matrix of rounding
## probabilities.  A cell of width w need not sit inside one rounding interval
## -- when w does not divide the level, a cell straddles two of them, and the
## mass is split in proportion to the overlap.  That is exact whenever the
## density is constant within a cell, which is the same local-smoothness
## condition the divisibility argument of Section 2 already requires, and
## which is innocuous at the widths used here.  When w divides the level the
## overlaps are 0 or 1 and nothing is assumed at all.
mech_matrix <- function(levels, pi_fun, K, w, mode = c("nearest", "floor")) {
  mode <- match.arg(mode)
  centres <- (seq_len(K) - 0.5) * w
  P <- pi_fun(centres)
  stopifnot(ncol(P) == length(levels))
  A <- matrix(0, K, K)
  for (j in seq_along(levels)) {
    L <- levels[j]
    for (k in seq_len(K)) {
      a <- (k - 1) * w; b <- k * w          # the cell
      ## the rounding intervals of level L that meet [a, b)
      lo <- if (mode == "floor") floor(a / L) else floor(a / L + 0.5)
      hi <- if (mode == "floor") floor((b - 1e-9) / L) else floor((b - 1e-9) / L + 0.5)
      for (m in lo:hi) {
        ia <- if (mode == "floor") m * L else (m - 0.5) * L
        ib <- ia + L
        ov <- max(0, min(b, ib) - max(a, ia)) / w
        if (ov <= 0) next
        v <- m * L                          # the report this interval produces
        vi <- min(max(as.integer(floor(v / w)) + 1L, 1L), K)
        A[vi, k] <- A[vi, k] + P[k, j] * ov
      }
    }
  }
  A
}

kappa_of <- function(A) tryCatch(kappa(A, exact = TRUE), error = function(e) NA_real_)

## min ||A p - q||^2 over the simplex, by projected gradient.  The problem is
## small, convex and well conditioned, so a few hundred steps are plenty.
simplex_ls <- function(A, q, iter = 3000, tol = 1e-12) {
  K <- ncol(A)
  AtA <- crossprod(A); Atq <- crossprod(A, q)
  step <- 1 / max(eigen(AtA, symmetric = TRUE, only.values = TRUE)$values)
  p <- rep(1 / K, K); prev <- Inf
  for (it in seq_len(iter)) {
    p <- proj_simplex(p - step * (AtA %*% p - Atq))
    if (it %% 50 == 0) {
      obj <- sum((A %*% p - q)^2)
      if (abs(prev - obj) < tol * max(1, obj)) break
      prev <- obj
    }
  }
  as.numeric(p)
}

## Euclidean projection onto the probability simplex (Duchi et al. 2008).
proj_simplex <- function(v) {
  u <- sort(as.numeric(v), decreasing = TRUE)
  cs <- cumsum(u) - 1
  rho <- max(which(u - cs / seq_along(u) > 0))
  pmax(as.numeric(v) - cs[rho] / rho, 0)
}

## The deconvolved histogram of reports `y` at width w.
deconv_hist <- function(y, w, levels, pi_fun, weights = NULL,
                        mode = "nearest", K = NULL) {
  if (is.null(weights)) weights <- rep(1, length(y))
  if (is.null(K)) K <- as.integer(ceiling(max(y) / w)) + 1L
  idx <- pmin(pmax(as.integer(floor(y / w)) + 1L, 1L), K)
  q <- as.numeric(tapply(weights, factor(idx, levels = seq_len(K)), sum))
  q[is.na(q)] <- 0; q <- q / sum(q)
  A <- mech_matrix(levels, pi_fun, K, w, mode)
  list(p = simplex_ls(A, q), q = q, A = A, w = w, K = K,
       kappa = kappa_of(A), n_eff = sum(weights)^2 / sum(weights^2))
}

## The mechanism-efficient test.  Identical in form to the cell test of
## Section 3 -- excess mass from the cone of unimodal sequences, calibrated by
## a projection bootstrap -- but applied to the deconvolved histogram at a
## width the mechanism buys rather than to the raw cells at the certified one.
## The null is resampled on the report scale, where the sampling is, and then
## pushed back through the same inversion, so the bootstrap carries the
## variance the inversion adds.
deconv_test <- function(y, w, levels, pi_fun, weights = NULL, B = 499,
                        mode = "nearest", K = NULL, tail = 0.01, top_tail = 0.05) {
  d <- deconv_hist(y, w, levels, pi_fun, weights, mode, K)
  keep <- which(cumsum(d$p) > tail & rev(cumsum(rev(d$p))) > top_tail)
  if (length(keep) < 4) return(list(K = length(keep), p_value = NA_real_, w = w))
  p <- d$p[keep]; p <- p / sum(p)
  obs <- em_cells(p, w)
  fit <- unimodal_fit(p)
  ## the report distribution the null implies, on the full support
  pnull <- d$p; pnull[keep] <- fit$f * sum(d$p[keep]); pnull <- pnull / sum(pnull)
  qnull <- as.numeric(d$A %*% pnull); qnull <- pmax(qnull, 0); qnull <- qnull / sum(qnull)
  ne <- max(30, round(d$n_eff))
  st <- replicate(B, {
    qs <- rmultinom(1, ne, qnull)[, 1] / ne
    ps <- simplex_ls(d$A, qs)
    ps <- ps[keep]; s <- sum(ps)
    if (s <= 0) NA_real_ else em_cells(ps / s, w)$delta
  })
  list(w = w, K = length(keep), kappa = d$kappa, delta = obs$delta,
       p_value = (1 + sum(st >= obs$delta, na.rm = TRUE)) / (sum(!is.na(st)) + 1),
       pstar_lo = sum(p[1:obs$j1]), pstar_hi = sum(p[1:(obs$i2 - 1)]),
       n_eff = d$n_eff)
}

## The mechanism as a function of the amount, from a fitted logit.  The fit
## returns the rounding probabilities at the observed reports; A(w) needs them
## at the cell centres, so they are carried across by interpolation in log
## amount, which is the scale the logit is linear on.
pi_from_fit <- function(fit, y) {
  ok <- is.finite(y) & y > 0
  ly <- log(y[ok]); P <- fit$p_x
  if (nrow(P) != length(ly)) {
    m <- min(nrow(P), length(ly)); ly <- ly[seq_len(m)]; P <- P[seq_len(m), , drop = FALSE]
  }
  o <- order(ly); ly <- ly[o]; P <- P[o, , drop = FALSE]
  function(x) {
    lx <- log(pmax(x, exp(min(ly))))
    out <- vapply(seq_len(ncol(P)), function(j)
      approx(ly, P[, j], xout = lx, rule = 2)$y, numeric(length(lx)))
    out <- matrix(pmax(out, 0), nrow = length(lx))
    out / pmax(rowSums(out), 1e-12)
  }
}
