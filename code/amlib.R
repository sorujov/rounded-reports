## ======================================================================
## amlib.R -- identification and inference for distributional shape from
## rounded reports.  Final implementation used for every table in the paper.
## ======================================================================
suppressWarnings(suppressMessages({
  if (!exists("em_exact_cpp")) Rcpp::sourceCpp(file.path(Sys.getenv("AMLIB_DIR", "."), "emcpp.cpp"))
}))
gcd2 <- function(a, b) { while (b > 0) { t <- b; b <- a %% b; a <- t }; a }
lcm2 <- function(a, b) a / gcd2(a, b) * b
gcdv <- function(a, b) { a <- abs(a); b <- rep_len(b, length(a)); while (any(b > 0)) { i <- b > 0; t <- b[i]; b[i] <- a[i] %% b[i]; a[i] <- t }; a }
totient <- function(n) { r <- n; p <- 2; m <- n
  while (p * p <= m) { if (m %% p == 0) { while (m %% p == 0) m <- m %/% p; r <- r - r %/% p }; p <- p + 1 }
  if (m > 1) r <- r - r %/% m; r }
divisors <- function(n) { d <- which(n %% seq_len(floor(sqrt(n))) == 0); sort(unique(c(d, n %/% d))) }
pattern_id <- function(v, levels) { id <- numeric(length(v)); for (j in seq_along(levels)) id <- id + 2^(j - 1) * (v %% levels[j] == 0); id }
pattern_law <- function(L, levels, W) { k <- seq_len(W / L) - 1; tab <- table(pattern_id(k * L, levels)); tab / sum(tab) }
em_mix <- function(counts, laws, iter = 3000, tol = 1e-11) {
  p <- rep(1 / ncol(laws), ncol(laws))
  for (it in seq_len(iter)) { num <- sweep(laws, 2, p, "*"); den <- rowSums(num); den[den <= 0] <- 1
    pn <- as.numeric(crossprod(counts, num / den)) / sum(counts)
    if (max(abs(pn - p)) < tol) { p <- pn; break }; p <- pn }
  p }

## ---- the mechanism (Proposition 1) ------------------------------------------
fit_mech <- function(y, weights = NULL, levels = c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000),
                     delta = 0.02, min_periods = 4) {
  y <- round(y); if (is.null(weights)) weights <- rep(1, length(y))
  levels <- sort(unique(levels))
  span <- diff(as.numeric(quantile(y, c(0.01, 0.99))))
  dropped <- levels[levels > span / min_periods & levels > 1]; lv <- setdiff(levels, dropped)
  W <- Reduce(lcm2, lv)
  pats <- pattern_id(y, lv); cnt <- tapply(weights, pats, sum)
  laws <- sapply(lv, function(L) { pl <- pattern_law(L, lv, W); out <- setNames(numeric(length(cnt)), names(cnt))
    cm <- intersect(names(pl), names(cnt)); out[cm] <- pl[cm]; out })
  keep <- rowSums(laws) > 0
  p <- em_mix(as.numeric(cnt[keep]), laws[keep, , drop = FALSE]); names(p) <- lv
  ## a level too coarse for the support may still be in use: compare its observed share with the fitted one
  sh <- function(l) sum(weights[y %% l == 0]) / sum(weights)
  fit_sh <- function(l) sum(p * sapply(lv, function(L) gcd2(L, l) / l))
  unident <- dropped[vapply(dropped, function(l) sh(l) - fit_sh(l) > delta, TRUE)]
  active <- lv[p >= delta & lv > 1]
  list(p = p, levels = lv, W = W, active = active, w0 = if (length(active)) Reduce(lcm2, active) else 1,
       dropped = dropped, unidentified = unident, share_below_W = sum(weights[y < W]) / sum(weights))
}

## ---- specification test on residues (Theorem 1) -----------------------------
res_setup <- function(levels, M = 10) {
  levels <- sort(unique(levels)); W <- Reduce(lcm2, levels); r <- 0:(W - 1)
  g <- gcdv(r, W); g[r == 0] <- W
  pat <- rep(0, W); for (j in seq_along(levels)) pat <- pat + 2^(j - 1) * (g %% levels[j] == 0)
  arc <- pmin(floor(r / (W / M)), M - 1); cell <- paste(pat, arc); tab <- table(cell)
  list(levels = levels, W = W, M = M, pat_of_cell = sub(" .*", "", names(tab)), cells = names(tab),
       ref = as.numeric(tab), pat_r = pat, arc_r = arc)
}
spec_test <- function(y, levels, weights = NULL, B = 999, seed = 1, min_class = 20, M = 10) {
  st <- res_setup(levels, M); y <- round(y); if (is.null(weights)) weights <- rep(1, length(y))
  r <- y %% st$W; key <- paste(st$pat_r[r + 1], st$arc_r[r + 1])
  neff <- sum(weights)^2 / sum(weights^2)
  cnt <- tapply(weights, factor(key, levels = st$cells), sum); cnt[is.na(cnt)] <- 0
  cnt <- as.numeric(cnt) / sum(weights) * neff
  cls <- split(seq_along(st$cells), st$pat_of_cell)
  stat <- function(cn) { T <- 0; df <- 0
    for (ix in cls) { NS <- sum(cn[ix]); if (length(ix) < 2 || NS < min_class) next
      e <- NS * st$ref[ix] / sum(st$ref[ix]); T <- T + sum((cn[ix] - e)^2 / e); df <- df + length(ix) - 1 }
    c(T, df) }
  s0 <- stat(cnt); set.seed(seed); Tb <- numeric(B)
  for (b in seq_len(B)) { cb <- cnt
    for (ix in cls) { NS <- round(sum(cnt[ix])); if (length(ix) < 2 || NS < min_class) next
      cb[ix] <- rmultinom(1, NS, st$ref[ix] / sum(st$ref[ix])) }
    Tb[b] <- stat(cb)[1] }
  list(W = st$W, T = s0[1], df = s0[2], T_df = if (s0[2] > 0) s0[1] / s0[2] else NA,
       p = if (s0[2] > 0) (1 + sum(Tb >= s0[1])) / (B + 1) else NA, n_eff = neff,
       share_below_W = sum(weights[y < st$W]) / sum(weights))
}

## ---- cells with boundary reports split (Lemma 2) ----------------------------
cell_split <- function(y, w, weights = NULL, tail = 0.01, top_q = 0.95) {
  if (is.null(weights)) weights <- rep(1, length(y))
  top <- as.numeric(quantile(y, top_q)); keep <- y <= top; y <- y[keep]; weights <- weights[keep]
  u <- y / w + 0.5; k <- floor(u); onb <- abs(u - round(u)) < 1e-9     # exactly on a boundary k*w - w/2
  kmax <- max(k) + 1; tab <- numeric(kmax + 2)
  add <- function(idx, wt) { t <- tapply(wt, factor(idx, levels = 0:(kmax + 1)), sum); t[is.na(t)] <- 0; as.numeric(t) }
  tab <- add(k[!onb], weights[!onb]) + add(k[onb], weights[onb] / 2) + add(pmax(k[onb] - 1, 0), weights[onb] / 2)
  p <- tab / sum(tab); ctr <- (0:(kmax + 1)) * w
  lo <- max(1, which(cumsum(p) >= tail)[1]); hi <- max(which(p > 0))
  below <- if (lo > 1) sum(p[1:(lo - 1)]) else 0; mid <- sum(p[lo:hi])
  p <- p[lo:hi]; ctr <- ctr[lo:hi]
  ## population shares: below the retained range, and in it (the top 1 - top_q was removed first)
  list(p = p / sum(p), centres = ctr, w = w, n = length(y), n_eff = sum(weights)^2 / sum(weights^2),
       share_below = top_q * below, share_mid = top_q * mid)
}

## ---- unimodal projection and the test (Theorem 3) ---------------------------
am_pava_inc <- function(y) { n <- length(y); if (!n) return(y); v <- numeric(n); ww <- numeric(n); m <- 0
  for (i in 1:n) { m <- m + 1; v[m] <- y[i]; ww[m] <- 1
    while (m > 1 && v[m - 1] > v[m]) { v[m - 1] <- (ww[m - 1] * v[m - 1] + ww[m] * v[m]) / (ww[m - 1] + ww[m]); ww[m - 1] <- ww[m - 1] + ww[m]; m <- m - 1 } }
  rep(v[1:m], ww[1:m]) }
unimodal_proj <- function(p) { K <- length(p); best <- NULL
  for (m in 1:K) { f <- numeric(K)
    if (m > 1) f[1:(m - 1)] <- am_pava_inc(p[1:(m - 1)])
    if (m < K) f[(m + 1):K] <- rev(am_pava_inc(rev(p[(m + 1):K])))
    f[m] <- max(p[m], if (m > 1) f[m - 1] else 0, if (m < K) f[m + 1] else 0)
    sse <- sum((p - f)^2); if (is.null(best) || sse < best$sse) best <- list(f = f, sse = sse, mode = m) }
  f <- pmax(best$f, 0); f / sum(f) }
pool_adjacent <- function(p, n, kappa) { K <- length(p); blk <- integer(K); b <- 1; blk[1] <- 1
  for (k in 2:K) { if (abs(p[k] - p[k - 1]) > kappa * sqrt((p[k] + p[k - 1]) / n)) b <- b + 1; blk[k] <- b }
  ave(p, blk) }
gap_at <- function(p, lam) { K <- length(p); g <- cumsum(p) - lam * seq_len(K); g0 <- c(0, g)
  left <- cummax(g0 - cummin(g0))[-1]; rg <- rev(g0); right <- rev(cummax(cummax(rg) - rg))[-1]
  t <- which.max(left[-K] + right[-1]); gl <- g[1:t]; j1 <- which.max(gl - cummin(gl))
  gr <- g[(t + 1):K]; i2 <- t + which.max(rev(cummax(rev(gr))) - gr); c(j1, i2) }
cell_test <- function(p, n_eff, B = 499, method = c("gms", "proj"), kappa = NULL, seed = NULL) {
  method <- match.arg(method); if (!is.null(seed)) set.seed(seed)
  K <- length(p); ne <- max(30, round(n_eff)); if (is.null(kappa)) kappa <- sqrt(log(ne))
  o <- em_exact_cpp(p); obs <- o[1]
  base <- if (method == "gms") pool_adjacent(p, ne, kappa) else p
  f <- unimodal_proj(base)
  st <- vapply(seq_len(B), function(b) em_exact_cpp(rmultinom(1, ne, f)[, 1] / ne)[1], 0)
  gp <- if (obs > 0) gap_at(p, o[2]) else c(NA, NA)
  list(delta = obs, p_value = (1 + sum(st >= obs - 1e-14)) / (B + 1), K = K, j1 = gp[1], i2 = gp[2])
}
## ---- the whole procedure at the certified width -----------------------------
modality_test <- function(y, weights = NULL, levels = c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000),
                          w = NULL, B = 499, delta = NULL, method = "gms", spec = TRUE, seed = 1, Kmax = 1000) {
  y <- round(y); if (is.null(weights)) weights <- rep(1, length(y))
  if (is.null(delta)) delta <- min(0.02, 2 / sqrt(sum(weights)^2 / sum(weights^2)))   # delta_n of Section 2.4
  fm <- fit_mech(y, weights, levels, delta)
  out <- list(active = paste(fm$active, collapse = "/"), w0 = fm$w0, unidentified = length(fm$unidentified) > 0,
              share_below_W = fm$share_below_W, p_hat = fm$p)

  span <- diff(as.numeric(quantile(y, c(0.01, 0.95))))
  # smallest multiple of w0 of the form {1,2,5}x10^k with at most Kmax cells over the bulk (always lattice-compatible)
  nice_w <- function(w0, lo) { if (w0 >= lo) return(w0); m <- ceiling(lo / w0); cand <- as.vector(outer(c(1, 2, 5), 10^(0:12)))
    cand <- sort(cand[cand >= m]); w0 * cand[1] }
  ww <- if (is.null(w)) nice_w(fm$w0, span / Kmax) else w
  if (spec) { lv_s <- sort(unique(c(1, levels[ww %% levels == 0]))); s <- if (length(lv_s) > 1) spec_test(y, lv_s, weights, B = 499) else NULL
    out$spec_T_df <- if (is.null(s)) NA else s$T_df; out$spec_p <- if (is.null(s)) NA else s$p
    out$spec_dist <- if (is.null(s)) NA else s$T / s$n_eff; out$spec_df <- if (is.null(s)) NA else s$df }
  if (out$unidentified && is.null(w)) { out$status <- "no grid"; out$p_value <- NA; out$K <- NA; return(out) }
  cp <- cell_split(y, ww, weights)
  if (length(cp$p) < 4) { out$status <- "too coarse"; out$p_value <- NA; out$K <- length(cp$p); return(out) }
  t <- cell_test(cp$p, cp$n_eff, B, method, seed = seed)
  out$status <- "ok"; out$w <- ww; out$K <- t$K; out$p_value <- t$p_value; out$delta <- t$delta
  if (!is.na(t$j1)) { out$gap_lo <- max(0, cp$centres[t$j1] - ww / 2); out$gap_hi <- cp$centres[t$i2] + ww / 2   # Section 3.3: cells j1..i2
    out$pstar_lo <- cp$share_below + cp$share_mid * (if (t$j1 > 1) sum(cp$p[1:(t$j1 - 1)]) else 0)
    out$pstar_hi <- cp$share_below + cp$share_mid * sum(cp$p[1:t$i2]) }
  out
}
