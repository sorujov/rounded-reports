## =========================================================================
## heaping.R -- modality inference for heaped (rounded) data.
##
## Survey income and earnings are reported in round numbers, so a kernel
## density estimate has spurious bumps and the usual modality tests reject
## for reasons that have nothing to do with the shape of the latent
## distribution.  The identified object at the coarsest heaping scale is the
## vector of cell probabilities.  We therefore
##   (1) choose the coarsest grid w that carries a substantial share of the
##       reports;
##   (2) compute weighted cell probabilities on cells centred at multiples
##       of w;
##   (3) measure the excess-mass difference of the cell vector;
##   (4) calibrate it by resampling from the projection of the cell vector
##       on the cone of unimodal sequences.  This is a projection bootstrap:
##       it is valid pointwise, not uniformly over the cone, and the place
##       where it can fail is a near-flat latent density.
## =========================================================================

## ---- grid choice ---------------------------------------------------------
if (!exists("fit_rounding")) {
  .d <- if (exists("here")) here else "R"
  source(file.path(.d, "rounding.R"))
}
heap_shares <- function(y, levels = c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000),
                        w = NULL) {
  if (is.null(w)) w <- rep(1, length(y))
  sapply(levels, function(l) sum(w[abs(y / l - round(y / l)) < 1e-8]) / sum(w))
}

## The grid must absorb every rounding lattice that carries mass.  If reports
## pile up both on multiples of 10000 and on multiples of 25000, a cell width
## of 10000 still leaves a saw-tooth, because consecutive cells contain
## different numbers of 25000-points; only a common multiple of all active
## lattices gives every cell the same heap structure.  The grid is therefore
## the least common multiple of the lattices the mechanism actually uses,
## which `fit_rounding` in rounding.R estimates from the divisibility
## patterns of the reports.  This is the resolution at which the data can
## speak: valleys narrower than it are not identified, which
## `low_resolution` records.  A grid finer than that resolution is never
## used, even when it would leave more cells: the saw-tooth would come back
## and the test would reject an artefact.
choose_grid <- function(y, weights = NULL, delta = 0.02, min_cells = 15,
                        default_cells = 40, max_cells = 200, unit = 1,
                        levels = c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000)) {
  fit <- tryCatch(fit_rounding(y, weights, levels, delta = delta, unit = unit),
                  error = function(e) NULL)
  act <- if (is.null(fit)) numeric(0) else fit$active
  span <- diff(as.numeric(quantile(y, c(0.01, 0.99))))
  heaped <- length(act) > 0
  unident <- if (is.null(fit)) numeric(0) else fit$unidentified
  if (!is.null(fit) && isTRUE(fit$misfit)) unident <- c(unident, NA)
  w <- if (heaped) Reduce(lcm2, act) else span / default_cells
  if (span / w > max_cells) w <- span / max_cells
  list(w = w, p_round = if (is.null(fit)) NULL else fit$p,
       shares = if (is.null(fit)) NULL else fit$shares_obs,
       shares_fit = if (is.null(fit)) NULL else fit$shares_fit,
       lattices = act, span = span, heaped = heaped, cells = span / w,
       gof = if (is.null(fit)) NA else fit$gof,
       misfit = !is.null(fit) && isTRUE(fit$misfit),
       unidentified = unident, low_resolution = span / w < min_cells)
}

## ---- cell probabilities --------------------------------------------------
## Cells are centred at multiples of w: C_k = [kw - w/2, kw + w/2).
## The top cell pools everything above `top`; it is excluded from the
## modality analysis (its width is different) unless pool_top = FALSE.
## The modality analysis is carried out on the bulk of the distribution:
## cells are kept while the cumulative weighted mass lies between `tail`
## and 1 - `tail`, so that a handful of observations in a sparse tail
## cannot manufacture a valley.
cell_probs <- function(y, w, weights = NULL, top = NULL, tail = 0.01, top_tail = 0.05) {
  if (is.null(weights)) weights <- rep(1, length(y))
  if (is.null(top)) top <- as.numeric(quantile(y, 1 - top_tail))
  k <- floor(y / w + 0.5)
  keep <- y <= top
  tab <- tapply(weights[keep], factor(k[keep], levels = 0:max(k[keep])), sum)
  tab[is.na(tab)] <- 0
  cnt <- tapply(rep(1, sum(keep)), factor(k[keep], levels = 0:max(k[keep])), sum)
  cnt[is.na(cnt)] <- 0
  p <- as.numeric(tab) / sum(tab)
  ctr <- (0:(length(p) - 1)) * w
  lo <- max(1, which(cumsum(p) >= tail)[1])
  hi <- length(p)
  p <- p[lo:hi] / sum(p[lo:hi]); ctr <- ctr[lo:hi]
  list(p = p, counts = as.numeric(cnt)[lo:hi], centres = ctr, w = w,
       n = sum(keep), n_eff = sum(weights[keep])^2 / sum(weights[keep]^2),
       trimmed = mean(!keep))
}

## ---- unimodal projection (least squares, PAVA) --------------------------
pava_inc <- function(y, wt) {                    # isotonic (non-decreasing)
  n <- length(y); if (!n) return(y)
  val <- y; wgt <- wt; idx <- 1:n; lvl <- rep(1, n); j <- 1
  v <- numeric(n); ww <- numeric(n); sz <- integer(n); m <- 0
  for (i in 1:n) {
    m <- m + 1; v[m] <- y[i]; ww[m] <- wt[i]; sz[m] <- 1
    while (m > 1 && v[m - 1] > v[m]) {
      vv <- (ww[m - 1] * v[m - 1] + ww[m] * v[m]) / (ww[m - 1] + ww[m])
      ww[m - 1] <- ww[m - 1] + ww[m]; sz[m - 1] <- sz[m - 1] + sz[m]
      v[m - 1] <- vv; m <- m - 1
    }
  }
  rep(v[1:m], sz[1:m])
}
pava_dec <- function(y, wt) rev(pava_inc(rev(y), rev(wt)))

unimodal_fit <- function(p) {
  K <- length(p); best <- NULL
  for (m in 1:K) {
    f <- numeric(K)
    if (m > 1) f[1:(m - 1)] <- pava_inc(p[1:(m - 1)], rep(1, m - 1))
    if (m < K) f[(m + 1):K] <- pava_dec(p[(m + 1):K], rep(1, K - m))
    f[m] <- max(p[m], if (m > 1) f[m - 1] else 0, if (m < K) f[m + 1] else 0)
    sse <- sum((p - f)^2)
    if (is.null(best) || sse < best$sse) best <- list(f = f, sse = sse, mode = m)
  }
  best$f <- pmax(best$f, 0); best$f <- best$f / sum(best$f)
  best
}

## ---- excess-mass difference of a cell vector ----------------------------
## E(C, lambda) = sum_{k in C} (p_k - lambda * w).  With equal cell widths
## this is the discrete analogue of the excess mass functional; the
## two-interval maximiser gives the gap between the two modal blocks.
em_cells <- function(p, w = 1, L = 200) {
  K <- length(p); cp <- cumsum(p)
  lams <- seq(0.002, 1.05, length.out = L) * max(p) / w
  best <- -Inf; lam_b <- NA
  for (lam in lams) {
    g <- cp - lam * w * seq_len(K)
    g0 <- c(0, g)                                  # allow empty prefix
    left <- cummax(g0 - cummin(g0))[-1]
    rg <- rev(g0); right <- rev(cummax(cummax(rg) - rg))[-1]
    e1 <- left[K]
    e2 <- if (K > 1) max(left[-K] + right[-1]) else e1
    if (e2 - e1 > best) { best <- e2 - e1; lam_b <- lam }
  }
  g <- cumsum(p) - lam_b * w * seq_len(K); g0 <- c(0, g)
  left <- cummax(g0 - cummin(g0))[-1]
  rg <- rev(g0); right <- rev(cummax(cummax(rg) - rg))[-1]
  t <- which.max(left[-K] + right[-1])
  gl <- g[1:t]; j1 <- which.max(gl - cummin(gl))
  gr <- g[(t + 1):K]; i2 <- t + which.max(rev(cummax(rev(gr))) - gr)
  list(delta = max(best, 0), lambda = lam_b, j1 = j1, i2 = i2)
}

## ---- the test ------------------------------------------------------------
## H0: the latent density is unimodal.  Because binning a unimodal density
## gives a unimodal cell vector, H0 implies that p is a unimodal sequence;
## the observed p is projected on that cone and the projection is resampled.
## The projection is not the least favourable point of the cone -- for a
## cone-constrained problem that is usually the flat vector -- so the
## calibration is pointwise rather than uniform; see the discussion in the
## paper and the near-flat designs in the simulation.
heaped_modality_test <- function(y, weights = NULL, w = NULL, B = 500,
                                 top = NULL, seed = NULL,
                                 tail = 0.01, top_tail = 0.05) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(weights)) weights <- rep(1, length(y))
  gr <- if (is.null(w)) choose_grid(y, weights) else list(w = w, shares = NA)
  if (length(gr$unidentified))
    stop("no grid can be certified: either a rounding level is coarser than the support can carry, or the reports do not fit a rounding model on this ladder")
  cp <- cell_probs(y, gr$w, weights, top, tail = tail, top_tail = top_tail)
  if (length(cp$p) < 4) stop("too few cells: the data are too coarse for this test")
  p <- cp$p; K <- length(p)
  obs <- em_cells(p, cp$w)
  fit <- unimodal_fit(p)
  ne <- max(30, round(cp$n_eff))
  stat0 <- numeric(B)
  for (b in seq_len(B)) {
    pb <- rmultinom(1, ne, fit$f)[, 1] / ne
    stat0[b] <- em_cells(pb, cp$w)$delta
  }
  lo <- obs$j1; hi <- obs$i2
  list(w = gr$w, shares = gr$shares, lattices = gr$lattices,
       low_resolution = K < 15, heaped = isTRUE(gr$heaped), K = K, n = cp$n, n_eff = cp$n_eff,
       delta = obs$delta, p_value = (1 + sum(stat0 >= obs$delta)) / (B + 1),
       mode_cell = fit$mode, gap_cells = c(lo, hi),
       gap_lo = cp$centres[lo] + cp$w / 2, gap_hi = cp$centres[hi] - cp$w / 2,
       p_star_lo = sum(p[1:lo]), p_star_hi = sum(p[1:(hi - 1)]),
       valley_cells = max(0, hi - lo - 1),
       degenerate = (hi - lo <= 1),
       centres = cp$centres, p = p, fit = fit$f, trimmed = cp$trimmed)
}
