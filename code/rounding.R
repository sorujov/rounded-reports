## =====================================================================
## Identifying the rounding mechanism from the divisibility pattern.
##
## A reported value is the latent value rounded to one of a known ladder
## of levels L (1 = no rounding beyond the recording unit).  Which level
## was used is not observed; what is observed, for each report, is which
## candidate levels divide it.  Given L the report is L*k, and over a
## support covering many periods of W = lcm(levels) the multiplier k is
## close to uniform modulo W/L, so the distribution of the divisibility
## pattern given L is known exactly and free of the latent density.  The
## mixing weights p_L are then the rounding probabilities and are
## estimated by maximum likelihood (EM) over the pattern counts.
##
## This replaces share rules and residue-class ratios, both of which
## confuse arithmetic coincidence with rounding: once reports pile up on
## multiples of 5000, a fifth of them are multiples of 25000 whatever the
## respondents did, and the mixture accounts for that by construction
## because the 5000 component already places that mass there.
## =====================================================================

gcd2 <- function(a, b) { while (b > 0) { t <- b; b <- a %% b; a <- t }; a }
lcm2 <- function(a, b) a / gcd2(a, b) * b

## The pattern of an integer: which candidate levels divide it, encoded as
## the integer whose bits are those indicators.  (Encoding rather than
## pasting keeps this linear in the number of reports, which matters for
## administrative files with hundreds of thousands of records.)
pattern_id <- function(v, levels) {
  id <- numeric(length(v))
  for (j in seq_along(levels)) id <- id + 2^(j - 1) * (v %% levels[j] == 0)
  id
}

## Distribution of the pattern under rounding to level L, exactly.
pattern_law <- function(L, levels, W) {
  k <- seq_len(W / L) - 1
  tab <- table(pattern_id(k * L, levels))
  tab / sum(tab)
}

## EM for the mixing weights.  `counts` is a named vector of observed
## pattern frequencies (weights allowed); `laws` the matrix of component
## pattern laws, patterns in rows, levels in columns.
em_mix <- function(counts, laws, iter = 500, tol = 1e-10) {
  p <- rep(1 / ncol(laws), ncol(laws))
  for (it in seq_len(iter)) {
    num <- sweep(laws, 2, p, "*")
    den <- rowSums(num); den[den <= 0] <- 1
    post <- num / den
    pn <- as.numeric(crossprod(counts, post)) / sum(counts)
    if (max(abs(pn - p)) < tol) { p <- pn; break }
    p <- pn
  }
  p
}

## Fit the mechanism.  Returns the estimated rounding probabilities, the
## implied and observed shares (a fit diagnostic), and the levels whose
## probability exceeds `delta`.
fit_rounding <- function(y, weights = NULL,
                         levels = c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000),
                         delta = 0.02, unit = 1, max_W = 2e6,
                         min_periods = 4) {
  if (is.null(weights)) weights <- rep(1, length(y))
  levels <- sort(unique(levels))
  ## Rounding to a level the data barely span is not identified: the
  ## multiplier k would take a handful of values and its distribution would
  ## be nothing like uniform, so such levels are dropped from the ladder.
  span <- diff(as.numeric(quantile(y, c(0.01, 0.99))))
  dropped <- levels[levels > span / min_periods]
  levels <- levels[levels <= span / min_periods]
  if (!length(levels)) stop("no candidate level is identified on this support")
  v <- round(y / unit)
  lv <- levels / unit
  if (any(abs(lv - round(lv)) > 1e-8)) stop("levels must be multiples of `unit`")
  lv <- round(lv)
  ## The identified resolution is the least common multiple of the ladder,
  ## and for a fine recording unit that can be enormous.  Rather than refuse
  ## the whole ladder, drop its coarsest members until the multiple is
  ## tractable; a dropped level is not ignored, it is moved to `dropped` and
  ## its observed divisibility share is checked against the fit below, so a
  ## level that is genuinely in use still blocks certification.
  W <- Reduce(lcm2, lv)
  while (W > max_W && length(lv) > 1L) {
    dropped <- c(dropped, levels[length(lv)])
    lv <- lv[-length(lv)]; levels <- levels[-length(levels)]
    W <- Reduce(lcm2, lv)
  }
  if (W > max_W) stop("candidate ladder too wide: lcm exceeds max_W")
  pats <- pattern_id(v, lv)
  cnt <- tapply(weights, pats, sum)
  laws <- sapply(lv, function(L) {
    pl <- pattern_law(L, lv, W)
    out <- setNames(numeric(length(cnt)), names(cnt))
    common <- intersect(names(pl), names(cnt))
    out[common] <- pl[common]
    out
  })
  dimnames(laws) <- list(names(cnt), as.character(levels))
  ## patterns the model cannot produce (rounding to a level outside the
  ## ladder, or a latent atom) carry no information about the mixture
  keep <- rowSums(laws) > 0
  p <- em_mix(as.numeric(cnt[keep]), laws[keep, , drop = FALSE])
  names(p) <- as.character(levels)
  ## fitted against observed shares, for every candidate including those the
  ## support cannot identify
  all_lv <- sort(unique(c(levels, dropped))) / unit
  s_obs <- heap_shares(y, all_lv * unit, weights)
  A <- outer(all_lv, lv, Vectorize(function(l, L) gcd2(L, l) / l))
  s_fit <- as.numeric(A %*% p)
  names(s_obs) <- names(s_fit) <- as.character(all_lv * unit)
  ## A level too coarse to be identified on this support may still be in use.
  ## If its observed share exceeds what the fitted mechanism accounts for by
  ## more than `delta`, the reporting cannot be pinned down and no grid can
  ## be certified; the caller is told so rather than given a grid that may be
  ## finer than the rounding.
  unident <- if (!length(dropped)) numeric(0) else {
    j <- match(as.character(dropped), names(s_obs))
    dropped[s_obs[j] - s_fit[j] > delta]
  }
  ## Goodness of fit.  If the mixture cannot reproduce the divisibility
  ## shares, the reports are not a rounded version of a continuous latent
  ## variable on this ladder at all -- equivalised household income, a sum of
  ## rounded amounts divided by the square root of household size, is the
  ## leading example -- and no grid should be certified from it.
  j <- match(as.character(levels), names(s_obs))
  gof <- max(abs(s_obs[j] - s_fit[j]))
  tol <- max(0.02, 3 / sqrt(sum(weights)^2 / sum(weights^2)))
  list(p = p, unexplained = 1 - sum(cnt[keep]) / sum(cnt),
       shares_obs = s_obs, shares_fit = s_fit, dropped = dropped,
       gof = gof, gof_tol = tol, misfit = gof > tol,
       unidentified = unident,
       active = levels[p >= delta & levels > unit], delta = delta, W = W * unit)
}

## ---- magnitude-dependent reporting ---------------------------------------
## Rounding is coarser at higher amounts: a respondent earning 180 thousand
## answers to the nearest 10 or 25 thousand, one earning 22 thousand to the
## nearest 500.  Because the divisibility law of Section 2 conditions only on
## the level, the mixture can simply be fitted band by band, and the bands
## are quantile bands so that each carries the same number of reports.  The
## grid must then absorb every level active in any band.
fit_rounding_bands <- function(y, weights = NULL, levels = NULL, n_bands = 5,
                               delta = 0.02, ...) {
  if (is.null(weights)) weights <- rep(1, length(y))
  if (is.null(levels)) levels <- default_levels(y)
  qs <- as.numeric(quantile(y, seq(0, 1, length.out = n_bands + 1)))
  qs[1] <- -Inf; qs[length(qs)] <- Inf
  b <- cut(y, qs, labels = FALSE, right = FALSE)
  fits <- lapply(seq_len(n_bands), function(k) {
    i <- which(b == k)
    if (length(i) < 200) return(NULL)
    tryCatch(fit_rounding(y[i], weights[i], levels, delta = delta, ...),
             error = function(e) NULL)
  })
  ok <- !vapply(fits, is.null, TRUE)
  act <- sort(unique(unlist(lapply(fits[ok], `[[`, "active"))))
  P <- matrix(NA_real_, n_bands, length(levels),
              dimnames = list(NULL, as.character(levels)))
  for (k in which(ok)) P[k, names(fits[[k]]$p)] <- fits[[k]]$p
  list(bands = qs, band = b, fits = fits, p_band = P, active = act,
       gof = if (any(ok)) max(vapply(fits[ok], `[[`, 0, "gof")) else NA,
       misfit = any(vapply(fits[ok], function(f) isTRUE(f$misfit), TRUE)),
       levels = levels)
}

## A decimal ladder anchored on the magnitude of the median: powers of ten
## and their halves and quarters, four orders below the median up to it.
default_levels <- function(y) {
  med <- median(y[y > 0]); e <- floor(log10(med))
  lv <- sort(unique(round(as.numeric(outer(c(1, 2.5, 5), 10^(seq(max(0, e - 4), e)))))))
  unique(c(1, lv[lv >= 1 & lv <= med]))
}

## ---- the mechanism as a function of the amount ---------------------------
## Fitting band by band runs into an identification problem of its own: a
## band spanning eighteen thousand cannot speak about rounding to fifty
## thousand.  The whole sample can, if the level probabilities are allowed to
## move smoothly with the amount.  We take a multinomial logit,
##   pi_l(x) = exp(a_l + b_l log x) / sum_m exp(a_m + b_m log x),
## and fit it by EM: the E-step forms the posterior over levels from the
## divisibility pattern, which is exactly the known law A of Section 2, and
## the M-step is a weighted multinomial logit on those posteriors.  The
## latent density never enters, because A does not depend on it.
## ---- a small multinomial logit -------------------------------------------
## The M-step needs a weighted multinomial logit with soft targets and two
## covariates, which is a handful of Newton steps on at most thirty
## parameters.  Doing it here rather than calling a package keeps the whole
## estimator in base R and avoids the package's handling of a
## proportion-valued response, which silently returns zero coefficients on
## the sparse targets this EM produces.
##   eta_ik = alpha_k + beta_k * z_i   for k = 2..K,   eta_i1 = 0
##   p_i = softmax(eta_i)
## maximising sum_i w_i sum_k t_ik log p_ik.
mnlogit <- function(target, z, w = NULL, iter = 25, tol = 1e-8, ridge = 1e-8) {
  n <- nrow(target); K <- ncol(target)
  if (is.null(w)) w <- rep(1, n)
  Z <- cbind(1, z)                        # n x 2
  d <- ncol(Z); m <- K - 1L
  B <- matrix(0, m, d)                    # rows: classes 2..K
  loglik <- function(BB) {
    eta <- cbind(0, Z %*% t(BB))
    mx <- apply(eta, 1, max)
    e <- exp(eta - mx); P <- e / rowSums(e)
    list(ll = sum(w * rowSums(target * log(pmax(P, 1e-300)))), P = P)
  }
  cur <- loglik(B); ll_prev <- -Inf
  for (it in seq_len(iter)) {
    if (is.finite(ll_prev) && abs(cur$ll - ll_prev) <= tol * (abs(ll_prev) + 1)) break
    ll_prev <- cur$ll
    P <- cur$P
    R <- target[, -1, drop = FALSE] - P[, -1, drop = FALSE]
    gm <- numeric(m * d)
    for (k in seq_len(m)) gm[((k - 1) * d + 1):(k * d)] <- crossprod(Z, R[, k] * w)
    H <- matrix(0, m * d, m * d)
    for (k in seq_len(m)) for (l in seq_len(k)) {
      a <- w * P[, k + 1] * ((k == l) - P[, l + 1])
      blk <- crossprod(Z * a, Z)
      ik <- ((k - 1) * d + 1):(k * d); il <- ((l - 1) * d + 1):(l * d)
      H[ik, il] <- blk
      if (k != l) H[il, ik] <- t(blk)
    }
    diag(H) <- diag(H) + ridge
    step <- tryCatch(solve(H, gm), error = function(e) NULL)
    if (is.null(step) || anyNA(step)) break
    moved <- FALSE
    for (sc in c(1, .5, .25, .1, .02)) {
      Bn <- B + matrix(sc * step, m, d, byrow = TRUE)
      nx <- loglik(Bn)
      if (is.finite(nx$ll) && nx$ll > cur$ll) { B <- Bn; cur <- nx; moved <- TRUE; break }
    }
    if (!moved) break
  }
  list(B = B, loglik = cur$ll)
}

## The fitted logit evaluated at arbitrary points.
mnlogit_probs <- function(B, z) {
  eta <- cbind(0, cbind(1, z) %*% t(B))
  mx <- apply(eta, 1, max)
  e <- exp(eta - mx)
  e / rowSums(e)
}

fit_rounding_logit <- function(y, weights = NULL, levels = NULL, delta = 0.02,
                               unit = 1, iter = 15, tol = 1e-4, trace = FALSE,
                               max_n = 4000, maxit = 25) {
  if (is.null(weights)) weights <- rep(1, length(y))
  ## the covariate is log income, so non-positive reports (which rounding to a
  ## coarse level can create out of small amounts) are dropped
  ok <- is.finite(y) & y > 0 & is.finite(weights) & weights > 0
  y <- y[ok]; weights <- weights[ok]
  if (length(y) < 100) stop("too few positive reports to fit a mechanism")
  if (is.null(levels)) levels <- default_levels(y)
  ## The mechanism is a low-dimensional object; a large administrative file
  ## does not need all of its records to pin it down, and the multinomial
  ## logit of the M-step is the expensive part.
  y_all <- y; w_all <- weights
  if (length(y) > max_n) {
    i <- sample.int(length(y), max_n)
    y <- y[i]; weights <- weights[i]
  }
  levels <- sort(unique(levels))
  lv <- round(levels / unit)
  W <- Reduce(lcm2, lv)
  v <- round(y / unit)
  pat <- pattern_id(v, lv)
  up <- sort(unique(pat))
  laws <- sapply(lv, function(L) {
    pl <- pattern_law(L, lv, W)
    out <- setNames(numeric(length(up)), as.character(up))
    cm <- intersect(names(pl), names(out)); out[cm] <- pl[cm]; out
  })
  dimnames(laws) <- list(as.character(up), as.character(levels))
  A <- laws[match(pat, up), , drop = FALSE]        # n x J, the known law
  colnames(A) <- as.character(levels)
  lx <- log(y); lx <- (lx - mean(lx)) / sd(lx)
  P <- matrix(1 / length(levels), length(y), length(levels))
  ll_old <- -Inf
  for (it in seq_len(iter)) {
    num <- A * P
    den <- rowSums(num)
    keep <- den > 0
    ll <- sum(weights[keep] * log(den[keep]))
    post <- num[keep, , drop = FALSE] / den[keep]
    ## levels that carry no posterior mass would give the multinomial logit
    ## an empty class and NA coefficients; they are dropped and restored at
    ## probability zero afterwards
    alive <- colSums(post * weights[keep]) / sum(weights[keep]) > 1e-6
    if (sum(alive) < 2) break
    fit <- mnlogit(post[, alive, drop = FALSE], lx[keep], weights[keep],
                   iter = maxit)
    P <- matrix(0, length(lx), ncol(A), dimnames = list(NULL, colnames(A)))
    P[, alive] <- mnlogit_probs(fit$B, lx)
    if (trace) cat(sprintf("  it %2d  loglik %.3f\n", it, ll))
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * abs(ll_old)) break
    ll_old <- ll
  }
  ## evaluate on the whole sample, so that p_x lines up with the data
  mu <- mean(log(y)); sg <- sd(log(y))
  lxa <- (log(y_all) - mu) / sg
  P <- matrix(0, length(y_all), length(levels),
              dimnames = list(NULL, as.character(levels)))
  P[, alive] <- mnlogit_probs(fit$B, lxa)
  y <- y_all; weights <- w_all
  pbar <- colSums(P * weights) / sum(weights)
  ## fit diagnostic: the divisibility shares the fitted mechanism implies
  G <- outer(lv, lv, Vectorize(function(l, L) gcd2(L, l) / l))   # rows l, cols L
  s_fit <- as.numeric(G %*% pbar)
  s_obs <- heap_shares(y, levels, weights)
  gof <- max(abs(s_obs - s_fit))
  ne <- sum(weights)^2 / sum(weights^2)
  qb <- cut(y, unique(quantile(y, seq(0, 1, .2))), labels = FALSE,
            include.lowest = TRUE)
  pmaxband <- apply(sapply(sort(unique(qb)), function(k)
    colSums(P[qb == k, , drop = FALSE] * weights[qb == k]) / sum(weights[qb == k])),
    1, max)
  ## Which levels are active, and hence the grid, is decided by the marginal
  ## mixture of `fit_rounding`: that is the estimator the identification
  ## result is about, and the one used everywhere else in the paper.  The
  ## logit adds the dependence on the amount, which the jittering and
  ## re-rounding procedures need but the grid does not.  If it fits the
  ## divisibility shares no better than the marginal mixture, the marginal
  ## probabilities are returned instead.
  base <- tryCatch(fit_rounding(y, weights, levels, delta = delta, unit = unit),
                   error = function(e) NULL)
  tol_ne <- max(0.02, 3 / sqrt(ne))
  if (!is.null(base) && gof > tol_ne && base$gof < gof) {
    P <- matrix(base$p, length(y_all), length(levels), byrow = TRUE,
                dimnames = list(NULL, as.character(levels)))
    pbar <- base$p; s_fit <- base$shares_fit[as.character(levels)]
    gof <- base$gof
  }
  list(p = pbar, p_x = P, p_band_max = pmaxband, ll = ll_old, levels = levels,
       shares_obs = setNames(s_obs, levels), shares_fit = setNames(s_fit, levels),
       gof = gof, misfit = gof > tol_ne,
       active = if (!is.null(base)) base$active else levels[pbar >= delta & levels > unit],
       delta = delta)
}
