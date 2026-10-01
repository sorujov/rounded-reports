## =====================================================================
## The obvious repairs, so that they can be compared with the two
## procedures of the paper rather than left for a referee to propose.
##
##  (i)   smooth: excess mass on a kernel estimate whose bandwidth is not
##        allowed below a fraction of the coarsest lattice.  A Gaussian
##        kernel of bandwidth h attenuates a lattice of spacing l by
##        2 exp(-2 pi^2 h^2 / l^2), which is below 1e-4 once h >= 0.7 l,
##        so the heaps are smoothed away by construction.
##  (ii)  jitter: replace each report by a uniform draw over the rounding
##        interval it is consistent with, and test the jittered sample.
##        This is the interval-data reading of rounding.
##  (iii) re-rounding bootstrap: keep the usual statistic on the reported
##        values, but calibrate it against a null that has been put
##        through the same reporting mechanism -- resample from a unimodal
##        fit of the latent scale, round the draws with the estimated
##        mechanism, recompute the statistic.
## =====================================================================

## ---- (i) bandwidth floor -------------------------------------------------
## `h_min` is on the scale of `x`.  For income the heaps are equally spaced
## in levels and the density is usually analysed in logs, where the spacing
## l/x is not constant; the floor that matches a level-scale lattice l near
## the centre of the data is then 0.7 * l / median(x), and we report the test
## on both scales because neither is obviously the right one.
crit_bw_test <- function(x, h_min = 0, k = 1, B = 199, alpha = 0.05, n_grid = 512) {
  x <- as.numeric(x); n <- length(x); s <- sd(x)
  nmodes <- function(h) {
    d <- density(x, bw = h, n = n_grid)$y
    sum(diff(sign(diff(d))) == -2)
  }
  lo <- max(h_min, s / 200); hi <- 2 * s
  if (nmodes(lo) <= k) return(list(h_crit = lo, p_value = 1, floored = TRUE))
  for (i in 1:60) { mid <- sqrt(lo * hi); if (nmodes(mid) > k) lo <- mid else hi <- mid }
  hc <- hi
  lam <- (0.94029 * alpha^3 - 1.59914 * alpha^2 + 0.17695 * alpha + 0.48971) /
         (alpha^3 - 1.77793 * alpha^2 + 0.36162 * alpha + 0.42423)
  hb <- lam * hc
  cnt <- 0
  for (b in 1:B) {
    xb <- sample(x, n, TRUE)
    xb <- mean(xb) + (xb - mean(xb) + hb * rnorm(n)) / sqrt(1 + hb^2 / var(x))
    if (sum(diff(sign(diff(density(xb, bw = hb, n = n_grid)$y))) == -2) > k) cnt <- cnt + 1
  }
  list(h_crit = hc, p_value = (1 + cnt) / (B + 1), floored = FALSE)
}

## ---- (ii) jitter over the rounding interval ------------------------------
## Rounding read as interval data: the latent value is somewhere in the
## interval the report is consistent with.  Which interval that is depends on
## the level used, which is not observed, so we draw it from its posterior
## given the report's divisibility pattern and the fitted mechanism -- the
## fair version of the naive rule "jitter by the coarsest lattice that
## divides the report", which over-assigns coarse levels to every round
## number that happened to come out of a fine rounding.
jitter_reports <- function(y, fit, unit = 1) {
  levels <- if (!is.null(fit$levels)) fit$levels else as.numeric(names(fit$p))
  lv <- round(levels / unit)
  W <- Reduce(lcm2, lv)
  v <- round(y / unit)
  pat <- pattern_id(v, lv)
  up <- sort(unique(pat))
  laws <- sapply(lv, function(L) {
    pl <- pattern_law(L, lv, W)
    out <- setNames(numeric(length(up)), as.character(up)); cm <- intersect(names(pl), names(out))
    out[cm] <- pl[cm]; out
  })
  A <- laws[match(pat, up), , drop = FALSE]
  ## the prior over levels may depend on the amount; `p_x` from the logit fit
  ## is used when it is there, and the marginal otherwise
  P <- if (!is.null(fit$p_x)) fit$p_x else
    matrix(fit$p, length(y), length(levels), byrow = TRUE)
  post <- A * P
  den <- rowSums(post)
  bad <- den <= 0
  post[bad, ] <- 1 / length(levels); den[bad] <- 1
  post <- post / den
  u <- runif(length(y)); cp <- t(apply(post, 1, cumsum))
  L <- levels[max.col(u < cp, ties.method = "first")]
  y + runif(length(y), -L / 2, L / 2)
}

## ---- (iii) calibrate the usual statistic against a re-rounded null -------
## `stat_fun` is applied to the reported values; the null draws are made by
## resampling the unimodal projection on a fine latent grid and putting the
## draws through `round_fun`.
reround_test <- function(y, fit, B = 199, n_cells = 60, stat_fun = NULL) {
  n <- length(y)
  levels <- if (!is.null(fit$levels)) fit$levels else as.numeric(names(fit$p))
  ## the default statistic is the dip on log values, the statistic the
  ## applied literature uses; any statistic of the reported sample will do
  if (is.null(stat_fun))
    stat_fun <- function(z) as.numeric(diptest::dip(log(z[z > 0])))
  qs <- as.numeric(quantile(y, c(0.005, 0.995)))
  br <- seq(qs[1], qs[2], length.out = n_cells + 1)
  ctr <- (br[-1] + br[-length(br)]) / 2; hw <- diff(br)[1] / 2
  cnt <- as.numeric(table(cut(y, br)))
  f0 <- unimodal_fit(cnt / sum(cnt))$f            # unimodal latent fit
  ## the reporting probabilities at a given amount: from the logit fit when
  ## it is available, so that the null draws are rounded the way reports of
  ## that size actually are
  lev_at <- if (is.null(fit$p_x)) {
    function(z) sample(levels, length(z), TRUE, prob = fit$p)
  } else {
    o <- order(y); ys <- y[o]; Ps <- fit$p_x[o, , drop = FALSE]
    function(z) {
      j <- pmin(pmax(findInterval(z, ys), 1), length(ys))
      pr <- Ps[j, , drop = FALSE]
      u <- runif(length(z)); cp <- t(apply(pr, 1, cumsum))
      levels[max.col(u < cp, ties.method = "first")]
    }
  }
  draw <- function() {
    j <- sample.int(n_cells, n, TRUE, prob = f0)
    z <- ctr[j] + runif(n, -hw, hw)
    L <- lev_at(z)
    round(z / L) * L
  }
  obs <- stat_fun(y)
  st <- replicate(B, stat_fun(draw()))
  list(stat = obs, p_value = (1 + sum(st >= obs)) / (B + 1))
}
