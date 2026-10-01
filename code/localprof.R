## Local pattern-EM estimate of the latent-conditional profile P(L = l | X in band).
local_profile <- function(y, levels, h, min_n = 1500, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(y))
  W <- Reduce(lcm2, levels)
  laws_all <- lapply(levels, function(L) pattern_law(L, levels, W))
  br <- seq(0, max(y) + h, by = h)
  b <- findInterval(y, br)
  ## merge sparse bands upwards/downwards
  cnt <- tabulate(b, nbins = length(br))
  grp <- integer(length(cnt)); g <- 1; acc <- 0
  for (i in seq_along(cnt)) { grp[i] <- g; acc <- acc + cnt[i]; if (acc >= min_n) { g <- g + 1; acc <- 0 } }
  if (acc > 0 && acc < min_n) grp[grp == g] <- max(1, g - 1)
  G <- grp[b]; ug <- sort(unique(G))
  mids <- numeric(length(ug)); P <- matrix(NA, length(ug), length(levels))
  for (j in seq_along(ug)) {
    i <- which(G == ug[j]); pats <- pattern_id(round(y[i]), levels)
    cnt2 <- tapply(weights[i], pats, sum)
    laws <- sapply(laws_all, function(pl) { out <- setNames(numeric(length(cnt2)), names(cnt2))
      cm <- intersect(names(pl), names(cnt2)); out[cm] <- pl[cm]; out })
    keep <- rowSums(laws) > 0
    P[j, ] <- em_mix(as.numeric(cnt2[keep]), laws[keep, , drop = FALSE], iter = 3000)
    mids[j] <- weighted.mean(y[i], weights[i])
  }
  okr <- which(rowSums(is.na(P)) == 0); mids <- mids[okr]; P <- P[okr, , drop = FALSE]
  if (length(mids) == 1) { P <- rbind(P, P); mids <- c(mids, mids + 1) }
  function(x) {
    out <- sapply(seq_along(levels), function(k) approx(mids, P[, k], xout = x, rule = 2)$y)
    out <- matrix(pmax(out, 0), nrow = length(x)); out / rowSums(out)
  }
}
