suppressMessages({ library(data.table); source("../code/amlib.R") })
LEV <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
d <- fread("../data/emp/acs22_person_ca.csv.gz"); ok <- is.finite(d$WAGP) & d$WAGP > 0 & d$PWGTP > 0 & d$FWAGP != 1
y <- round(d$WAGP[ok]); wt <- d$PWGTP[ok]
fm <- fit_mech(y, wt, LEV); print(round(fm$p[fm$p > 0.001], 4)); print(fm$active)
A_centred <- function(lv, pv, w, K) {   # cells centred at k*w, k=1..K; boundary reports split
  A <- matrix(0, K, K)
  for (j in seq_along(lv)) { L <- lv[j]; if (pv[j] <= 0) next
    for (k in seq_len(K)) { a <- k * w - w / 2; b <- k * w + w / 2
      for (m in floor(a / L + 0.5):ceiling(b / L + 0.5)) {
        ov <- max(0, min(b, m * L + L / 2) - max(a, m * L - L / 2)) / w; if (ov <= 0) next
        r <- m * L / w; if (abs(r - floor(r) - 0.5) < 1e-9) tg <- c(floor(r), ceiling(r)) else tg <- round(r)
        for (t in tg) { t <- min(max(t, 1), K); A[t, k] <- A[t, k] + pv[j] * ov / length(tg) } } } }
  A }
act <- fm$active; pv <- fm$p[as.character(act)]; pv <- pv / sum(pv)
top <- 250000
for (w in c(25000, 10000, 5000, 1000, 500, 100)) { K <- top / w
  A <- A_centred(act, pv, w, K); cat(w, K, signif(kappa(A, exact = TRUE), 3), "\n") }
