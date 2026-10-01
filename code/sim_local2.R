## Inverted test under a mechanism whose profile varies strongly with income.
source("rounding.R"); source("heaping.R"); source("dgps_heap.R"); source("deconv.R"); source("localprof.R")
lv <- c(1, 100, 1000, 5000, 10000); XMAX <- 400000; W5 <- 5000
pi_strong <- function(x) { b <- pmin(pmax((log(pmax(x, 1)) - log(15000)) / (log(150000) - log(15000)), 0), 1)
  P <- cbind(0.30 - 0.25 * b, 0.25 - 0.20 * b, 0.30 - 0.05 * b, 0.10 + 0.15 * b, 0.05 + 0.35 * b); P / rowSums(P) }
heap_strong <- function(x) { P <- pi_strong(x); u <- runif(length(x)); cp <- t(apply(P, 1, cumsum))
  l <- lv[max.col(u < cp, ties.method = "first")]; round(x / l) * l }
DES <- c("LN", "GA", "PLAT", "M1"); NS <- c(50000, 200000); REPS <- 60
grid <- expand.grid(rep = seq_len(REPS), n = NS, design = DES, stringsAsFactors = FALSE)
N_UNITS <- nrow(grid)
start <- as.integer(Sys.getenv("UNIT_START", "0")); end <- as.integer(Sys.getenv("UNIT_END", as.character(N_UNITS)))
workers <- as.integer(Sys.getenv("WORKERS", "1")); out_dir <- Sys.getenv("RESULTS_DIR", "results2"); dir.create(out_dir, FALSE, TRUE)
gen <- function(d) if (d %in% names(LATENT)) LATENT[[d]]$r else ALT[[d]]$r
run_unit <- function(u) {
  g <- grid[u + 1, ]; set.seed(30260927 + u)
  x <- gen(g$design)(g$n); x <- x[x > 0 & x < XMAX]; y <- heap_strong(x); K <- as.integer(XMAX / W5)
  pl8 <- local_profile(y, lv, 80000); pl16 <- local_profile(y, lv, 160000)
  fl <- tryCatch(fit_rounding_logit(y, levels = lv), error = function(e) NULL); plog <- if (!is.null(fl)) pi_from_fit(fl, y) else NULL
  fm <- fit_rounding(y, levels = lv); pm <- function(x) matrix(fm$p, length(x), length(lv), byrow = TRUE)
  f <- function(pf) if (is.null(pf)) NA_real_ else tryCatch(deconv_test(y, W5, lv, pf, B = 99, K = K)$p_value, error = function(e) NA_real_)
  xs <- c(12e3, 3e4, 6e4, 12e4)
  data.frame(design = g$design, n = g$n, rep = g$rep, p_true = f(pi_strong), p_local8 = f(pl8), p_local16 = f(pl16),
             p_logit = f(plog), p_marg = f(pm), err8 = max(abs(pl8(xs) - pi_strong(xs))), err16 = max(abs(pl16(xs) - pi_strong(xs))),
             errm = max(abs(pm(xs) - pi_strong(xs))))
}
safe_unit <- function(u) tryCatch(run_unit(u), error = function(e) { cat('ERR unit', u, conditionMessage(e), '\n'); NULL })
rows <- parallel::mclapply(start:(end - 1), safe_unit, mc.cores = workers, mc.preschedule = FALSE)
ok <- vapply(rows, is.data.frame, TRUE)
write.csv(do.call(rbind, rows[ok]), file.path(out_dir, sprintf("chunk_%s.csv", Sys.getenv("CHUNK_ID", "0"))), row.names = FALSE)
cat(sprintf("UNITS_DONE=%d\n", end - start))
