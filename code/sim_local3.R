## Inverted test with the true profile, on the correct ladder and on a ladder that omits the 5000 level.
source("rounding.R"); source("heaping.R"); source("dgps_heap.R"); source("deconv.R"); source("localprof.R")
lv <- c(1, 100, 1000, 5000, 10000); XMAX <- 400000; W5 <- 5000
pi_strong <- function(x) { b <- pmin(pmax((log(pmax(x, 1)) - log(15000)) / (log(150000) - log(15000)), 0), 1)
  P <- cbind(0.30 - 0.25 * b, 0.25 - 0.20 * b, 0.30 - 0.05 * b, 0.10 + 0.15 * b, 0.05 + 0.35 * b); P / rowSums(P) }
heap_strong <- function(x) { P <- pi_strong(x); u <- runif(length(x)); cp <- t(apply(P, 1, cumsum))
  l <- lv[max.col(u < cp, ties.method = "first")]; round(x / l) * l }
lv_om <- c(1, 100, 1000, 10000)
pi_om <- function(x) { P <- pi_strong(x)[, -4, drop = FALSE]; P / rowSums(P) }
DES <- c("LN", "GA", "PLAT", "SM", "SHLD"); NS <- c(50000, 200000); REPS <- 60
grid <- expand.grid(rep = seq_len(REPS), n = NS, design = DES, stringsAsFactors = FALSE)
start <- as.integer(Sys.getenv("UNIT_START", "0")); end <- as.integer(Sys.getenv("UNIT_END", as.character(nrow(grid))))
workers <- as.integer(Sys.getenv("WORKERS", "1")); out_dir <- Sys.getenv("RESULTS_DIR", "results3"); dir.create(out_dir, FALSE, TRUE)
run_unit <- function(u) {
  g <- grid[u + 1, ]; set.seed(40260927 + u)
  x <- LATENT[[g$design]]$r(g$n); x <- x[x > 0 & x < XMAX]; y <- heap_strong(x); K <- as.integer(XMAX / W5)
  f <- function(L, pf) tryCatch(deconv_test(y, W5, L, pf, B = 99, K = K)$p_value, error = function(e) NA_real_)
  data.frame(design = g$design, n = g$n, rep = g$rep, p_true = f(lv, pi_strong), p_omit = f(lv_om, pi_om))
}
safe_unit <- function(u) tryCatch(run_unit(u), error = function(e) { cat('ERR unit', u, conditionMessage(e), '\n'); NULL })
rows <- parallel::mclapply(start:(end - 1), safe_unit, mc.cores = workers, mc.preschedule = FALSE)
ok <- vapply(rows, is.data.frame, TRUE)
write.csv(do.call(rbind, rows[ok]), file.path(out_dir, "chunk_0.csv"), row.names = FALSE)
cat(sprintf("UNITS_DONE=%d\n", end - start))
