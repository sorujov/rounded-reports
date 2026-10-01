## Main Monte Carlo: size and power of every procedure, by reporting mechanism.
## Unit = (design, mechanism, n, replication).  Seeded by unit.
suppressMessages({ source("rounding.R"); source("heaping.R"); source("benchmarks.R"); source("dgps_heap.R"); source("flatdgp.R"); source("amlib.R") })
library(multimode); library(diptest)
LATENT$FLATDEC <- list(type = "null", label = "Flat, decreasing side", r = r_flatdec)
LATENT$FLATTOP <- list(type = "null", label = "Flat top", r = r_flattop)
env <- function(k, d) Sys.getenv(k, d)
NULLS <- names(LATENT); ALTS <- names(ALT); MECHS <- c("none", "acs", "k25", "mis"); NS <- c(2000, 10000)
RN <- as.integer(env("SIM_RN", "400")); RA <- as.integer(env("SIM_RA", "200")); B <- as.integer(env("SIM_B", "299"))
units <- rbind(expand.grid(rep = seq_len(RN) - 1L, n = NS, mech = MECHS, design = NULLS, stringsAsFactors = FALSE),
               expand.grid(rep = seq_len(RA) - 1L, n = NS, mech = MECHS, design = ALTS, stringsAsFactors = FALSE))
set.seed(20260927); units <- units[sample.int(nrow(units)), ]
N_UNITS <- nrow(units)
start <- as.integer(env("UNIT_START", "0")); end <- as.integer(env("UNIT_END", as.character(N_UNITS)))
workers <- as.integer(env("WORKERS", "1")); out_dir <- env("RESULTS_DIR", "results"); dir.create(out_dir, FALSE, TRUE)
LEV <- c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000)
NAr <- function(e) NA_real_
run_unit <- function(u) {
  un <- units[u + 1, ]; set.seed(700000L + u)
  gen <- if (un$design %in% NULLS) LATENT[[un$design]]$r else ALT[[un$design]]$r
  y <- MECH[[un$mech]](gen(un$n)); y <- round(y[y > 0]); t0 <- proc.time()[[3]]
  yl <- log(y)
  p_acr <- tryCatch(suppressWarnings(modetest(yl, method = "ACR", B = B))$p.value, error = NAr)
  p_ch  <- tryCatch(suppressWarnings(modetest(yl, method = "CH",  B = B))$p.value, error = NAr)
  fit <- tryCatch(fit_rounding(y, levels = LEV), error = function(e) NULL)
  fm <- fit_mech(y, levels = LEV); lmax <- if (length(fm$active)) max(fm$active) else 0
  p_fl <- tryCatch(crit_bw_test(yl, h_min = 0.7 * lmax / median(y), B = B)$p_value, error = NAr)
  p_jt <- if (is.null(fit)) NA_real_ else tryCatch({ yj <- jitter_reports(y, fit)
            suppressWarnings(modetest(log(yj[yj > 0]), method = "ACR", B = B))$p.value }, error = NAr)
  p_rr <- if (is.null(fit)) NA_real_ else tryCatch(reround_test(y, fit, B = B)$p_value, error = NAr)
  g <- tryCatch(modality_test(y, levels = LEV, B = 499, method = "gms", seed = u), error = function(e) list(status = "error"))
  gp <- if (identical(g$status, "ok")) tryCatch(cell_test(cell_split(y, g$w)$p, cell_split(y, g$w)$n_eff, B = 499, method = "proj", seed = u)$p_value, error = NAr) else NA_real_
  data.frame(unit = u, design = un$design, mech = un$mech, n = un$n, rep = un$rep,
             type = if (un$design %in% NULLS) "null" else "alt",
             p_acr = p_acr, p_ch = p_ch, p_floor = p_fl, p_jitter = p_jt, p_reround = p_rr,
             p_cell = if (identical(g$status, "ok")) g$p_value else NA_real_, p_cell_proj = gp,
             status = g$status, w0 = if (is.null(g$w0)) NA else g$w0, K = if (is.null(g$K)) NA else g$K,
             active = if (is.null(g$active)) NA else g$active,
             spec_T_df = if (is.null(g$spec_T_df)) NA else g$spec_T_df, spec_p = if (is.null(g$spec_p)) NA else g$spec_p,
             share_below_W = if (is.null(g$share_below_W)) NA else g$share_below_W, spec_dist = if (is.null(g$spec_dist)) NA else g$spec_dist,
             pi25000 = if ("25000" %in% names(fm$p)) unname(fm$p["25000"]) else NA,
             pstar_lo = if (is.null(g$pstar_lo)) NA else g$pstar_lo, pstar_hi = if (is.null(g$pstar_hi)) NA else g$pstar_hi,
             secs = proc.time()[[3]] - t0)
}
safe <- function(u) tryCatch(run_unit(u), error = function(e) data.frame(unit = u, error = conditionMessage(e)))
rows <- parallel::mclapply(start:(end - 1), safe, mc.cores = workers, mc.preschedule = FALSE)
cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(z) { z[setdiff(cols, names(z))] <- NA; z[cols] }))
write.csv(res, file.path(out_dir, sprintf("chunk_%s.csv", env("CHUNK_ID", "0"))), row.names = FALSE)
cat(sprintf("UNITS_DONE=%d\n", end - start))
