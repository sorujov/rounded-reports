## Auxiliary Monte Carlo.
## Part "spec": size and power of the residue specification test (Theorem 1):
##   correctly specified ladder, and the same ladder with its coarsest active level omitted.
## Part "delta": size of the cell test at the certified width when a level of probability q is ignored.
suppressMessages({ source("rounding.R"); source("heaping.R"); source("dgps_heap.R"); source("flatdgp.R"); source("amlib.R") })
LATENT$FLATDEC <- list(type = "null", label = "Flat, decreasing side", r = r_flatdec)
LATENT$FLATTOP <- list(type = "null", label = "Flat top", r = r_flattop)
env <- function(k, d) Sys.getenv(k, d)
NULLS <- names(LATENT); R <- as.integer(env("SIM_R", "200"))
LEV <- c(1, 10, 100, 500, 1000, 2500, 5000, 10000, 25000)
uS <- expand.grid(rep = seq_len(R) - 1L, n = c(2000, 5000, 20000), q = NA, mech = c("acs", "mis", "k25"),
                  design = NULLS, part = "spec", stringsAsFactors = FALSE)
uD <- expand.grid(rep = seq_len(R) - 1L, n = c(2000, 10000, 50000), q = c(0, .01, .02, .03, .05, .10), mech = "hid",
                  design = NULLS, part = "delta", stringsAsFactors = FALSE)
units <- rbind(uS, uD); set.seed(20260928); units <- units[sample.int(nrow(units)), ]
start <- as.integer(env("UNIT_START", "0")); end <- min(nrow(units), as.integer(env("UNIT_END", as.character(nrow(units)))))
workers <- as.integer(env("WORKERS", "1")); out_dir <- env("RESULTS_DIR", "results_aux"); dir.create(out_dir, FALSE, TRUE)
run_unit <- function(u) {
  un <- units[u + 1, ]; set.seed(900000L + u); t0 <- proc.time()[[3]]
  x <- LATENT[[un$design]]$r(un$n)
  if (un$part == "spec") {
    y <- MECH[[un$mech]](x); y <- round(y[y > 0])
    fm <- fit_mech(y, levels = LEV); act <- fm$active
    lv <- sort(unique(c(1, LEV[fm$w0 %% LEV == 0])))            # candidates dividing the certified width
    s1 <- if (length(lv) > 1) spec_test(y, lv, B = 299, seed = u) else NULL
    ## misspecified ladder: the analyst uses decimal levels only (1, 10, 100, ...), omitting 500, 2500, 5000, 25000
    lv2 <- lv[lv %in% 10^(0:6)]; om <- paste(setdiff(lv, lv2), collapse = "/")
    s2 <- if (nchar(om) > 0 && length(lv2) > 1) spec_test(y, lv2, B = 299, seed = u) else NULL
    data.frame(unit = u, part = un$part, design = un$design, mech = un$mech, n = un$n, rep = un$rep, q = NA,
               active = paste(act, collapse = "/"), w0 = fm$w0, share_below_W = if (is.null(s1)) NA else s1$share_below_W,
               spec_p = if (is.null(s1)) NA else s1$p, spec_T_df = if (is.null(s1)) NA else s1$T_df,
               omit = om, spec_p_omit = if (is.null(s2)) NA else s2$p,
               p_cell = NA, K = NA, secs = proc.time()[[3]] - t0)
  } else {
    y <- heap_hid(x, q = un$q, hidden = 25000); y <- round(y[y > 0])
    cp <- cell_split(y, 10000); t <- cell_test(cp$p, cp$n_eff, B = 499, method = "gms", seed = u)
    data.frame(unit = u, part = un$part, design = un$design, mech = un$mech, n = un$n, rep = un$rep, q = un$q,
               active = NA, w0 = 10000, share_below_W = NA, spec_p = NA, spec_T_df = NA, omit = NA, spec_p_omit = NA,
               p_cell = t$p_value, K = t$K, secs = proc.time()[[3]] - t0)
  }
}
safe <- function(u) tryCatch(run_unit(u), error = function(e) data.frame(unit = u, error = conditionMessage(e)))
rows <- parallel::mclapply(start:(end - 1), safe, mc.cores = workers, mc.preschedule = FALSE)
cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(z) { z[setdiff(cols, names(z))] <- NA; z[cols] }))
write.csv(res, file.path(out_dir, sprintf("chunk_%s.csv", env("CHUNK_ID", "0"))), row.names = FALSE)
cat(sprintf("UNITS_DONE=%d of %d\n", end - start, nrow(units)))
