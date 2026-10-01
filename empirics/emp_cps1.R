suppressMessages({ library(data.table); source("../code/amlib.R") })
LEV <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
out <- list()
for (yr in 2016:2020) {
  d <- fread(sprintf("../data/cps/onep_%d.csv", yr)); d <- d[dhi > 0 & hsup_wgt > 0 & htotval > 0]
  d[, filt := !(pubtrans > 0.5 * htotval)]
  d[, er5 := earnings > 0 & earnings %% 5000 == 0]
  subs <- list(full = rep(TRUE, nrow(d)), filtered = d$filt, filt_noalloc = d$filt & d$any_alloc == 0,
               filt_earn5000 = d$filt & d$er5, filt_not5000 = d$filt & !d$er5)
  for (s in names(subs)) { i <- which(subs[[s]]); set.seed(yr)
    t0 <- Sys.time(); m <- modality_test(d$dhi[i], d$hsup_wgt[i], levels = LEV, B = 999); cat('mt', s, format(Sys.time() - t0), '\n')
    t0 <- Sys.time(); g <- fit_mech(d$htotval[i], d$hsup_wgt[i], LEV); cat('fm', format(Sys.time() - t0), '\n')
    out[[length(out) + 1]] <- data.table(year = yr, sample = s, n = length(i), active = m$active, w0 = m$w0, w = if (is.null(m$w)) NA else m$w, spec_df = m$spec_df, w = m$w,
      pi1 = round(m$p_hat[1], 3), spec_T_df = round(m$spec_T_df, 2), spec_p = m$spec_p, status = m$status, K = m$K,
      p = m$p_value, delta = if (is.null(m$delta)) NA else round(m$delta, 4),
      gap = if (is.null(m$gap_lo)) NA else paste0(m$gap_lo, "-", m$gap_hi),
      pstar = if (is.null(m$pstar_lo)) NA else paste0(round(m$pstar_lo, 3), "-", round(m$pstar_hi, 3)),
      gross_active = paste(g$active, collapse = "/"))
    print(out[[length(out)]]) } }
fwrite(rbindlist(out), "emp_cps1.csv")
