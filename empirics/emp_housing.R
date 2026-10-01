suppressMessages({ library(data.table); source("../code/amlib.R") })
LEV <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000)
d <- fread("../data/emp/acs22_housing_ca.csv.gz")
tasks <- list(list(var = "owner cost (SMOCP)", y = d$SMOCP, alloc = d$FSMOCP == 1),
              list(var = "contract rent (RNTP)", y = d$RNTP, alloc = d$FRNTP == 1),
              list(var = "gross rent (GRNTP)", y = d$GRNTP, alloc = d$FGRNTP == 1))
res <- rbindlist(lapply(tasks, function(tk) rbindlist(lapply(c("all", "reported"), function(sub) {
  ok <- is.finite(tk$y) & tk$y > 0 & d$WGTP > 0; if (sub == "reported") ok <- ok & !tk$alloc
  y <- round(tk$y[ok]); w <- d$WGTP[ok]
  m <- modality_test(y, w, levels = LEV, B = 999)
  data.table(var = tk$var, sample = sub, n = length(y), active = m$active, w0 = m$w0, share_below_W = round(sum(w[y < m$w0]) / sum(w), 3), spec_df = m$spec_df, spec_dist = signif(m$spec_dist, 3),
             spec_T_df = round(m$spec_T_df, 2), spec_p = m$spec_p, status = m$status, K = m$K, p = m$p_value,
             delta = if (is.null(m$delta)) NA else m$delta, gap_lo = if (is.null(m$gap_lo)) NA else m$gap_lo,
             gap_hi = if (is.null(m$gap_hi)) NA else m$gap_hi, pstar_lo = if (is.null(m$pstar_lo)) NA else m$pstar_lo,
             pstar_hi = if (is.null(m$pstar_hi)) NA else m$pstar_hi)
}))))
print(res); fwrite(res, "emp_housing.csv")
