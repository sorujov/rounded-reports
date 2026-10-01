suppressMessages({ library(data.table); source("../code/amlib.R") })
LEV <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
for (yr in 2016:2020) { d <- fread(sprintf("../data/cps/onep_%d.csv", yr)); d <- d[dhi > 0 & hsup_wgt > 0 & htotval > 0]
  d[, filt := !(pubtrans > 0.5 * htotval)]; d[, er5 := earnings > 0 & earnings %% 5000 == 0]
  for (s in c("filtered", "filt_earn5000")) { i <- if (s == "filtered") which(d$filt) else which(d$filt & d$er5)
    for (ww in c(500, 1000, 2000)) { set.seed(yr); m <- modality_test(d$dhi[i], d$hsup_wgt[i], levels = LEV, w = ww, B = 999, spec = FALSE)
      cat(yr, s, ww, m$K, m$p_value, if (!is.null(m$gap_lo)) paste(m$gap_lo, m$gap_hi, round(m$pstar_lo, 3), round(m$pstar_hi, 3)) else "", "\n") } } }
