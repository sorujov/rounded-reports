## Earnings and income in the ACS 2022 and CPS ASEC 2023: the whole procedure.
suppressMessages({ library(data.table); library(multimode); library(diptest)
  source("../code/rounding.R"); source("../code/heaping.R"); source("../code/benchmarks.R"); source("../code/amlib.R") })
LEV <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
D <- "../data/emp/"
tasks <- list()
for (st in c("ca", "ms", "pr")) { d <- fread(paste0(D, "acs22_person_", st, ".csv.gz"))
  tasks[[length(tasks) + 1]] <- list(src = paste("ACS", toupper(st)), var = "wage", y = d$WAGP, w = d$PWGTP, alloc = d$FWAGP == 1)
  tasks[[length(tasks) + 1]] <- list(src = paste("ACS", toupper(st)), var = "personal income", y = d$PINCP, w = d$PWGTP, alloc = d$FPINCP == 1) }
d <- fread(paste0(D, "cps23_person.csv.gz"))
tasks[[length(tasks) + 1]] <- list(src = "CPS 2023", var = "wage", y = d$WSAL_VAL, w = d$MARSUPWT, alloc = d$I_WSVAL != 0 | d$I_ERNVAL != 0)
B <- 499
res <- rbindlist(lapply(tasks, function(tk) rbindlist(lapply(c("all", "reported"), function(sub) {
  ok <- is.finite(tk$y) & tk$y > 0 & tk$w > 0; if (sub == "reported") ok <- ok & !tk$alloc
  y <- round(tk$y[ok]); w <- tk$w[ok]; set.seed(1)
  m <- modality_test(y, w, levels = LEV, B = B)
  idx <- sample(length(y), min(6000, length(y)), TRUE, prob = w); ys <- y[idx]
  p_acr <- tryCatch(suppressWarnings(modetest(log(ys), method = "ACR", B = B))$p.value, error = function(e) NA)
  p_ch <- tryCatch(suppressWarnings(modetest(log(ys), method = "CH", B = B))$p.value, error = function(e) NA)
  data.table(source = tk$src, var = tk$var, sample = sub, n = length(y), median = round(median(rep(y, 1))),
             wmedian = y[order(y)][which(cumsum(w[order(y)]) >= sum(w) / 2)[1]],
             m1000 = round(sum(w[y %% 1000 == 0]) / sum(w), 3), m5000 = round(sum(w[y %% 5000 == 0]) / sum(w), 3),
             m10000 = round(sum(w[y %% 10000 == 0]) / sum(w), 3), m25000 = round(sum(w[y %% 25000 == 0]) / sum(w), 3),
             active = m$active, w0 = m$w0, share_below_W = round(sum(w[y < m$w0]) / sum(w), 3),
             spec_T_df = round(m$spec_T_df, 1), spec_df = m$spec_df, spec_p = m$spec_p, spec_dist = signif(m$spec_dist, 3),
             status = m$status, K = m$K, p_cell = m$p_value, p_acr = p_acr, p_ch = p_ch)
}))))
print(res); fwrite(res, "emp_earn.csv")
