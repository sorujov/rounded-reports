suppressMessages(library(data.table))
f3 <- function(x) formatC(x, format = "f", digits = 3); f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) ifelse(is.na(p), "--", ifelse(p <= 0.002, "$<0.01$", f2(p)))
big <- function(x) formatC(x, format = "d", big.mark = "{,}")
e <- fread("emp_earn.csv")[sample == "all"]
lab <- c("wage" = "wages", "personal income" = "personal income")
cat("%% earnings table body\n")
for (i in seq_len(nrow(e))) with(e[i], cat(sprintf("%s & %s & %s & %s & %.0f & %.0f & %.0f & %s & %s & %.2f & %s & %s & %d & %s \\\\\n",
  source, lab[var], big(n), big(wmedian), 100*m1000, 100*m5000, 100*m10000, big(w0), f2(share_below_W), spec_T_df, big(spec_df), f3(spec_dist), K, f2(p_cell))))
h <- fread("emp_housing.csv")[sample == "all"]
cat("%% housing table body\n")
for (i in seq_len(nrow(h))) with(h[i], cat(sprintf("%s & %s & %s & %s & %.2f & %d & %s & %d & %s & %s & %s \\\\\n",
  var, big(n), gsub("/", ", ", active), big(w0), spec_T_df, spec_df, f3(spec_dist), K, f2(p),
  ifelse(is.na(gap_lo), "--", sprintf("%s--%s", big(gap_lo), big(gap_hi))),
  ifelse(is.na(pstar_lo), "--", sprintf("%s--%s", f2(pstar_lo), f2(pstar_hi))))))
c1 <- fread("emp_cps1.csv")
cat("%% cps table body\n")
for (s in c("full", "filtered", "filt_earn5000", "filt_not5000")) for (yr in 2016:2020) with(c1[sample == s & year == yr],
  cat(sprintf("%s & %d & %s & %s & %s & %.2f & %s & %d & %s & %s & %s \\\\\n", s, yr, big(n), ifelse(active == "", "none", gsub("/", ", ", active)), big(w),
  spec_T_df, fp(spec_p), K, fp(p), gsub("-", "--", gsub("(\\d)(\\d{3})(?!\\d)", "\\1{,}\\2", gap, perl = TRUE)), gsub("-", "--", pstar))))
