## Latent income distributions and a reporting-heaping mechanism calibrated
## to ACS 2022 wages (84% multiples of 1000, 43% of 5000, 28% of 10000).

## ---- latent distributions (all positive, income-like) -------------------
LATENT <- list(
  LN   = list(type = "null", label = "Lognormal", r = function(n) rlnorm(n, 10.6, 0.9)),
  GA   = list(type = "null", label = "Gamma",     r = function(n) rgamma(n, 2, scale = 25000)),
  SM   = list(type = "null", label = "Singh--Maddala",
              r = function(n) { u <- runif(n); 45000 * ((1 - u)^(-1 / 3) - 1)^(1 / 2.8) }),
  SHLD = list(type = "null", label = "Shoulder",
              r = function(n) { z <- runif(n) < .75
                ifelse(z, rlnorm(n, 10.4, .55), rlnorm(n, 11.25, .30)) })
)

## Bimodal alternatives: mixture of two lognormals; `sep` controls the
## separation and therefore the width of the valley relative to the grid.
bimodal <- function(mu1, mu2, s1 = 0.35, s2 = 0.35, pi1 = 0.5)
  function(n) { z <- runif(n) < pi1; ifelse(z, rlnorm(n, mu1, s1), rlnorm(n, mu2, s2)) }

## Valley widths are chosen relative to the coarsest heaping grid
## (w = 10000 in the calibrated mechanism): N1--N3 have valleys narrower
## than or comparable to the grid and are therefore at or beyond the limit
## of what the reported data can identify.
ALT <- list(
  N1 = list(label = "valley 6k (below grid)",  x_star = 52932, p_star = 0.501, width =  5985,
            r = bimodal(log(50000), log(56000), .03, .03)),
  N2 = list(label = "valley 10k (at grid)",    x_star = 52814, p_star = 0.501, width =  9959,
            r = bimodal(log(48000), log(58000), .05, .05)),
  N3 = list(label = "valley 16k",              x_star = 53477, p_star = 0.501, width = 15918,
            r = bimodal(log(46000), log(62000), .07, .07)),
  M1 = list(label = "valley 25k",              x_star = 56464, p_star = 0.504, width = 24573,
            r = bimodal(log(45000), log(70000), .12, .12)),
  M2 = list(label = "valley 29k",              x_star = 58792, p_star = 0.507, width = 29091,
            r = bimodal(log(45000), log(75000), .15, .15)),
  M3 = list(label = "valley 52k",              x_star = 56886, p_star = 0.502, width = 52380,
            r = bimodal(log(35000), log(90000), .22, .22)),
  M4 = list(label = "valley 52k, weights 30/70", x_star = 53195, p_star = 0.298, width = 51573,
            r = bimodal(log(35000), log(90000), .25, .25, .3)),
  M5 = list(label = "low-income cluster",      x_star = 23592, p_star = 0.365, width = 33929,
            r = bimodal(log(12000), log(55000), .30, .45, .35)),
  M6 = list(label = "valley 108k",             x_star = 73120, p_star = 0.498, width = 107715,
            r = bimodal(log(30000), log(150000), .45, .35))
)

## ---- heaping mechanism ---------------------------------------------------
## Rounding level depends on magnitude: larger reports are rounded more
## coarsely.  Probabilities are chosen so that the marginal shares match ACS.
heap <- function(x, strength = 1) {
  n <- length(x)
  lv <- c(1, 100, 1000, 5000, 10000)
  big <- pmin(pmax((log(x) - log(20000)) / (log(200000) - log(20000)), 0), 1)
  P <- cbind(0.08 + 0.10 * (1 - big), 0.08 + 0.08 * (1 - big),
             0.40 - 0.05 * big,       0.20 + 0.05 * big,
             0.24 + 0.02 * big)
  if (strength != 1) {                      # shrink towards exact reporting
    P <- P * strength; P[, 1] <- P[, 1] + (1 - strength)
  }
  P <- P / rowSums(P)
  u <- runif(n); cp <- t(apply(P, 1, cumsum))
  k <- max.col(u < cp, ties.method = "first")
  l <- lv[k]
  round(x / l) * l
}

## ---- additional latent shapes for the difficult cases --------------------
## PLAT is the design the projection bootstrap should find hardest: a genuine
## flat stretch, so that many cells sit exactly on the boundary of the
## unimodal cone.  The flat component is a uniform smoothed by a normal, so
## that the density is continuous and has no spurious turning point at the
## edges of the plateau; over 50--95 thousand its coefficient of variation
## is under 0.1.  S1 and S2 are shallow valleys between the same two modes
## as M3, so that power can be reported against depth and not only width.
LATENT$PLAT <- list(type = "null", label = "Plateau",
  r = function(n) { u <- runif(n)
    ifelse(u < .45, rlnorm(n, 10.35, .40),
      ifelse(u < .80, runif(n, 40000, 105000) + rnorm(n, 0, 7000),
             rlnorm(n, 11.35, .35))) })
shallow <- function(a) function(n) {
  z <- runif(n) < a
  ifelse(z, rlnorm(n, log(55000), .55),
         ifelse(runif(n) < .5, rlnorm(n, log(35000), .22), rlnorm(n, log(90000), .22)))
}
ALT$S1 <- list(label = "valley 52k, depth 0.85", x_star = 60469, p_star = 0.539,
               width = 52000, depth = 0.845, r = shallow(.45))
ALT$S2 <- list(label = "valley 52k, depth 0.66", x_star = 58539, p_star = 0.519,
               width = 52000, depth = 0.657, r = shallow(.30))

## ---- reporting mechanisms ------------------------------------------------
## `heap` is the magnitude-dependent mechanism calibrated on ACS wages.
## The variants below are used to check that nothing depends on the
## mechanism being the one the de-heaping model assumes.
##  heap_25k : the same plus genuine rounding to 25000, so that the grid
##             must widen to 50000.
##  heap_mis : respondent types (a third of respondents round coarsely
##             whatever their income, the rest almost never), and an
##             off-ladder level of 2500 -- neither feature is in the model
##             fitted by the de-heaping algorithm.
##  heap_hid : the ACS mechanism plus a small amount of rounding at `hidden`
##             (default 25000), used to calibrate the activity cut-off.
heap_25k <- function(x, q = 0.25) {
  y <- heap(x); k <- runif(length(x)) < q
  y[k] <- round(x[k] / 25000) * 25000; y
}
heap_mis <- function(x) {
  n <- length(x)
  coarse <- runif(n) < 1/3
  lv <- c(1, 100, 1000, 2500, 5000, 10000)
  P <- matrix(0, n, length(lv))
  P[coarse, ] <- matrix(c(.02, .03, .20, .15, .30, .30), sum(coarse), 6, byrow = TRUE)
  P[!coarse, ] <- matrix(c(.55, .20, .18, .04, .02, .01), sum(!coarse), 6, byrow = TRUE)
  u <- runif(n); cp <- t(apply(P, 1, cumsum))
  l <- lv[max.col(u < cp, ties.method = "first")]
  round(x / l) * l
}
heap_hid <- function(x, q = 0.03, hidden = 25000) {
  y <- heap(x); k <- runif(length(x)) < q
  y[k] <- round(x[k] / hidden) * hidden; y
}
MECH <- list(none = function(x) x, acs = heap, k25 = heap_25k, mis = heap_mis)
