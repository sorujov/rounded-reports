## Exactly flat unimodal designs (piecewise-linear densities, sampled by rejection).
rpiece <- function(n, dens, lo, hi, M) { out <- numeric(0)
  while (length(out) < n) { x <- runif(2 * n, lo, hi); keep <- runif(2 * n) * M < dens(x); out <- c(out, x[keep]) }
  out[seq_len(n)] }
## flat on the decreasing side: rise 0->2c on [0,30k], fall 2c->c on [30k,40k], flat c on [40k,100k], exp tail
r_flatdec <- function(n) { m <- c(45, 60, 30) / 135; z <- sample(1:3, n, TRUE, m)
  d1 <- function(x) ifelse(x < 30000, 2 * x / 30000, 2 - (x - 30000) / 10000)
  x <- numeric(n); i <- z == 1; x[i] <- rpiece(sum(i), d1, 0, 40000, 2)
  i <- z == 2; x[i] <- runif(sum(i), 40000, 100000); i <- z == 3; x[i] <- 100000 + rexp(sum(i), 1 / 30000); x }
## flat top: rise 0->c on [0,40k], flat c on [40k,100k], exp tail
r_flattop <- function(n) { m <- c(20, 60, 30) / 110; z <- sample(1:3, n, TRUE, m); x <- numeric(n)
  i <- z == 1; x[i] <- 40000 * sqrt(runif(sum(i))); i <- z == 2; x[i] <- runif(sum(i), 40000, 100000)
  i <- z == 3; x[i] <- 100000 + rexp(sum(i), 1 / 30000); x }
