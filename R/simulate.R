#' @export
sim_genotypes <- function(n, R) {
  X <- mvtnorm::rmvnorm(n, sigma = R)
  X <- round(X)
  pmin(pmax(X, 0), 2)
}

#' @export
prune_ld <- function(R, thresh) {
  keep <- rep(TRUE, nrow(R))
  for (i in seq_len(nrow(R) - 1)) {
    if (!keep[i]) next
    for (j in (i + 1):nrow(R)) {
      if (keep[j] && R[i, j]^2 > thresh) keep[j] <- FALSE
    }
  }
  keep
}

#' @export
get_marginal_stats <- function(X, y) {
  yc <- y - mean(y)
  Xc <- sweep(X, 2, colMeans(X), "-")
  sxx <- colSums(Xc^2)
  sxy <- as.numeric(crossprod(Xc, yc))
  beta <- sxy / sxx
  syy <- sum(yc^2)
  rss <- pmax(syy - beta * sxy, 0)
  se <- sqrt((rss / (length(y) - 2)) / sxx)
  z <- beta / se
  z[!is.finite(z)] <- 0
  list(beta = beta, z = z)
}

#' @export
get_marginal_stats_bin <- function(X, y) {
  p0 <- mean(y)
  Xc <- sweep(X, 2, colMeans(X), "-")
  U <- as.numeric(crossprod(Xc, y - p0))
  V <- p0 * (1 - p0) * colSums(Xc^2)
  z <- U / sqrt(pmax(V, 1e-12))
  z[!is.finite(z)] <- 0
  list(beta = U / pmax(V, 1e-12), z = z)
}

#' @export
make_pheno <- function(X, beta, trait, seed) {
  set.seed(seed)
  L <- as.numeric(X %*% beta) + stats::rnorm(nrow(X), 0, 3)
  if (trait == "continuous") L else as.numeric(L > stats::qnorm(0.8))
}