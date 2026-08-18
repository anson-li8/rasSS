#' Simulate genotypes from an LD matrix
#'
#' @param n Number of individuals.
#' @param R LD correlation matrix.
#' @return Integer genotype matrix with values 0, 1, or 2.
#' @export
sim_genotypes <- function(n, R) {
  X <- mvtnorm::rmvnorm(n, sigma = R)
  X <- round(X)
  pmin(pmax(X, 0), 2)
}

#' Greedy LD pruning
#'
#' @param R LD correlation matrix.
#' @param thresh r-squared pruning threshold.
#' @return Logical vector indicating retained SNPs.
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

#' Marginal statistics for continuous trait
#'
#' @param X Genotype matrix.
#' @param y Phenotype vector.
#' @return List with `beta` and `z`.
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

#' Rao score marginal statistics for binary trait
#'
#' @param X Genotype matrix.
#' @param y Binary phenotype vector.
#' @return List with `beta` and `z`.
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

#' Simulate phenotype
#'
#' @param X Genotype matrix.
#' @param beta True effect-size vector.
#' @param trait Either `"continuous"` or `"binary"`.
#' @param seed Random seed.
#' @return Phenotype vector.
#' @export
make_pheno <- function(X, beta, trait, seed) {
  set.seed(seed)
  L <- as.numeric(X %*% beta) + stats::rnorm(nrow(X), 0, 3)
  if (trait == "continuous") L else as.numeric(L > stats::qnorm(0.8))
}