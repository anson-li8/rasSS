#' Simulate genotypes from an LD matrix
#'
#' Using a multivariate normal distribution constructed by the provided LD correlation matrix,
#' this function simulates genotypes in a population, then rounds and clips to
#' integer genotypic dosages (0, 1, 2).
#'
#' @param n Integer. Number of individuals to simulate.
#' @param R Numeric matrix. LD correlation matrix.
#' @return An `n` x `ncol(R)` integer matrix of simulated genotype dosages (values 0, 1, or 2).
#' @export
sim_genotypes <- function(n, R) {
  X <- mvtnorm::rmvnorm(n, sigma = R)
  X <- round(X)
  pmin(pmax(X, 0), 2)
}

#' Greedy LD pruning
#'
#' Removes redundant SNPs based on the criterion of squared correlation (\eqn{r^2})
#' exceeding a specified threshold. Traverses SNPs in order and drops any
#' following SNPs that are redundant with a retained one.
#'
#' @param R Numeric matrix. LD correlation matrix.
#' @param thresh Numeric. Squared correlation threshold (e.g., 0.2). SNPs with
#'   \eqn{r^2 > thresh} compared to a retained SNP are eliminated.
#' @return Logical vector of length `nrow(R)`; retained SNPs are indicated by `TRUE`.
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

#' Marginal GWAS statistics for continuous trait
#' 
#' Calculates per-SNP marginal effect sizes and Z-scores via ordinary least squares (OLS).
#'
#' @param X Numeric matrix. Genotypic dosage matrix (individuals x SNPs).
#' @param y Numeric vector. Continuous phenotype vector.
#' 
#' @return List with two elements:
#'   \itemize{
#'     \item `beta`: Marginal effect sizes.
#'     \item `z`: Marginal Z-scores (non-finite values set to 0).
#'   }
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

#' Marginal GWAS statistics for binary trait
#'
#' Calculates per-SNP Rao score test statistics for case/control phenotypes.
#'
#' @param X Numeric matrix. Genotype dosage matrix (individuals x SNPs).
#' @param y Numeric vector. Binary phenotype vector (0/1).
#' 
#' @return List with two elements:
#'   \itemize{
#'     \item `beta`: Score-based effect sizes (\eqn{U/V}).
#'     \item `z`: Rao score Z-scores (non-finite values set to 0).
#'   }
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

#' Simulate phenotype from causal effects
#' 
#' Generates phenotype based on linear combination of genotypes and true effect,
#' with added Gaussian noise.
#' 
#' @param X Numeric matrix. Genotype dosage matrix (individuals x SNPs).
#' @param beta Numeric vector. True effect sizes (0 for non-causal SNPs).
#' @param trait Character. Either `"continuous"` or `"binary"`.
#' @param seed Integer. Random seed for reproducibility.
#' 
#' @return Phenotype vector. (numeric liability for continuous, 0/1 status for binary).
#' @export
make_pheno <- function(X, beta, trait, seed) {
  set.seed(seed)
  L <- as.numeric(X %*% beta) + stats::rnorm(nrow(X), 0, 3)
  if (trait == "continuous") L else as.numeric(L > stats::qnorm(0.8))
}