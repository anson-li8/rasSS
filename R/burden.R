#' Summary-Statistic Weighted Burden Statistic
#'
#' Calculates the 1-degree-of-freedom weighted burden statistic 
#' \eqn{T_\text{burden} = w^T Z / \sqrt{w^T Rw}}{T_\text{burden} = w^T Z / \sqrt{w'Rw}} of a given window of SNPs. Under the null 
#' hypothesis of no association in the window, \eqn{T}{T} asymptotically follows a  
#' standard normal distribution \eqn{N(0,1)}{N(0,1)}.
#'
#' @param w Numeric vector. Discovery weights (e.g., marginal effect sizes 
#'   from the discovery cohort) for the SNPs in the window.
#' @param Z Numeric vector. Marginal Z-scores from the target cohort for 
#'   the SNPs in the window. Must be the same length as `w`.
#' @param R Numeric matrix. Linkage Disequilibrium (LD) correlation matrix 
#'   for the SNPs in the window. Dimensions must be `length(w)` x `length(w)`.
#'
#' @return A single numeric value representing the burden test statistic. 
#'   Returns `0` if the variance denominator is non-positive or non-finite.
#' @export
t_burden <- function(w, Z, R) {
  denom <- sqrt(as.numeric(t(w) %*% R %*% w))
  if (!is.finite(denom) || denom <= 0) return(0)
  as.numeric((t(w) %*% Z) / denom)
}