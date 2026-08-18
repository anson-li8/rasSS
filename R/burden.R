#' Summary-statistic weighted burden statistic
#'
#' @param w Weight vector.
#' @param Z Marginal Z-score vector.
#' @param R LD correlation matrix.
#' @return Numeric burden statistic.
#' @export
#' @examples
#' w <- c(1, 1); Z <- c(1, 1); R <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' t_burden(w, Z, R)
t_burden <- function(w, Z, R) {
  denom <- sqrt(as.numeric(t(w) %*% R %*% w))
  if (!is.finite(denom) || denom <= 0) return(0)
  as.numeric((t(w) %*% Z) / denom)
}