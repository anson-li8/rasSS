#' Summary-Statistic RAS Profile Generation
#'
#' Generates a Regional Association Score (RAS) profile using only 
#' summary statistics. For each pivotal SNP, it runs an adaptive window 
#' scan, calculates the burden statistic \code{\link{t_burden}}, 
#' and keeps the minimum p-value (highest RAS) across the grid.
#'
#' @param b_disc Numeric vector. Genome-wide discovery marginal effect sizes.
#' @param z_targ Numeric vector. Genome-wide target marginal Z-scores.
#' @param R Numeric matrix. Genome-wide LD correlation matrix (e.g., from an external reference panel).
#' @param mask Logical vector. Indicates which SNPs are retained after LD pruning.
#' @param skip1 Integer. Step size for pivotal SNPs (e.g., evaluate every 10th SNP).
#' @param skip2 Integer. Step size for the adaptive window scan at each SNP.
#' @param min_window_size Integer. Minimum half-window size (in SNPs).
#' @param max_window_size Integer. Maximum half-window size (in SNPs).
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item `x`: Integer vector of the genomic indices of the pivotal SNPs evaluated.
#'     \item `y`: Numeric vector of the \eqn{-\log_{10}(p)}-values from the most 
#'       significant adaptive window at each pivotal SNP.
#'   }
#' @export
scan_ss <- function(b_disc, z_targ, R, mask, skip1 = 10, skip2 = 3, 
                    min_window_size = 3, max_window_size = 30) {
  stopifnot(length(b_disc) == length(z_targ), length(mask) == length(b_disc))
  n_snps <- length(b_disc)
  sites <- seq(1, n_snps, by = skip1)
  sub_windows <- c(0, seq(min_window_size, max_window_size, by = skip2))
  if (sub_windows[length(sub_windows)] != max_window_size) {
    sub_windows <- c(sub_windows, max_window_size)
  }
  
  y_profile <- vapply(sites, function(s) {
    best_p <- 1
    for (ws in sub_windows) {
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws)
      win_snps <- win_snps[mask[win_snps]]
      if (length(win_snps) < 1) next
      
      t_val <- t_burden(b_disc[win_snps], z_targ[win_snps], R[win_snps, win_snps, drop = FALSE])
      best_p <- min(best_p, 2 * stats::pnorm(-abs(t_val)))
    }
    -log10(max(best_p, .Machine$double.xmin))
  }, numeric(1))
  
  list(x = sites, y = y_profile)
}