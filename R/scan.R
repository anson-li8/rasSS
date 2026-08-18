#' Summary-statistic RAS scan
#'
#' @param b_disc Discovery marginal effect-size weights.
#' @param z_targ Target marginal Z-scores.
#' @param R LD correlation matrix.
#' @param mask Logical SNP mask.
#' @param skip1 Pivotal SNP step size.
#' @param skip2 Window-size step size.
#' @param min_window_size Minimum half-window size.
#' @param max_window_size Maximum half-window size.
#' @return List with `x` and `y`.
#' @export
scan_ss <- function(b_disc, z_targ, R, mask, skip1 = 10, skip2 = 3, min_window_size = 3, max_window_size = 30) {
  n_snps <- length(b_disc)
  sites <- seq(1, n_snps, by = skip1)
  sub_windows <- c(0, seq(min_window_size, max_window_size, by = skip2))
  if (sub_windows[length(sub_windows)] != max_window_size) {
    sub_windows <- c(sub_windows, max_window_size)
  }
  
  y_profile <- sapply(sites, function(s) {
    best_p <- 1
    for (ws in sub_windows) {
      win_snps <- if (ws == 0) s else max(1, s - ws):min(n_snps, s + ws)
      win_snps <- win_snps[mask[win_snps]]
      if (length(win_snps) < 1) next
      
      w_sub <- b_disc[win_snps]
      z_sub <- z_targ[win_snps]
      R_sub <- R[win_snps, win_snps, drop = FALSE]
      
      num <- sum(w_sub * z_sub)
      denom <- sqrt(as.numeric(w_sub %*% R_sub %*% w_sub))
      if (!is.finite(denom) || denom <= 0) next
      
      t_val <- num / denom
      best_p <- min(best_p, 2 * stats::pnorm(-abs(t_val)))
    }
    -log10(max(best_p, .Machine$double.xmin))
  })
  list(x = sites, y = y_profile)
}