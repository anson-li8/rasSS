library(data.table); library(rasSS)

scan_ss_offset <- function(w, Z, R, offset, skip1 = 10, skip2 = 5,
                           min_ws = 5, max_ws = 20) {
  n <- length(w)
  sites <- seq(1, n, by = skip1)
  sites <- sites[sites - max_ws >= 1 & sites + max_ws + offset <= n]
  sw <- c(0, seq(min_ws, max_ws, by = skip2))
  if (tail(sw, 1) != max_ws) sw <- c(sw, max_ws)
  y <- vapply(sites, function(s) {
    best <- 1
    for (ws in sw) {
      iw <- if (ws == 0) s else (s-ws):(s+ws)
      jw <- iw + offset # Z and R co-shifted!
      num <- sum(w[iw] * Z[jw])
      den <- sqrt(as.numeric(crossprod(w[iw], R[jw, jw, drop = FALSE] %*% w[iw])))
      if (!is.finite(den) || den <= 0) next
      best <- min(best, 2 * pnorm(-abs(num / den)))
    }
    -log10(max(best, .Machine$double.xmin))
  }, numeric(1))
  list(x = sites, y = y)
}

DISC <- "data/sumstats/SPARK_update_model1QC_EUR_only.tsv.gz"
TARG <- "data/sumstats/iPSYCH-PGC_ASD_Nov2017.gz"
REF  <- "data/1000G"; PLINK <- "plink2.exe"

jobs <- fread("data/jobs.csv")
qc <- readRDS("data/results/asd_peaks_qc.rds")

set.seed(7); sel <- jobs[sample(.N, 100)]
null_hs <- list()
null_peaks_per_shift <- numeric()
cat("Running corrected shift null on 100 blocks...\n")

for (i in 1:nrow(sel)) {
  j <- sel[i]; chr <- j$chr; start <- j$core_start; end <- j$core_end
  cat(sprintf("Block %d/100 (chr%d:%d-%d) ... ", i, chr, start, end))
  
  disc <- read_sumstats(DISC, chr = chr); targ <- read_sumstats(TARG, chr = chr)
  disc$snp <- sub("^.*_(rs[0-9]+)$", "\\1", disc$snp)
  targ$snp <- sub("^.*_(rs[0-9]+)$", "\\1", targ$snp)
  bim  <- read_bim(file.path(REF, sprintf("chr%d_eur.bim", chr)))
  keep <- bim[between(bp, start - 1e6, end + 1e6)]$id
  keep <- keep[keep %in% disc$snp & keep %in% targ$snp]
  if (length(keep) < 500 || length(keep) > 20000) { cat("skip\n"); next }
  
  snpl <- tempfile(fileext = ".txt"); writeLines(keep, snpl); blk <- tempfile("block")
  system2(PLINK, c("--bfile", file.path(REF, sprintf("chr%d_eur", chr)),
                   "--extract", snpl, "--make-bed", "--out", blk), stdout = FALSE)
  system2(PLINK, c("--bfile", blk, "--export", "A", "--out", blk), stdout = FALSE)
  
  ld <- ld_matrix(paste0(blk, ".raw")); bimb <- read_bim(paste0(blk, ".bim"))
  h  <- harmonize_two(disc[snp %in% keep], targ[snp %in% keep], bimb, ld)
  if (h$n_kept / h$n_in < 0.6) { cat("skip ret\n"); next }
  
  n <- length(h$w)
  offsets <- round(n * c(0.15, 0.30, 0.45, 0.60, 0.75, 0.90))
  
  block_hs <- numeric()
  for (off in offsets) {
    s_off <- scan_ss_offset(h$w, h$Z, h$R, off)
    p_off <- call_peaks(s_off, h$bp)
    
    p_off_dt <- data.table(snp = p_off$snp_index, bp = p_off$bp)
    p_off_dt <- p_off_dt[between(bp, start, end - 1)]
    
    if (nrow(p_off_dt) > 0) {
      block_hs <- c(block_hs, s_off$y[match(p_off_dt$snp, s_off$x)])
    }
  }
  null_hs[[i]] <- block_hs
  null_peaks_per_shift <- c(null_peaks_per_shift, length(block_hs))
  cat(sprintf("found %d null peaks\n", length(block_hs)))

  rm(ld, h, bimb, disc, targ, bim, keep); gc()
}

# Ground Truth 2: Null rate must be < 0.05 per block-shift
mean_null_rate <- mean(null_peaks_per_shift) / 6
cat("\nMean null peaks per block-shift:", round(mean_null_rate, 4), "\n")
if (mean_null_rate > 0.05) {
  cat("WARNING: Null is still too noisy! Co-shift fix might not be working.\n")
} else {
  cat("SUCCESS: Null rate is clean.\n")
}

# Calculate FDR Curve for Shift Null (B)
all_null_hs <- unlist(null_hs)
grid <- seq(2.5, 8, by = 0.25)
tabB <- rbindlist(lapply(grid, function(y) {
  # Denominator is now 100 blocks * 6 shifts = 600
  E <- 725 * (sum(all_null_hs >= y) / (100 * 6)) 
  O <- sum(qc$ras >= y)
  data.table(y, E_false = round(E, 2), O_obs = O, FDR = ifelse(O > 0, E/O, NA))
}))

cat("\n=== CORRECTED SHIFT NULL FDR ===\n")
print(tabB)

valid_B <- tabB[!is.na(FDR) & FDR <= 0.10]
ystar_B <- if (nrow(valid_B) > 0) min(valid_B$y) else Inf
cat("\nystar_B (Corrected Shift Null FDR <= 0.10) =", ystar_B, "\n")
cat("Final Threshold = max(4.5,", ystar_B, ") =", max(4.5, ystar_B), "\n")

# Save the results just in case
saveRDS(list(tabB = tabB, ystar_B = ystar_B, null_hs = all_null_hs), "data/results/shift_null_final.rds")