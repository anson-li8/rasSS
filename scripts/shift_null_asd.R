library(data.table); library(rasSS)

DISC  <- "data/sumstats/SPARK_update_model1QC_EUR_only.tsv.gz"
TARG  <- "data/sumstats/iPSYCH-PGC_ASD_Nov2017.gz"
REF   <- "data/1000G"; PLINK <- "plink2.exe"
OUTD  <- "data/results/nullcal"
FRAC  <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # circular-shift fractions
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

jobs <- fread("data/jobs.csv")
brk  <- quantile(jobs$n_snp, probs = 0:4/4)
jobs[, stratum := cut(n_snp, brk, include.lowest = TRUE)]
saveRDS(brk, "data/results/asd_null_strata.rds")
set.seed(42)
sel <- jobs[, .SD[sample(.N, min(.N, 10))], by = stratum] # 40 blocks

done_keys <- sub("^nc_(.*)\\.rds$", "\\1", list.files(OUTD, pattern = "^nc_.*\\.rds$"))
t0 <- proc.time()[3]

for (i in seq_len(nrow(sel))) {
  j <- sel[i]; chr <- j$chr; start <- j$core_start; end <- j$core_end
  key <- sprintf("%d_%d_%d", chr, start, end)
  if (key %in% done_keys) { cat(sprintf("[%2d/%d] chr%d:%d-%d skip (done)\n", i, nrow(sel), chr, start, end)); next }
  cat(sprintf("[%2d/%d] chr%d:%d-%d | harmonizing... ", i, nrow(sel), chr, start, end))
  disc <- read_sumstats(DISC, chr = chr); targ <- read_sumstats(TARG, chr = chr)
  disc$snp <- sub("^.*_(rs[0-9]+)$", "\\1", disc$snp)
  targ$snp <- sub("^.*_(rs[0-9]+)$", "\\1", targ$snp)
  bim  <- read_bim(file.path(REF, sprintf("chr%d_eur.bim", chr)))
  keep <- bim[between(bp, start - 1e6, end + 1e6)]$id
  keep <- keep[keep %in% disc$snp & keep %in% targ$snp]
  if (length(keep) < 500 || length(keep) > 20000) {
    cat("skip: keep =", length(keep), "\n")
    saveRDS(list(stratum = j$stratum, chr = chr, core_start = start, core_end = end,
                 n_shifts = 0L, skipped = TRUE), file.path(OUTD, paste0("nc_", key, ".rds")))
    next
  }
  snpl <- tempfile(fileext = ".txt"); writeLines(keep, snpl); blk <- tempfile("block")
  system2(PLINK, c("--bfile", file.path(REF, sprintf("chr%d_eur", chr)),
                   "--extract", snpl, "--make-bed", "--out", blk), stdout = FALSE, stderr = FALSE)
  system2(PLINK, c("--bfile", blk, "--export", "A", "--out", blk), stdout = FALSE, stderr = FALSE)
  ld <- ld_matrix(paste0(blk, ".raw")); bimb <- read_bim(paste0(blk, ".bim"))
  h  <- harmonize_two(disc[snp %in% keep], targ[snp %in% keep], bimb, ld)
  unlink(c(snpl, paste0(blk, c(".bed", ".bim", ".fam", ".raw"))))
  if (h$n_kept / h$n_in < 0.6) {
    cat("skip: retention", round(h$n_kept / h$n_in, 2), "\n")
    saveRDS(list(stratum = j$stratum, chr = chr, core_start = start, core_end = end,
                 n_shifts = 0L, skipped = TRUE), file.path(OUTD, paste0("nc_", key, ".rds")))
    rm(ld, h); gc(); next
  }
  n <- length(h$w)
  chi2      <- mean(h$Z^2) # background magnitude
  conc_obs  <- mean(h$w * h$Z) # real cross-cohort alignment
  conc_null <- mean(vapply(FRAC, function(f) { # chance alignment under shift
    k <- round(f * n); mean(h$w * c(h$Z[(k+1):n], h$Z[1:k])) }, numeric(1)))
  cat("shifts: ")
  null_hs <- list()
  for (s in seq_along(FRAC)) {
    k  <- round(FRAC[s] * n)
    Zs <- c(h$Z[(k+1):n], h$Z[1:k]) # circular shift: background kept, alignment destroyed
    s0 <- scan_ss(h$w, Zs, h$R, rep(TRUE, n), 10, 5, 5, 20)
    p0 <- call_peaks(s0, h$bp)
    pd <- data.table(bp = p0$bp, snp = p0$snp_index)[between(bp, start, end - 1)]
    null_hs[[s]] <- if (nrow(pd)) s0$y[match(pd$snp, s0$x)] else numeric(0)
    cat(if (length(null_hs[[s]])) "+" else ".")
  }
  cat(sprintf(" | kept %d/%d | chi2 %.2f | conc %.2f vs %.2f | null peaks %d/%d\n",
              h$n_kept, h$n_in, chi2, conc_obs, conc_null, sum(lengths(null_hs)), length(FRAC)))
  saveRDS(list(stratum = j$stratum, chr = chr, core_start = start, core_end = end,
               n_kept = h$n_kept, n_in = h$n_in, n_shifts = length(FRAC), skipped = FALSE,
               chi2 = chi2, conc_obs = conc_obs, conc_null = conc_null,
               nulls = data.table(shift = rep(seq_along(FRAC), lengths(null_hs)),
                                  height = unlist(null_hs))),
          file.path(OUTD, paste0("nc_", key, ".rds")))
  rm(ld, h, bimb, disc, targ, bim, keep); gc()
  el <- proc.time()[3] - t0
  cat(sprintf("    elapsed %.0f min | ETA %.0f min\n", el/60, el/60 * (nrow(sel) - i)/i))
}
cat("DONE.\n")