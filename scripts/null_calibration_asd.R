library(data.table); library(rasSS)
jobs <- fread("data/jobs.csv")
brk <- quantile(jobs$n_snp, probs = 0:4/4)
jobs[, stratum := cut(n_snp, brk, include.lowest = TRUE)]
saveRDS(brk, "data/results/asd_null_strata.rds")
set.seed(42)
sel <- jobs[, .SD[sample(.N, min(.N, 10))], by = stratum]   # 40 blocks
B <- 40                                                     # sims per block -> 1600 block-sims
out <- list(); sims_done <- 0L
for (i in seq_len(nrow(sel))) {
  j <- sel[i]; chr <- j$chr; start <- j$core_start; end <- j$core_end
  cat(sprintf("--- block %d/40 chr%d:%d-%d ---\n", i, chr, start, end))
  disc <- read_sumstats("data/sumstats/SPARK_update_model1QC_EUR_only.tsv.gz", chr = chr)
  targ <- read_sumstats("data/sumstats/iPSYCH-PGC_ASD_Nov2017.gz", chr = chr)
  disc$snp <- sub("^.*_(rs[0-9]+)$", "\\1", disc$snp)
  targ$snp <- sub("^.*_(rs[0-9]+)$", "\\1", targ$snp)
  bim  <- read_bim(sprintf("data/1000G/chr%d_eur.bim", chr))
  keep <- bim[between(bp, start - 1e6, end + 1e6)]$id
  keep <- keep[keep %in% disc$snp & keep %in% targ$snp]
  if (length(keep) < 500) next
  snpl <- tempfile(fileext = ".txt"); writeLines(keep, snpl); blk <- tempfile("block")
  system2("plink2.exe", c("--bfile", sprintf("data/1000G/chr%d_eur", chr),
                          "--extract", snpl, "--make-bed", "--out", blk), stdout = FALSE)
  system2("plink2.exe", c("--bfile", blk, "--export", "A", "--out", blk), stdout = FALSE)
  ld <- ld_matrix(paste0(blk, ".raw")); bimb <- read_bim(paste0(blk, ".bim"))
  h  <- harmonize_two(disc[snp %in% keep], targ[snp %in% keep], bimb, ld)
  Gk <- ld$G[, match(h$ids, ld$ids), drop = FALSE]
  Gc <- sweep(Gk, 2, colMeans(Gk), "-")                  # exact null draws, no chol
  heights <- list()
  for (b in seq_len(B)) {
    set.seed(1000 * i + b)
    Z0 <- as.numeric(crossprod(Gc, rnorm(nrow(Gc)))) / sqrt(nrow(Gc) - 1)
    s0 <- scan_ss(h$w, Z0, h$R, rep(TRUE, length(h$w)), 10, 5, 5, 20)
    p0 <- call_peaks(s0, h$bp)
    pd <- data.table(bp = p0$bp, snp = p0$snp_index)[between(bp, start, end - 1)]
    heights[[b]] <- if (nrow(pd)) s0$y[match(pd$snp, s0$x)] else numeric(0)
  }
  sims_done <- sims_done + B
  out[[i]] <- data.table(stratum = j$stratum, sim = rep(seq_len(B), lengths(heights)),
                         height = unlist(heights))
  cat(sprintf(">>> null peaks: %d in %d sims\n", sum(lengths(heights)), B))
  rm(Gk, Gc, h, ld); gc()
}
nulls <- rbindlist(out)
saveRDS(list(peaks = nulls, sims = sims_done, B = B), "data/results/asd_null_heights.rds")
cat("DONE:", sims_done, "block-sims\n")