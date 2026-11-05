# Rscript scripts/run_chr.R CHR JOBS_CSV DISC TARG REFDIR OUTDIR PLINK2
library(data.table); library(rasSS)
a <- commandArgs(TRUE); stopifnot(length(a) == 7)
ci   <- as.integer(a[1])
jobs <- fread(a[2])
jobs <- jobs[jobs$chr == ci]
stopifnot(nrow(jobs) > 0)
dir.create(a[6], showWarnings = FALSE, recursive = TRUE)

disc <- read_sumstats(a[3], chr = ci); targ <- read_sumstats(a[4], chr = ci)
disc$snp <- sub("^.*_(rs[0-9]+)$", "\\1", disc$snp)
targ$snp <- sub("^.*_(rs[0-9]+)$", "\\1", targ$snp)

bim  <- read_bim(file.path(a[5], sprintf("chr%d_eur.bim", ci)))
fail <- file.path(a[6], sprintf("fail_chr%d.tsv", ci))

# Because of the mismatch between iPSYCH's and SPARK's genome assemblies (hg19 vs hg38),
# rsID-based matching is used to ensure matching SNPs.
for (i in seq_len(nrow(jobs))) {
  j <- jobs[i]
  out <- file.path(a[6], sprintf("blk_chr%d_%d_%d.rds", ci, j$core_start, j$core_end))
  if (file.exists(out)) { cat("skip", out, "\n"); next }
  
  # FIX: Subset SNPs using the BIM coordinates (hg19) to avoid build mismatch.
  keep <- bim[between(bp, j$pad_start, j$pad_end)]$id
  keep <- keep[keep %in% disc$snp & keep %in% targ$snp]
  
  if (length(keep) < 500 || length(keep) > 20000) {
    fwrite(data.table(chr = ci, j$core_start, j$core_end, n = length(keep), why = "bounds"), fail, append = TRUE); next 
  }
  
  snpl <- tempfile(fileext = ".txt"); writeLines(keep, snpl); blk <- tempfile("block")
  
  r1 <- system2(a[7], c("--bfile", file.path(a[5], sprintf("chr%d_eur", ci)),
                        "--extract", snpl, "--make-bed", "--out", blk), stdout = TRUE, stderr = TRUE)
  r2 <- system2(a[7], c("--bfile", blk, "--export", "A", "--out", blk), stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(r1, "status")) || !is.null(attr(r2, "status"))) {
    fwrite(data.table(chr = ci, j$core_start, j$core_end, n = length(keep), why = "plink"), fail, append = TRUE); next 
  }
  
  ld <- ld_matrix(paste0(blk, ".raw")); bimb <- read_bim(paste0(blk, ".bim"))
  
  h <- harmonize_two(disc[snp %in% keep], targ[snp %in% keep], bimb, ld)
  
  if (h$n_kept / h$n_in < 0.6) {
    fwrite(data.table(chr = ci, j$core_start, j$core_end, n = h$n_kept, why = "retention"), fail, append = TRUE); next 
  }
  
  scan <- scan_ss(h$w, h$Z, h$R, mask = rep(TRUE, length(h$w)),
                  skip1 = 10, skip2 = 5, min_window_size = 5, max_window_size = 20)
  pk <- call_peaks(scan, h$bp)
  pkd <- data.table(snp = pk$snp_index, bp = pk$bp)[between(bp, j$core_start, j$core_end - 1)]
  
  saveRDS(list(chr = ci, core = c(j$core_start, j$core_end), n_kept = h$n_kept,
               n_in = h$n_in, scan = scan, peaks = pkd), out)
  cat("done", out, "| peaks:", nrow(pkd), "\n")
}