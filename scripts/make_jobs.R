# Rscript scripts/make_jobs.R DISC TARG REFDIR OUT_CSV
library(data.table); library(rasSS)
a <- commandArgs(TRUE); stopifnot(length(a) == 4)
disc <- read_sumstats(a[1])[, .(snp, chr, bp)][snp != "."]
targ <- read_sumstats(a[2])[, .(snp, chr, bp)][snp != "."]

# Fix: Force chr to numeric and extract rsIDs directly here
disc$chr <- as.integer(disc$chr); targ$chr <- as.integer(targ$chr)
disc$snp <- sub("^.*_(rs[0-9]+)$", "\\1", disc$snp)
targ$snp <- sub("^.*_(rs[0-9]+)$", "\\1", targ$snp)

MAXN <- 20000; CORE <- 4e6; PAD <- 1e6
rows <- list()
for (ci in 1:22) {
  bim <- fread(file.path(a[3], sprintf("chr%d_eur.bim", ci)), header = FALSE, select = c(2, 4))
  setnames(bim, c("id", "bp"))
  # Find intersecting rsIDs across all three datasets
  ids <- Reduce(intersect, list(disc[chr == ci]$snp, targ[chr == ci]$snp, bim$id))
  # Pull coordinates EXCLUSIVELY from the 1000G reference panel (hg19)
  com <- bim[id %in% ids]                      
  cat(sprintf("chr%02d: intersect=%d\n", ci, nrow(com)))
  if (nrow(com) == 0) next
  emit <- function(cs, ce) {
    n <- com[between(bp, cs - PAD, ce + PAD), .N]
    if (n > MAXN) {
      if (ce - cs <= 2e6) { cat("  SKIP dense:", cs, ce, n, "\n"); return(invisible()) }
      mid <- cs + (ce - cs) / 2; emit(cs, mid); emit(mid, ce)
    } else if (n >= 500) {                     # matches scan_block's lower gate
      rows[[length(rows) + 1]] <<- data.table(chr = ci, core_start = cs, core_end = ce,
                                              pad_start = cs - PAD, pad_end = ce + PAD, n_snp = n)
    }
  }
  s <- 0; span <- max(bim$bp)
  while (s < span) { emit(s, min(s + CORE, span)); s <- s + CORE }
}
if (length(rows) == 0) {
  stop("Zero jobs returned: SNP-ID intersection empty on all chromosomes.")
}
jobs <- rbindlist(rows); setkey(jobs, chr, core_start)
fwrite(jobs, a[4]); cat("jobs:", nrow(jobs), "\n")

# results 725 blocks, skipping excessively dense regions due to compute restrictions
# skip MHC (chr 7 30Mb-34Mb)
# skip ch8 2Mb-6Mb (8p23.1)
# skip ch16 6Mb-8Mb (16p13.11/16p12.3 region)