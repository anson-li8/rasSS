library(data.table)

f <- list.files("data/results/asd", pattern="\\.rds$", full.names=TRUE)
qc <- rbindlist(lapply(f, function(p) {
  r <- readRDS(p)
  y <- r$scan$y[match(r$peaks$snp, r$scan$x)]
  if (any(is.na(y))) y[is.na(y)] <- r$scan$y[vapply(r$peaks$snp[is.na(y)],
                                                    function(i) which.min(abs(r$scan$x - i)), 1L)]
  data.table(chr=r$chr, core_start=r$core[1], core_end=r$core[2],
             bp=r$peaks$bp, ras=y, n_kept=r$n_kept, n_in=r$n_in, n_pk=nrow(r$peaks))
}))
saveRDS(qc, "data/results/asd_peaks_qc.rds")
cat("Saved", nrow(qc), "peaks to data/results/asd_peaks_qc.rds\n")