# NOTE: Significant cohort overlap btw Sniekers 2017 and Savage 2018 -> winner's curse bias
# only used to demonstrate portability and continuous-trait support

library(rasSS)
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)

disc_path <- "data/sumstats/sniekers_2017_intelligence.txt"
targ_path <- "data/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt"

r <- scan_block(disc_path, targ_path, chr = 21, start_bp = 15e6, end_bp = 20e6,
                ref_prefix = "data/1000G/chr21_eur",
                plink2 = "plink2.exe",
                out_rds = "data/results/blk_int_chr21_15_20.rds")

targ <- read_sumstats(targ_path, chr = 21, start_bp = 15e6, end_bp = 20e6)
marg <- -log10(2 * pnorm(-abs(targ$z)))

plot(targ$bp / 1e6, marg, type = "p", pch = 16, cex = 0.25,
     col = rgb(0.55, 0.55, 0.55, 0.30),
     xlab = "Chr21 position (Mb, GRCh37)", ylab = expression(-log[10](p)),
     main = "Intelligence: disc = Sniekers 2017, targ = Savage 2018",
     ylim = c(0, max(marg, r$scan$y) * 1.05))
lines(r$bp[r$scan$x] / 1e6, r$scan$y, col = "darkblue", lwd = 1)
if (length(r$peaks$bp)) abline(v = r$peaks$bp / 1e6, col = "red", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("rasSS profile", "target marginal", "validated changepoint"),
       col = c("darkblue", rgb(0.55, 0.55, 0.55, 0.6), "red"),
       lty = c(1, NA, 2), pch = c(NA, 16, NA), lwd = c(1, NA, 1.5),
       bty = "o", bg = "white", cex = 0.8)

print(data.frame(snp = r$peaks$snp_index, bp = r$peaks$bp))
cat(sprintf("intersect %d | kept %d (%.1f%%) | peaks %d\n",
            r$n_intersect, r$n_kept, 100 * r$n_kept / r$n_in, length(r$peaks$bp)))