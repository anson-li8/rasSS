# Rscript scripts/verify_scan_block.R   (from C:/Users/Anson/Documents/Coding/rasSS)
library(rasSS)
ok_all <- TRUE
check <- function(label, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  if (!isTRUE(cond)) ok_all <<- FALSE
}
plink2 <- if (Sys.info()[["sysname"]] == "Windows") "data/1000G/plink2.exe" else "plink2"
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)

# --- A: SCZ MHC positive control, fresh run (re-anchors the baseline) ---
f <- "data/results/verify_scz_mhc.rds"; if (file.exists(f)) file.remove(f)
r <- scan_block("data/sumstats/scz2018_clozuk.tsv.gz", "data/sumstats/scz2022.tsv.gz",
                6, 29e6, 34e6, "data/1000G/chr6_eur", plink2 = plink2, out_rds = f)
check("A1 intersect == 21816",            r$n_intersect == 21816)
check("A2 kept == 18720 of 21816",        r$n_kept == 18720 && r$n_in == 21816)
check("A3 MHC class II peak 32.0-33.0Mb", any(r$peaks$bp >= 32e6 & r$peaks$bp <= 33e6))
check("A4 all peaks inside block",        all(r$peaks$bp >= 29e6 & r$peaks$bp <= 34e6))

# --- A5: idempotency (second call must skip, not recompute) ---
out2 <- capture.output(r2 <- scan_block(
  "data/sumstats/scz2018_clozuk.tsv.gz", "data/sumstats/scz2022.tsv.gz",
  6, 29e6, 34e6, "data/1000G/chr6_eur", plink2 = plink2, out_rds = f))
check("A5 idempotent skip + same peaks", any(grepl("skip: exists", out2)) &&
        identical(r2$peaks$bp, r$peaks$bp))

# --- B: chr21 mechanics parity vs the validated pilot fingerprint ---
f21 <- "data/results/verify_chr21_pilot.rds"; if (file.exists(f21)) file.remove(f21)
sv  <- "data/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt"
p <- scan_block(sv, sv, 21, 15e6, 20e6, "data/1000G/chr21_eur",
                plink2 = plink2, out_rds = f21)
check("B1 intersect == 12466", p$n_intersect == 12466)
check("B2 kept == 10747",      p$n_kept == 10747)
check("B3 pilot peaks reproduced exactly",
      identical(as.integer(p$peaks$snp_index), c(981L, 5901L, 8511L, 8791L, 9461L)))

# --- C: CLI parity and failure semantics (HPC contract) ---
fc <- "data/results/verify_cli.rds"; if (file.exists(fc)) file.remove(fc)
rc <- system2("Rscript", c("scripts/run_block.R",
                           "data/sumstats/scz2018_clozuk.tsv.gz", "data/sumstats/scz2022.tsv.gz",
                           "6", "29000000", "34000000", "data/1000G/chr6_eur", fc, plink2),
              stdout = TRUE, stderr = TRUE)
cat("---- CLI stdout/stderr ----\n")
cat(rc, sep = "\n")
cat("\n---- end ----\n")
cat("Rscript:", Sys.which("Rscript"), "\nwd:", getwd(), "\n")
check("C1 CLI exits 0", is.null(attr(rc, "status")))
check("C2 CLI peaks == package peaks", identical(readRDS(fc)$peaks$bp, r$peaks$bp))
bad <- system2("Rscript", c("scripts/run_block.R",
                            "nope.tsv", "nope2.tsv",
                            "6", "1", "2",
                            "data/1000G/chr6_eur",
                            plink2,
                            "data/results/verify_bad.rds"),
               stdout = TRUE, stderr = TRUE)
cat(bad, sep = "\n")
check("C3 bad input exits non-zero", !is.null(attr(bad, "status")))

cat(if (ok_all) "\nALL CHECKS PASS\n" else "\nSOME CHECKS FAILED\n")
if (!ok_all) quit(status = 1)