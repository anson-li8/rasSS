library(rasSS)
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)
r <- scan_block("data/sumstats/scz2018_clozuk.tsv.gz", "data/sumstats/scz2022.tsv.gz",
                6, 29e6, 34e6, "data/1000G/chr6_eur",
                plink2 = "data/1000G/plink2.exe",
                out_rds = "data/results/blk_scz_chr6_mhc.rds")
pk <- r$peaks
print(data.frame(snp = pk$snp_index, bp = pk$bp))
pass <- any(pk$bp >= 32e6 & pk$bp <= 33e6)
cat(if (pass) "POSITIVE CONTROL: PASS\n" else "POSITIVE CONTROL: FAIL\n")