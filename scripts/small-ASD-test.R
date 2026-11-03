library(data.table); library(rasSS)
print(names(fread("data/sumstats/iPSYCH-PGC_ASD_Nov2017.gz", nrows = 0)))
print(names(fread("data/sumstats/SPARK_update_model1QC_EUR_only.tsv.gz", nrows = 0)))
r <- scan_block("data/sumstats/SPARK_update_model1QC_EUR_only.tsv.gz",
                "data/sumstats/iPSYCH-PGC_ASD_Nov2017.gz",
                21, 15e6, 20e6, "data/1000G/chr21_eur",
                plink2 = "plink2.exe",
                out_rds = "data/results/asd_mech_check.rds")
cat("intersect:", r$n_intersect, "| kept:", r$n_kept, "/", r$n_in,
    "| peaks:", length(r$peaks$bp), "\n")