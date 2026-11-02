library(rasSS); library(RAS); library(data.table)

targ_path <- "data/sumstats/scz2022.tsv.gz"        # PGC3 SCZ 2022, EUR autosome
disc_path <- "data/sumstats/scz2018_clozuk.tsv.gz" # PGC2+CLOZUK 2018
bim_path  <- "data/1000G/chr6_block.bim"
raw_path  <- "data/1000G/chr6_block_dosage.raw"
stopifnot(file.exists(targ_path), file.exists(disc_path),
          file.exists(bim_path), file.exists(raw_path))
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)

# 1. Ingest MHC block (chr6:29-34Mb)
targ <- read_sumstats(targ_path, chr = 6, start_bp = 29e6, end_bp = 34e6)
disc <- read_sumstats(disc_path, chr = 6, start_bp = 29e6, end_bp = 34e6)
cat("targ:", nrow(targ), "| disc:", nrow(disc), "\n")   # expect ~34434 / ~31141

# 2. Reference + two-sample harmonization
bim <- read_bim(bim_path)
ld  <- ld_matrix(raw_path)
h   <- harmonize_two(disc, targ, bim, ld)
cat("kept after harmonization:", h$n_kept, "of", h$n_in, "\n")  # expect 15774 / 18349
stopifnot(h$n_kept > 10000)

# 3. Scan + peak detection (parameters validated on the chr21 sweep)
scan <- scan_ss(h$w, h$Z, h$R, mask = rep(TRUE, length(h$w)),
                skip1 = 10, skip2 = 5, min_window_size = 5, max_window_size = 20)
pk <- call_peaks(scan, h$bp, window_size = 20, slope_check_window_size = 8,
                 slope_threshold = 1e-4, davies_threshold = 1e-5, skip = 10)

png("data/results/scz_mhc_control.png", width = 1000, height = 600)
plot(scan$x, scan$y, type = "l", xlab = "SNP index (chr6:29-34Mb)",
     ylab = "-log10(p)", main = "SCZ positive control: disc=PGC2018, targ=PGC2022")
if (length(pk$bp)) abline(v = pk$snp_index, col = "red", lty = 2)
dev.off()
res <- data.frame(snp = pk$snp_index, bp = pk$bp)
print(res)
saveRDS(list(scan = scan, peaks = pk, h_n_kept = h$n_kept),
        "data/results/scz_mhc_control.rds")

# 4. Pre-registered success criterion: peak inside HLA class II (32.0-33.0 Mb)
pass <- any(pk$bp >= 32e6 & pk$bp <= 33e6)
cat(if (pass) "POSITIVE CONTROL: PASS\n" else "POSITIVE CONTROL: FAIL\n")