library(rasSS); 
library(RAS); 
library(data.table)

sumstat_path <- "data/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt"
raw_path     <- "data/1000G/pilot_dosage.raw"
bim_path     <- "data/1000G/pilot_block.bim"
stopifnot(file.exists(sumstat_path), file.exists(raw_path), file.exists(bim_path))

# 1. Ingest target block (Savage 2018, chr21:15-20Mb)
ss <- read_sumstats(sumstat_path, chr = 21, start_bp = 15e6, end_bp = 20e6)
cat("sumstats in block:", nrow(ss), "\n")                 # expect 15781

# 2. Reference side
bim <- read_bim(bim_path)
ld  <- ld_matrix(raw_path)
cat("bim SNPs:", nrow(bim), "| dosage SNPs:", length(ld$ids), "\n")  # expect 12466 / 12466

# 3. Harmonize (frequency-based orientation, palindromic drop)
h <- harmonize(ss, bim, ld)
cat("kept after harmonization:", h$n_kept, "of", h$n_in, "\n")      # expect 10747 / 12466 (counted-allele harmonizer)
stopifnot(h$n_kept > 1000)

# 4. Summary-stat scan
scan <- scan_ss(h$w, h$Z, h$R, mask = rep(TRUE, length(h$w)),
                skip1 = 10, skip2 = 5, min_window_size = 5, max_window_size = 20)

# 5. Peak detection (sweep-validated: ws=20, scw=8, sl=1e-4)
pk <- call_peaks(scan, h$bp, window_size = 20, slope_check_window_size = 8,
                 slope_threshold = 1e-4, davies_threshold = 1e-5, skip = 10)

plot(scan$x, scan$y, type = "l", xlab = "SNP index (chr21:15-20Mb)",
     ylab = "-log10(p)", main = "rasSS pilot via ingest.R")
if (length(pk$bp)) abline(v = pk$snp_index, col = "red", lty = 2)
print(data.frame(grid = pk$grid_index, snp = pk$snp_index, bp = pk$bp))