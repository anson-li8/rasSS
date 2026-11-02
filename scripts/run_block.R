args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) stop("usage: run_block.R DISC TARG CHR START END REF_PREFIX OUT_RDS PLINK2")
d   <- args[1]
t   <- args[2]
chr <- args[3]
s   <- args[4]
e   <- args[5]
ref <- args[6]
out <- args[7]
pl  <- args[8]
options(error = function() { traceback(2, max.lines = 20); quit(status = 1, save = "no") })
library(rasSS)
r <- scan_block(d, t, as.integer(chr), as.numeric(s), as.numeric(e), ref,
                plink2 = pl, out_rds = out)
cat("peaks bp:", paste(r$peaks$bp, collapse = ", "), "\n")