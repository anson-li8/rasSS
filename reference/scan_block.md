# Run one genomic block end-to-end

Intersects discovery/target sumstats with the reference, extracts and
exports dosages via PLINK2, harmonizes, scans, and calls peaks.
Idempotent when `out_rds` exists.

## Usage

``` r
scan_block(
  disc_path,
  targ_path,
  chr,
  start_bp,
  end_bp,
  ref_prefix,
  plink2 = "plink2",
  out_rds = NULL,
  skip1 = 10,
  skip2 = 5,
  min_window_size = 5,
  max_window_size = 20,
  window_size = 20,
  slope_check_window_size = 8,
  slope_threshold = 1e-04,
  davies_threshold = 1e-05,
  max_snps = 25000
)
```

## Arguments

- disc_path:

  Path to discovery GWAS summary statistics.

- targ_path:

  Path to target GWAS summary statistics.

- chr:

  Chromosome number.

- start_bp:

  Block start position (base pairs).

- end_bp:

  Block end position (base pairs).

- ref_prefix:

  Prefix of the PLINK binary reference panel (bfile stem).

- plink2:

  Path to the plink2 executable.

- out_rds:

  Output `.rds` path; if it exists, the call is skipped.

- skip1:

  Number of leading SNPs skipped during the scan (discovery side).

- skip2:

  Number of leading SNPs skipped during the scan (target side).

- min_window_size:

  Minimum CPD window size in SNPs.

- max_window_size:

  Maximum CPD window size in SNPs.

- window_size:

  CPD window size used for peak calling.

- slope_check_window_size:

  CPD slope-check window size.

- slope_threshold:

  CPD slope p-value threshold.

- davies_threshold:

  CPD Davies p-value threshold.

- max_snps:

  Maximum allowed SNPs after intersection.

## Value

A list with `peaks`, `scan`, `bp`, `n_intersect`, `n_kept`, `n_in`.
