# Read GWAS summary statistics into a standardized table

Auto-detects common column schemas (CNCR/Savage, PGC, GCST). Derives Z
from beta/SE, or from P and effect sign when Z is absent. Handles files
with or without an Effect Allele Frequency (EAF) column.

## Usage

``` r
read_sumstats(path, chr = NULL, start_bp = NULL, end_bp = NULL)
```

## Arguments

- path:

  Path to (optionally gzipped) delimited sumstat file.

- chr:

  Optional chromosome to subset.

- start_bp, end_bp:

  Optional base-pair window.

## Value

A `data.table` with standardized columns: snp, chr, bp, a1, a2, z, beta,
eaf.
