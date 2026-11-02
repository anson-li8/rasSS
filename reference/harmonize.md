# Harmonize sumstats against an LD reference (single study)

Orients Z and weights to the allele the dosage matrix counts. Uses the
exact counted allele from PLINK if available, falling back to the EAF
frequency match (XOR rule). Drops palindromic SNPs and
frequency-ambiguous SNPs.

## Usage

``` r
harmonize(sumstats, bim, ld, freq_tol = 0.15)
```

## Arguments

- sumstats:

  Output of
  [`read_sumstats`](https://anson-li8.github.io/rasSS/reference/read_sumstats.md).

- bim:

  Output of
  [`read_bim`](https://anson-li8.github.io/rasSS/reference/read_bim.md).

- ld:

  Output of
  [`ld_matrix`](https://anson-li8.github.io/rasSS/reference/ld_matrix.md).

- freq_tol:

  Maximum \|EAF - reference frequency\| tolerance.

## Value

List with w, Z, bp, R, ids, n_kept, all aligned to each other.
