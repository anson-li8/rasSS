# Marginal GWAS statistics for binary trait

Calculates per-SNP Rao score test statistics for case/control
phenotypes.

## Usage

``` r
get_marginal_stats_bin(X, y)
```

## Arguments

- X:

  Numeric matrix. Genotype dosage matrix (individuals x SNPs).

- y:

  Numeric vector. Binary phenotype vector (0/1).

## Value

List with two elements:

- `beta`: Score-based effect sizes (\\U/V\\).

- `z`: Rao score Z-scores (non-finite values set to 0).
