# Marginal GWAS statistics for continuous trait

Calculates per-SNP marginal effect sizes and Z-scores via ordinary least
squares (OLS).

## Usage

``` r
get_marginal_stats(X, y)
```

## Arguments

- X:

  Numeric matrix. Genotypic dosage matrix (individuals x SNPs).

- y:

  Numeric vector. Continuous phenotype vector.

## Value

List with two elements:

- `beta`: Marginal effect sizes.

- `z`: Marginal Z-scores (non-finite values set to 0).
