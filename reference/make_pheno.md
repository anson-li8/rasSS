# Simulate phenotype from causal effects

Generates phenotype based on linear combination of genotypes and true
effect, with added Gaussian noise.

## Usage

``` r
make_pheno(X, beta, trait, seed)
```

## Arguments

- X:

  Numeric matrix. Genotype dosage matrix (individuals x SNPs).

- beta:

  Numeric vector. True effect sizes (0 for non-causal SNPs).

- trait:

  Character. Either `"continuous"` or `"binary"`.

- seed:

  Integer. Random seed for reproducibility.

## Value

Phenotype vector. (numeric liability for continuous, 0/1 status for
binary).
